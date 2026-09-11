function R = t04_sweep_roundtrip()
%T04  The whole pipeline end to end, on a deliberately tiny study.
%
%   init -> two slices of a sweep -> resume -> merge -> rasterize -> refine.
%   The physics is meaningless (0.2 s runs); what is being tested is that the
%   pieces that replaced the copied-folder workflow actually fit together:
%     * two slices cover the grid exactly once between them
%     * a re-run of a finished slice does no work and loses nothing
%     * merge recovers every point and writes a readable CSV
%     * rasterize puts them back on a grid in the right cells
%     * refine appends without renumbering and is idempotent when dry-run
%
%   Takes about 2 minutes.

    R = struct('pass',true,'skip',false,'msg',{{}});
    function chk(cond, fmt, varargin)
        if ~cond
            R.pass = false;
            R.msg{end+1} = sprintf(fmt, varargin{:});
            fprintf('    FAIL  %s\n', R.msg{end});
        end
    end

    CASE  = 'unheated_c1_NoSBLI_Periodic';
    STUDY = 'test_tiny_v1';

    % ---------- clean slate ------------------------------------------------
    C = qs_config(CASE, STUDY);
    if exist(C.paths.points_dir,'dir')==7
        old = dir(fullfile(C.paths.points_dir,'partial_*.mat'));
        for k = 1:numel(old), delete(fullfile(old(k).folder, old(k).name)); end
    end

    % ---------- init --------------------------------------------------------
    qs_init_study(CASE, STUDY, true);
    P = qs_pointlist('load', C.paths.pointlist);
    n0 = numel(P.id);
    chk(n0 == 4, 'tiny study should have 4 points, got %d', n0);
    chk(exist(C.paths.manifest,'file')==2, 'manifest was not written');
    chk(exist(fullfile(C.paths.results_dir,'pointlist.csv'),'file')==2, 'pointlist.csv was not written');

    man = qs_jsonread(C.paths.manifest);
    chk(isfield(man,'config') && isfield(man.config,'duration'), ...
        'the manifest does not carry the resolved config -- provenance is lost');

    % ---------- two slices ---------------------------------------------------
    fprintf('    running slice 1 of 2 ...\n');  qs_run_sweep(CASE, STUDY, 1, 2);
    fprintf('    running slice 2 of 2 ...\n');  qs_run_sweep(CASE, STUDY, 2, 2);

    p1 = load(fullfile(C.paths.points_dir,'partial_1_of_2.mat'));
    p2 = load(fullfile(C.paths.points_dir,'partial_2_of_2.mat'));
    ids1 = cellfun(@(r) r.id, p1.Rs, 'UniformOutput', false);
    ids2 = cellfun(@(r) r.id, p2.Rs, 'UniformOutput', false);
    chk(isempty(intersect(ids1, ids2)), 'the two slices ran overlapping points');
    chk(numel(union(ids1, ids2)) == n0, 'the two slices covered %d of %d points', ...
        numel(union(ids1,ids2)), n0);

    % ---------- resume does nothing ------------------------------------------
    d1 = dir(fullfile(C.paths.points_dir,'partial_1_of_2.mat'));
    fprintf('    re-running slice 1 (should skip everything) ...\n');
    t0 = tic;  qs_run_sweep(CASE, STUDY, 1, 2);  el = toc(t0);
    p1b = load(fullfile(C.paths.points_dir,'partial_1_of_2.mat'));
    chk(numel(p1b.Rs) == numel(p1.Rs), 'resume changed the partial: %d records became %d', ...
        numel(p1.Rs), numel(p1b.Rs));
    chk(el < 60, 'resume took %.0f s -- finished points were not skipped', el);
    fprintf('    resume skipped all %d points in %.1f s\n', numel(p1.Rs), el);

    % ---------- merge ----------------------------------------------------------
    D = qs_merge(CASE, STUDY);
    chk(numel(D.id) == n0, 'merge recovered %d of %d points', numel(D.id), n0);
    chk(D.meta.n_missing == 0, 'merge reports %d missing points', D.meta.n_missing);
    chk(all(D.ok), '%d point(s) failed to run', sum(~D.ok));
    chk(exist(fullfile(C.paths.merged_dir,'points.csv'),'file')==2, 'points.csv was not written');
    chk(all(isfinite(D.Amp_wh(D.ok))), 'some merged amplitudes are not finite');
    chk(all(cellfun(@(w) ~isempty(w), D.W)), 'some merged points carry no window signature');

    % ---------- rasterize ------------------------------------------------------
    G = qs_rasterize(D);
    chk(numel(G.pc_kPa)==2 && numel(G.dT_K)==2, 'raster is %dx%d, expected 2x2', ...
        numel(G.pc_kPa), numel(G.dT_K));
    chk(all(G.filled(:)), 'raster has %d empty cells for a complete 2x2 study', sum(~G.filled(:)));
    for k = 1:numel(D.id)
        i = find(abs(G.pc_kPa - D.pc_Pa(k)/1e3) < 1e-6);
        j = find(abs(G.dT_K   - D.dT_K(k))      < 1e-9);
        chk(~isempty(i) && ~isempty(j), 'point %s has no raster cell', D.id{k});
        if ~isempty(i) && ~isempty(j)
            chk(abs(G.Amp_wh(i,j) - D.Amp_wh(k)) < 1e-12, ...
                'point %s landed in the wrong raster cell', D.id{k});
        end
    end

    % ---------- refine ----------------------------------------------------------
    Pd = qs_refine_points(CASE, STUDY, true);      % dry run
    chk(numel(Pd.id) == n0, 'a dry-run refine modified the point list (%d -> %d)', n0, numel(Pd.id));

    P2 = qs_refine_points(CASE, STUDY, false);
    chk(numel(P2.id) >= n0, 'refine shrank the point list');
    chk(all(strcmp(P2.id(1:n0), P.id(1:n0))), 'refine renumbered existing points -- old results would be orphaned');
    fprintf('    refine: %d -> %d points, originals untouched\n', n0, numel(P2.id));

    % refining again with no new results must add nothing
    n2 = numel(P2.id);
    P3 = qs_refine_points(CASE, STUDY, false);
    chk(numel(P3.id) == n2, 'a second refine with no new results added %d points', numel(P3.id)-n2);
end
