function R = t01_config_and_points()
%T01  Configuration resolution, grid construction, IDs and the point list.
%   No physics, no data files are read for their contents -- this is the test
%   that the plumbing which replaced the copied drivers actually works.

    R = struct('pass',true,'skip',false,'msg',{{}});
    function chk(cond, fmt, varargin)
        if ~cond
            R.pass = false;
            R.msg{end+1} = sprintf(fmt, varargin{:});
            fprintf('    FAIL  %s\n', R.msg{end});
        end
    end

    CASE  = 'unheated_c1_NoSBLI_Periodic';
    STUDY = 'regression_unheated_c1';

    % ---------- 1. config resolution ------------------------------------
    C = qs_config(CASE, STUDY);
    chk(strcmp(C.case.case_id, CASE),    'case_id wrong: %s', C.case.case_id);
    chk(strcmp(C.study.study_id, STUDY), 'study_id wrong: %s', C.study.study_id);
    chk(abs(C.case.pinf_Pa - 50287) < 1e-9,  'pinf not from case.json: %g', C.case.pinf_Pa);
    chk(abs(C.aero.gamma  - 1.4)    < 1e-12, 'gamma not from defaults: %g', C.aero.gamma);
    chk(abs(C.duration.base_t_end - 3.0) < 1e-12, 'study duration did not override defaults: %g', C.duration.base_t_end);
    chk(abs(C.solver.RelTol - 1e-7) < 1e-20, 'solver RelTol wrong: %g', C.solver.RelTol);
    chk(exist(C.paths.rom_file,'file')==2,  'rom_file does not resolve: %s',  C.paths.rom_file);
    chk(exist(C.paths.aero_file,'file')==2, 'aero_file does not resolve: %s', C.paths.aero_file);
    chk(exist(C.paths.ref_file,'file')==2,  'ref_file does not resolve: %s',  C.paths.ref_file);

    % override precedence: caller beats study beats case beats defaults
    C2 = qs_config(CASE, STUDY, struct('duration', struct('base_t_end', 7.5)));
    chk(abs(C2.duration.base_t_end - 7.5) < 1e-12, 'overrides did not win: %g', C2.duration.base_t_end);
    chk(abs(C2.duration.dt_out - 1e-4) < 1e-20, 'override clobbered a sibling field');

    % ---------- 2. the grid must reproduce the original driver exactly ----
    % These two vectors are copied verbatim from the ORIGINAL
    % run_detailed_stability_map.m for case 1. If the study file and qs_axis
    % do not reproduce them, every regression comparison below is meaningless.
    pc_old = unique(round([48.5:0.5:49.0, 49.5:0.1:53.5, 54.0:0.5:55.0]*1e3)/1e3);
    dT_old = unique(round([6.0:0.75:8.25, 9.0:0.5:11.0, 11.5:0.15:15.5, 16.0:0.5:18.0]*1e3)/1e3);

    pc_new = qs_axis(C.study.grid.pc_kPa, 'pc');
    dT_new = qs_axis(C.study.grid.dT_K,   'dT');

    chk(numel(pc_new)==numel(pc_old), 'p_c axis has %d values, original had %d', numel(pc_new), numel(pc_old));
    chk(numel(dT_new)==numel(dT_old), 'dT axis has %d values, original had %d',  numel(dT_new), numel(dT_old));
    if numel(pc_new)==numel(pc_old)
        d = max(abs(pc_new(:)-pc_old(:)));
        chk(d < 1e-9, 'p_c axis differs from the original by up to %g kPa', d);
    end
    if numel(dT_new)==numel(dT_old)
        d = max(abs(dT_new(:)-dT_old(:)));
        chk(d < 1e-9, 'dT axis differs from the original by up to %g K', d);
    end
    fprintf('    grid: %d x %d = %d points (original: %d)\n', ...
        numel(pc_new), numel(dT_new), numel(pc_new)*numel(dT_new), numel(pc_old)*numel(dT_old));

    % ---------- 3. point ordering matches ind2sub([npc ndT], idx) --------
    P = qs_make_points(C);
    npc = numel(pc_new);
    for probe = [1, 2, npc, npc+1, numel(P.id)]
        [ip, jd] = ind2sub([npc, numel(dT_new)], probe);
        chk(abs(P.pc_Pa(probe) - pc_new(ip)*1e3) < 1e-6, ...
            'point %d p_c is %g, expected %g (ordering changed)', probe, P.pc_Pa(probe), pc_new(ip)*1e3);
        chk(abs(P.dT_K(probe) - dT_new(jd)) < 1e-9, ...
            'point %d dT is %g, expected %g (ordering changed)', probe, P.dT_K(probe), dT_new(jd));
    end

    % ---------- 4. IDs ---------------------------------------------------
    chk(numel(unique(P.id))==numel(P.id), 'point IDs are not unique: %d ids, %d unique', ...
        numel(P.id), numel(unique(P.id)));
    id1 = qs_id(51927, 12.76, 3.0, 'flat');
    id2 = qs_id(51927, 12.76, 3.0, 'flat');
    id3 = qs_id(51927, 12.76, 10.0, 'flat');
    id4 = qs_id(51927, 12.76, 3.0, 'flat_neg');
    chk(strcmp(id1,id2), 'qs_id is not deterministic');
    chk(~strcmp(id1,id3), 'qs_id ignores t_end: %s vs %s', id1, id3);
    chk(~strcmp(id1,id4), 'qs_id ignores ic_tag: %s vs %s', id1, id4);
    chk(isempty(regexp(id1,'[^A-Za-z0-9_]','once')), 'qs_id produced an unsafe filename: %s', id1);
    fprintf('    id example: %s\n', id1);

    % ---------- 5. slicing is a complete, disjoint partition -------------
    for n = [1 3 7 48]
        seen = {};
        for k = 1:n
            S = qs_pointlist('slice', P, k, n);
            seen = [seen; S.id(:)]; %#ok<AGROW>
        end
        chk(numel(seen)==numel(P.id), 'slicing into %d loses/duplicates points: %d vs %d', ...
            n, numel(seen), numel(P.id));
        chk(numel(unique(seen))==numel(P.id), 'slicing into %d produces duplicates', n);
    end

    % ---------- 6. append is idempotent ----------------------------------
    A = qs_pointlist('rows', P, (1:10)');
    B = qs_pointlist('append', A, A);
    chk(numel(B.id)==10, 'appending a list to itself added %d rows', numel(B.id)-10);

    % a longer re-run of the same operating point is a NEW row, by design
    Q = qs_pointlist('make', A.pc_Pa(1), A.dT_K(1), struct('t_end',10.0,'t_trans',2.0));
    B2 = qs_pointlist('append', A, Q);
    chk(numel(B2.id)==11, 'a longer re-run of an existing point was swallowed (%d rows)', numel(B2.id));

    % ---------- 7. guard against the t_trans >= t_end mistake -------------
    threw = false;
    try
        qs_pointlist('make', 50000, 12, struct('t_end',1.0,'t_trans',1.5));
    catch
        threw = true;
    end
    chk(threw, 't_trans >= t_end was accepted; it should be rejected');

    if R.pass, fprintf('    all config / grid / point-list checks passed\n'); end
end
