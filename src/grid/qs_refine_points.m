function P = qs_refine_points(caseId, studyId, dry_run)
%QS_REFINE_POINTS  Add points where the regime map is undecided, and add time
%                  where the label is not yet trustworthy.
%
%   P = QS_REFINE_POINTS(caseId, studyId)            append and save
%   P = QS_REFINE_POINTS(caseId, studyId, true)      report only, change nothing
%
%   Run AFTER qs_merge. Two independent things happen, and they answer two
%   different questions:
%
%   (1) REFINE IN SPACE -- "where exactly is the boundary?"
%       Any two adjacent filled cells with different regime labels get a new
%       point at their midpoint, as long as the gap is still wider than
%       refine.min_dpc_kPa / refine.min_ddT_K. Repeated passes bisect the
%       boundary until it is resolved to that floor.
%
%   (2) REFINE IN TIME -- "is this label real?"
%       Any point flagged transition / nonstationary / near_threshold is
%       re-queued at the promoted run length. Because run length is part of
%       the point ID, the longer run is a NEW row: the short answer is kept
%       and the two can be compared. This is the mechanism that catches a
%       regime that flips at 8 s in a 3 s sweep.
%
%   New points are appended to the study's point list and the manifest is
%   updated. Resubmit the same fan-out afterwards: slices skip every ID that
%   already has a result, so only the new points run.

    if nargin < 3, dry_run = false; end
    C = qs_config(caseId, studyId);

    mf = fullfile(C.paths.merged_dir,'points.mat');
    assert(exist(mf,'file')==2, 'qs_refine_points: run qs_merge(''%s'',''%s'') first', caseId, studyId);
    S = load(mf);  D = S.D;
    P = qs_pointlist('load', C.paths.pointlist);

    R = getdef(C,'refine', struct());
    min_dpc = getdef(R,'min_dpc_kPa', 0.05);
    min_ddT = getdef(R,'min_ddT_K',   0.05);
    max_new = getdef(R,'max_new_points', 2000);

    nextLevel = max(P.level) + 1;
    G = qs_rasterize(D);

    % ================= (1) spatial bisection of label boundaries ===========
    newpc = [];  newdT = [];
    npc = numel(G.pc_kPa);  ndT = numel(G.dT_K);

    for j = 1:ndT                                  % along p_c
        fi = find(G.filled(:,j));
        for a = 1:numel(fi)-1
            i1 = fi(a);  i2 = fi(a+1);
            if G.label(i1,j) == G.label(i2,j), continue; end
            gap = G.pc_kPa(i2) - G.pc_kPa(i1);
            if gap <= min_dpc, continue; end
            newpc(end+1,1) = 0.5*(G.pc_kPa(i1)+G.pc_kPa(i2)); %#ok<AGROW>
            newdT(end+1,1) = G.dT_K(j);                        %#ok<AGROW>
        end
    end
    n_pc_side = numel(newpc);

    for i = 1:npc                                  % along dT
        fj = find(G.filled(i,:));
        for a = 1:numel(fj)-1
            j1 = fj(a);  j2 = fj(a+1);
            if G.label(i,j1) == G.label(i,j2), continue; end
            gap = G.dT_K(j2) - G.dT_K(j1);
            if gap <= min_ddT, continue; end
            newpc(end+1,1) = G.pc_kPa(i);                       %#ok<AGROW>
            newdT(end+1,1) = 0.5*(G.dT_K(j1)+G.dT_K(j2));       %#ok<AGROW>
        end
    end

    n_dT_side = numel(newpc) - n_pc_side;

    % de-duplicate the midpoints themselves (a cell can be a boundary in both
    % directions) and round onto the ID resolution
    if ~isempty(newpc)
        key = [round(newpc*1e3), round(newdT*1e3)];
        [~, ia] = unique(key, 'rows');
        newpc = newpc(ia);  newdT = newdT(ia);
    end

    % new boundary points get the promoted run length straight away: they are
    % by definition the points whose label is least certain
    DP = C.duration_policy;
    Pnew_space = qs_pointlist('make', newpc*1e3, newdT, struct( ...
        't_end',   getdef(DP,'promoted_t_end', 10.0), ...
        't_trans', getdef(DP,'promoted_t_transient', 2.0), ...
        'level',   nextLevel, ...
        'ic_tag',  getdef(getdef(C,'ic',struct()),'mode','flat'), ...
        'origin',  'refine'));

    % ================= (2) temporal promotion of shaky labels ==============
    promote = false(numel(D.id),1);
    want = getdef(DP,'promote_if', {'boundary','nonstationary','near_threshold'});
    if ischar(want), want = {want}; end
    for k = 1:numel(want)
        switch lower(want{k})
            case 'nonstationary',  promote = promote | D.nonstationary;
            case 'near_threshold', promote = promote | D.near_threshold;
            case 'transition',     promote = promote | D.trans_flag;
            case 'boundary'        % handled spatially above
            case 'divergent',      promote = promote | D.divergent;
            otherwise
                warning('qs_refine_points:badPolicy','unknown promote_if entry "%s"', want{k});
        end
    end
    promote = promote & D.ok;

    % a point that already transitioned inside a long run has nothing more to
    % learn from another run of the same length -- escalate it instead
    t_prom = getdef(DP,'promoted_t_end', 10.0);
    t_esc  = getdef(DP,'escalate_t_end', 20.0);
    tr_prom = getdef(DP,'promoted_t_transient', 2.0);
    tr_esc  = getdef(DP,'escalate_t_transient', 3.0);

    idx = find(promote);
    tnew  = zeros(numel(idx),1);
    trnew = zeros(numel(idx),1);
    for a = 1:numel(idx)
        k = idx(a);
        if D.t_end(k) >= t_prom - 1e-9
            tnew(a) = t_esc;   trnew(a) = tr_esc;
        else
            tnew(a) = t_prom;  trnew(a) = tr_prom;
        end
    end
    keep = tnew > (D.t_end(idx) + 1e-9);     % never re-queue at the same length
    idx = idx(keep);  tnew = tnew(keep);  trnew = trnew(keep);

    Pnew_time = qs_pointlist('empty');
    if ~isempty(idx)
        % Fields are assigned one at a time rather than through struct(...):
        % struct() expands a cell-array value (ic_tag here) into a STRUCT
        % ARRAY, which is not what is wanted and fails far from the cause.
        opts = struct();
        opts.t_end     = tnew;
        opts.t_trans   = trnew;
        opts.level     = nextLevel;
        opts.ic_tag    = D.ic_tag(idx);
        opts.origin    = 'refine';
        opts.save_hist = true;         % a promoted run is exactly the one to keep
        Pnew_time = qs_pointlist('make', D.pc_Pa(idx), D.dT_K(idx), opts);
    end

    % ================= append ==============================================
    before = numel(P.id);
    Q = qs_pointlist('append', qs_pointlist('empty'), Pnew_space);
    Q = qs_pointlist('append', Q, Pnew_time);
    Q = drop_existing(Q, P);
    nnew = numel(Q.id);

    fprintf('\n--- qs_refine_points  %s / %s  (level %d) ---\n', caseId, studyId, nextLevel);
    fprintf('  boundary midpoints  : %d  (%d along p_c, %d along dT, before dedup)\n', ...
        numel(Pnew_space.id), n_pc_side, n_dT_side);
    fprintf('  promoted for time   : %d  (longer re-runs of flagged points)\n', numel(Pnew_time.id));
    fprintf('  new after dedup     : %d\n', nnew);
    if nnew > 0
        fprintf('  estimated cost      : %.1f core-hours at 2.3 min per simulated second\n', ...
            sum(Q.t_end)*2.3/60);
    end

    if nnew > max_new
        error(['qs_refine_points: %d new points exceeds refine.max_new_points (%d).\n' ...
               '  Either raise the limit deliberately or coarsen the refinement floor.'], nnew, max_new);
    end

    if dry_run
        fprintf('  DRY RUN -- point list not modified.\n\n');
        return
    end
    if nnew == 0
        fprintf('  nothing to add: the boundary is resolved to the floor and no labels are flagged.\n\n');
        return
    end

    P = qs_pointlist('append', P, Q);
    qs_pointlist('save',  P, C.paths.pointlist);
    qs_pointlist('tocsv', P, fullfile(C.paths.results_dir,'pointlist.csv'));

    man = qs_jsonread(C.paths.manifest);
    man.n_points = numel(P.id);
    if ~isfield(man,'refine_history') || isempty(man.refine_history)
        man.refine_history = {};
    elseif ~iscell(man.refine_history)
        man.refine_history = num2cell(man.refine_history(:).');
    end
    man.refine_history{end+1} = struct('level',nextLevel, 'added',nnew, ...
        'total',numel(P.id), 'at',C.meta.resolved_at);
    qs_jsonwrite(C.paths.manifest, man);

    fprintf('  point list %d -> %d. Resubmit the fan-out; finished points are skipped.\n\n', ...
        before, numel(P.id));
end

% ======================================================================
function Q = drop_existing(Q, P)
    if isempty(Q.id), return; end
    keep = ~qs_pointlist('has', P, Q.id);
    Q = qs_pointlist('rows', Q, find(keep));
end

function v = getdef(S, f, d)
    if isstruct(S) && isfield(S,f) && ~isempty(S.(f)), v = S.(f); else, v = d; end
end
