function C = qs_init_study(caseId, studyId, force)
%QS_INIT_STUDY  Create the results folder, point list and manifest for a study.
%
%   qs_init_study(caseId, studyId)          create if absent, otherwise report
%   qs_init_study(caseId, studyId, true)    rebuild the point list from scratch
%
%   RUN THIS ONCE before submitting any fan-out. It is deliberately a separate
%   step from qs_run_sweep: if every array task tried to create the point list
%   on startup they would race each other and could disagree about which point
%   is which. With the list written first, each task only ever reads it.
%
%   force = true rewrites the point list. Existing partial results are NOT
%   deleted, but any result whose ID is no longer in the list will be dropped
%   at merge time -- so only force a rebuild when you meant to change the grid.

    if nargin < 3, force = false; end
    C = qs_config(caseId, studyId);

    for d = {C.paths.results_dir, C.paths.points_dir, C.paths.merged_dir, ...
             C.paths.hist_dir, C.paths.fig_dir, C.paths.log_dir}
        if exist(d{1},'dir')~=7, mkdir(d{1}); end
    end

    if exist(C.paths.pointlist,'file')==2 && ~force
        P = qs_pointlist('load', C.paths.pointlist);
        fprintf('qs_init_study: study already initialised -- %d points at %s\n', ...
            numel(P.id), C.paths.results_dir);
        fprintf('               (pass force=true to rebuild the point list)\n');
        return
    end

    mode = 'sweep';
    if isfield(C.study,'mode'), mode = C.study.mode; end

    switch lower(mode)
        case 'sweep'
            P = qs_make_points(C);
        case 'long'
            P = qs_make_long_points(C);
        otherwise
            error('qs_init_study: unknown study mode "%s" (expected "sweep" or "long")', mode);
    end

    qs_pointlist('save',  P, C.paths.pointlist);
    qs_pointlist('tocsv', P, fullfile(C.paths.results_dir,'pointlist.csv'));

    man = struct();
    man.case_id   = caseId;
    man.study_id  = studyId;
    man.mode      = mode;
    man.n_points  = numel(P.id);
    man.created   = C.meta.resolved_at;
    man.config    = strip_paths(C);
    qs_jsonwrite(C.paths.manifest, man);

    fprintf(['qs_init_study: %s / %s initialised\n' ...
             '  %d points, mode=%s\n' ...
             '  results  : %s\n' ...
             '  manifest : %s\n'], ...
        caseId, studyId, numel(P.id), mode, C.paths.results_dir, C.paths.manifest);
end

% ======================================================================
function P = qs_make_long_points(C)
%QS_MAKE_LONG_POINTS  Point list for a "long" study: named operating points
%   run for a long time. Reads C.study.points, a list of objects:
%       {"pc_kPa": 50.6, "dT_K": 9.5, "label": "case1_calibrated"}
%   If absent, the case nominal (pc_nom_Pa, dT_nom_K) is used.
    if isfield(C.study,'points') && ~isempty(C.study.points)
        Q = C.study.points;
        n = numel(Q);
        pc = zeros(n,1); dT = zeros(n,1);
        for k = 1:n
            if iscell(Q), q = Q{k}; else, q = Q(k); end
            pc(k) = q.pc_kPa*1e3;
            dT(k) = q.dT_K;
        end
    else
        pc = C.case.pc_nom_Pa;
        dT = C.case.dT_nom_K;
        fprintf('qs_init_study: no study.points given -- using the case nominal.\n');
    end

    tend = 20.0;  ttr = 2.0;
    if isfield(C.study,'duration')
        if isfield(C.study.duration,'base_t_end'),       tend = C.study.duration.base_t_end; end
        if isfield(C.study.duration,'base_t_transient'), ttr  = C.study.duration.base_t_transient; end
    end

    P = qs_pointlist('make', pc, dT, struct( ...
            't_end', tend, 't_trans', ttr, 'level', 0, ...
            'ic_tag', 'flat', 'origin', 'long', 'save_hist', true));
end

function C = strip_paths(C)
%STRIP_PATHS  Keep the manifest readable: absolute paths are environment,
%   not configuration, and they differ between the laptop and the cluster.
    if isfield(C,'paths')
        keep = struct('results_dir', C.paths.results_dir);
        C.paths = keep;
    end
end
