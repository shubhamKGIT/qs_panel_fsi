function qs_run_sweep(caseId, studyId, slice, nslice)
%QS_RUN_SWEEP  Run one slice of a study's point list.
%
%   qs_run_sweep(caseId, studyId)                 whole list, this process
%   qs_run_sweep(caseId, studyId, slice, nslice)  rows slice:nslice:end
%
%   Slicing is round-robin over the point list, the same scheme the original
%   sweep used over linear grid indices, so neighbouring points land on
%   different workers and every task finishes in comparable time.
%
%   RESUMABLE. Each slice writes partial_<slice>_of_<nslice>.mat under
%   results/<case>/<study>/points/, checkpointed every C.checkpoint_every
%   points. Re-running the same slice skips IDs that file already holds, so a
%   killed job is restarted by resubmitting the identical command.
%
%   The study must have been created first with qs_init_study.

    if nargin < 3 || isempty(slice),  slice  = 1; end
    if nargin < 4 || isempty(nslice), nslice = 1; end

    C = qs_config(caseId, studyId);
    assert(exist(C.paths.pointlist,'file')==2, ...
        ['qs_run_sweep: study not initialised.\n' ...
         '  Run  qs_init_study(''%s'',''%s'')  once before submitting.'], caseId, studyId);

    P  = qs_pointlist('load', C.paths.pointlist);
    S  = qs_pointlist('slice', P, slice, nslice);
    nP = numel(S.id);

    fprintf('=== qs_run_sweep  %s / %s  slice %d of %d ===\n', caseId, studyId, slice, nslice);
    fprintf('    %d of %d points on this task\n', nP, numel(P.id));
    if nP == 0
        fprintf('    nothing to do.\n');  return
    end

    % ---- model (built once, reused for every point) -----------------------
    M = qs_build_model(C);
    if isfield(C.preflight,'run_on_sweep') && C.preflight.run_on_sweep
        qs_preflight(M, C);
    end

    % ---- resume from an existing partial ----------------------------------
    partfile = fullfile(C.paths.points_dir, sprintf('partial_%d_of_%d.mat', slice, nslice));
    Rs = {};  done = {};
    if exist(partfile,'file')==2
        old = load(partfile);
        if isfield(old,'Rs')
            Rs   = old.Rs(:).';
            done = cellfun(@(r) r.id, Rs, 'UniformOutput', false);
            fprintf('    resuming: %d point(s) already in %s\n', numel(Rs), partfile);
        end
    end

    every = 5;
    if isfield(C,'checkpoint_every'), every = C.checkpoint_every; end
    t0 = tic;  nrun = 0;  nfail = 0;

    for k = 1:nP
        pt = qs_pointlist('rows', S, k);
        if any(strcmp(pt.id{1}, done))
            continue
        end

        R = qs_run_point(M, C, pt);
        Rs{end+1} = R; %#ok<AGROW>
        nrun = nrun + 1;
        if ~R.ok, nfail = nfail + 1; end

        % full centre history, when the storage policy kept one
        if ~isempty(R.wc)
            save_history(C, R);
            Rs{end}.wc = [];  Rs{end}.tvec = [];   % not duplicated in the partial
        end
        if ~isempty(R.yend)
            save_state(C, R);
            Rs{end}.yend = [];
        end

        if mod(nrun, every)==0 || k==nP
            saveRs(partfile, Rs, C, slice, nslice);
            fprintf('    %4d/%4d  p_c=%7.2f kPa  dT=%6.3f K  %-13s amp=%6.3f  %s[%.1f min]\n', ...
                k, nP, R.pc_Pa/1e3, R.dT_K, R.label_name, R.Amp_wh, ...
                trans_note(R), toc(t0)/60);
        end
    end

    saveRs(partfile, Rs, C, slice, nslice);
    fprintf('=== slice %d/%d done: %d run, %d failed, %.1f min -> %s\n', ...
        slice, nslice, nrun, nfail, toc(t0)/60, partfile);
end

% ======================================================================
function saveRs(partfile, Rs, C, slice, nslice) %#ok<INUSD>
    d = fileparts(partfile);
    if exist(d,'dir')~=7, mkdir(d); end
    meta = struct('slice',slice,'nslice',nslice,'case_id',C.case.case_id, ...
                  'study_id',C.study.study_id,'git',C.meta.git_commit);
    save(partfile, 'Rs', 'meta', qs_matver());
end

function save_history(C, R)
    if exist(C.paths.hist_dir,'dir')~=7, mkdir(C.paths.hist_dir); end
    f = fullfile(C.paths.hist_dir, [R.id '_hist.mat']);
    tvec = R.tvec;  wc = R.wc;  W = R.W; %#ok<NASGU>
    id = R.id;  pc_Pa = R.pc_Pa;  dT_K = R.dT_K; %#ok<NASGU>
    label = R.label;  trans_flag = R.trans_flag;  trans_time = R.trans_time; %#ok<NASGU>
    save(f, 'id','pc_Pa','dT_K','tvec','wc','W','label','trans_flag','trans_time',qs_matver());
end

function save_state(C, R)
    if exist(C.paths.hist_dir,'dir')~=7, mkdir(C.paths.hist_dir); end
    f = fullfile(C.paths.hist_dir, [R.id '_state.mat']);
    yend = R.yend;  id = R.id; %#ok<NASGU>
    save(f, 'yend', 'id', '-v7');
end

function s = trans_note(R)
    if R.trans_flag
        s = sprintf('TRANSITION %s->%s @ %.2f s  ', ...
            qs_label_name(R.trans_from), qs_label_name(R.trans_to), R.trans_time);
    elseif R.nonstationary
        s = 'nonstationary  ';
    else
        s = '';
    end
end
