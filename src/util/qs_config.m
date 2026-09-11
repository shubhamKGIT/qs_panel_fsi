function C = qs_config(caseId, studyId, overrides)
%QS_CONFIG  Resolve the full configuration for one (case, study) pair.
%
%   C = QS_CONFIG(caseId, studyId)
%   C = QS_CONFIG(caseId, studyId, overrides)
%
%   caseId    name of a folder under  <root>/cases/     e.g. 'unheated_c1_NoSBLI_Periodic'
%   studyId   name of a file   under  <root>/studies/   e.g. 'sweep_coarse_v1'
%             (with or without the .json extension)
%   overrides optional struct merged in LAST, for one-off changes without
%             editing a file, e.g. struct('duration',struct('base_t_end',5))
%
%   Resolution order (later wins):
%       shared/defaults.json  <-  cases/<caseId>/case.json  <-  studies/<studyId>.json  <-  overrides
%
%   The resolved struct is what every other function reads, and a copy of it
%   is written into the results manifest. That is the whole point: a result
%   folder carries the exact settings that produced it, so nothing has to be
%   reconstructed from memory or from which folder the code was copied into.
%
%   OUTPUT (struct C)
%     C.case        case block (physics of the CFD case)
%     C.study       study block (what to do with it)
%     C.solver, C.aero, C.duration, C.ic, C.windows, C.classify,
%     C.store, C.duration_policy, C.refine, C.preflight   (resolved defaults)
%     C.paths       every path the run needs, absolute and already resolved
%     C.meta        framework version, git commit, timestamp, host

    if nargin < 3, overrides = struct(); end
    root = qs_root();

    % ---------- locate the three inputs --------------------------------
    defFile = fullfile(root,'shared','defaults.json');

    caseDir = fullfile(root,'cases',caseId);
    assert(exist(caseDir,'dir')==7, ...
        'qs_config: no case folder\n  %s\nAvailable: %s', caseDir, qs_listdir(fullfile(root,'cases')));
    caseFile = fullfile(caseDir,'case.json');

    if numel(studyId)>5 && strcmp(studyId(end-4:end),'.json')
        studyId = studyId(1:end-5);
    end
    studyFile = fullfile(root,'studies',[studyId '.json']);
    assert(exist(studyFile,'file')==2, ...
        'qs_config: no study file\n  %s\nAvailable: %s', studyFile, qs_listdir(fullfile(root,'studies')));

    % ---------- merge ---------------------------------------------------
    C = qs_jsonread(defFile);
    Ccase  = qs_jsonread(caseFile);
    Cstudy = qs_jsonread(studyFile);

    C = qs_mergestruct(C, rmfield_if(Ccase,  {'case_id'}));
    C = qs_mergestruct(C, rmfield_if(Cstudy, {'study_id'}));
    C = qs_mergestruct(C, overrides);

    % case / study blocks are kept whole as well, so nothing is lost
    C.case  = Ccase;
    C.study = Cstudy;
    C.case.case_id   = caseId;
    C.study.study_id = studyId;

    if isfield(Ccase,'case_id') && ~strcmp(Ccase.case_id, caseId)
        warning('qs_config:caseIdMismatch', ...
            'case.json says case_id="%s" but the folder is "%s". Folder name wins.', ...
            Ccase.case_id, caseId);
    end
    if isfield(Cstudy,'study_id') && ~strcmp(Cstudy.study_id, studyId)
        warning('qs_config:studyIdMismatch', ...
            'study json says study_id="%s" but the file is "%s.json". File name wins.', ...
            Cstudy.study_id, studyId);
    end

    % ---------- paths ----------------------------------------------------
    C.paths.root       = root;
    C.paths.case_dir   = caseDir;
    C.paths.study_file = studyFile;
    C.paths.rom_file   = resolve(caseDir, req(C.case,'rom_file',  caseFile));
    C.paths.aero_file  = resolve(caseDir, req(C.case,'aero_file', caseFile));
    C.paths.ref_file   = resolve(caseDir, req(C.case,'ref_file',  caseFile));

    C.paths.results_dir = fullfile(root,'results',caseId,studyId);
    C.paths.points_dir  = fullfile(C.paths.results_dir,'points');
    C.paths.merged_dir  = fullfile(C.paths.results_dir,'merged');
    C.paths.hist_dir    = fullfile(C.paths.results_dir,'histories');
    C.paths.fig_dir     = fullfile(C.paths.results_dir,'figures');
    C.paths.log_dir     = fullfile(C.paths.results_dir,'logs');
    C.paths.manifest    = fullfile(C.paths.results_dir,'manifest.json');
    C.paths.pointlist   = fullfile(C.paths.results_dir,'pointlist.mat');

    % ---------- validate the parts that silently ruin a run --------------
    must(C,'case.pinf_Pa');     must(C,'case.Minf');
    must(C,'case.pc_nom_Pa');   must(C,'case.dT_nom_K');
    must(C,'case.panel');       must(C,'duration.dt_out');
    assert(C.duration.dt_out > 0, 'qs_config: duration.dt_out must be positive');
    assert(isfield(C.case.panel,'h_m') && C.case.panel.h_m > 0, ...
        'qs_config: case.panel.h_m (panel thickness) is required -- it sets the w/h amplitude scale');

    for k = {'rom_file','aero_file','ref_file'}
        p = C.paths.(k{1});
        if exist(p,'file')~=2
            warning('qs_config:missingData','%s does not exist yet:\n  %s', k{1}, p);
        end
    end

    % ---------- provenance ------------------------------------------------
    C.meta.framework_version = '1.0.0';
    C.meta.resolved_at       = datestr(now,'yyyy-mm-ddTHH:MM:SS');  %#ok<TNOW1,DATST>
    C.meta.git_commit        = qs_gitcommit(root);
    try, C.meta.host = char(java.net.InetAddress.getLocalHost.getHostName); catch, C.meta.host = getenv('HOSTNAME'); end
    C.meta.matlab_version    = version();
end

% ======================================================================
function s = qs_listdir(d)
    if exist(d,'dir')~=7, s = '(none)'; return; end
    L = dir(d);  names = {};
    for k = 1:numel(L)
        if L(k).name(1)=='.', continue; end
        names{end+1} = L(k).name; %#ok<AGROW>
    end
    if isempty(names), s = '(none)'; else, s = strjoin(names, ', '); end
end

function S = rmfield_if(S, flds)
    for k = 1:numel(flds)
        if isfield(S, flds{k}), S = rmfield(S, flds{k}); end
    end
end

function v = req(S, key, srcfile)
    assert(isfield(S,key), 'qs_config: "%s" is missing from %s', key, srcfile);
    v = S.(key);
end

function p = resolve(baseDir, p)
%RESOLVE  Absolute paths pass through; relative ones hang off the case folder.
    if isempty(p), return; end
    if p(1)=='/' || p(1)=='~' || (numel(p)>1 && p(2)==':')
        return
    end
    p = fullfile(baseDir, p);
end

function must(C, dotted)
    parts = strsplit(dotted,'.');
    S = C;
    for k = 1:numel(parts)
        assert(isstruct(S) && isfield(S,parts{k}), ...
            'qs_config: required setting "%s" is missing', dotted);
        S = S.(parts{k});
    end
end

function h = qs_gitcommit(root)
%QS_GITCOMMIT  Short commit hash of the framework, or 'nogit'.
    h = 'nogit';
    try
        [st, out] = system(sprintf('git -C "%s" rev-parse --short HEAD 2>/dev/null', root));
        if st==0 && ~isempty(strtrim(out))
            h = strtrim(out);
            [~, dirtyOut] = system(sprintf('git -C "%s" status --porcelain 2>/dev/null', root));
            if ~isempty(strtrim(dirtyOut)), h = [h '-dirty']; end
        end
    catch
    end
end
