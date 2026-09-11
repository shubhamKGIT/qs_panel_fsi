function D = qs_merge(caseId, studyId, quiet)
%QS_MERGE  Collect every slice's partial results into one table.
%
%   D = QS_MERGE(caseId, studyId)
%
%   Reads results/<case>/<study>/points/partial_*.mat, concatenates the
%   records, de-duplicates by point ID (last write wins) and writes
%     merged/points.mat    struct of columns, one row per point
%     merged/points.csv    the same thing readable without MATLAB
%
%   Safe to run while slices are still going: whatever exists is merged and
%   the report says how many of the study's points are still missing.
%
%   The merge is a JOIN ON ID, not an assignment into a grid. That is what
%   lets refinement points be added later without renumbering anything.

    if nargin < 3, quiet = false; end
    C = qs_config(caseId, studyId);

    parts = dir(fullfile(C.paths.points_dir,'partial_*.mat'));
    assert(~isempty(parts), ...
        'qs_merge: no partial_*.mat in\n  %s\nRun the sweep first.', C.paths.points_dir);

    all_R = {};
    for k = 1:numel(parts)
        S = load(fullfile(parts(k).folder, parts(k).name));
        if ~isfield(S,'Rs'), continue; end
        all_R = [all_R, S.Rs(:).']; %#ok<AGROW>
        if ~quiet
            fprintf('  merged %-28s %4d points\n', parts(k).name, numel(S.Rs));
        end
    end
    assert(~isempty(all_R), 'qs_merge: partial files contained no records');

    % ---- de-duplicate by ID, keeping the last occurrence -------------------
    ids = cellfun(@(r) r.id, all_R, 'UniformOutput', false);
    [uid, ia] = unique(ids(:), 'last');
    R = all_R(ia);
    ndup = numel(ids) - numel(uid);

    % ---- flatten to columns ------------------------------------------------
    n = numel(R);
    D = struct();
    D.id = uid(:);
    num_fields = {'pc_Pa','dT_K','t_end','t_trans','level','Amp_wh','MeanPeak', ...
                  'StdPeak','MeanRMSE','StdRMSE','Freq','pfrac','signflip', ...
                  'label','label_start','trans_time','trans_t_lo','trans_t_hi','trans_from','trans_to', ...
                  'n_changes','amp_first','amp_last','wall_s'};
    log_fields = {'nonstationary','divergent','near_threshold','short_record', ...
                  'trans_flag','ok'};
    str_fields = {'ic_tag','label_name','errmsg'};

    for f = num_fields, D.(f{1}) = nan(n,1); end
    for f = log_fields, D.(f{1}) = false(n,1); end
    for f = str_fields, D.(f{1}) = repmat({''},n,1); end
    D.W = cell(n,1);

    for k = 1:n
        r = R{k};
        for f = num_fields, if isfield(r,f{1}) && ~isempty(r.(f{1})), D.(f{1})(k) = r.(f{1}); end, end
        for f = log_fields, if isfield(r,f{1}) && ~isempty(r.(f{1})), D.(f{1})(k) = logical(r.(f{1})); end, end
        for f = str_fields, if isfield(r,f{1}) && ~isempty(r.(f{1})), D.(f{1}){k} = r.(f{1}); end, end
        if isfield(r,'W'), D.W{k} = r.W; end
    end

    % ---- sort by (dT, pc) so the table reads like the map -------------------
    [~, ord] = sortrows([D.dT_K, D.pc_Pa]);
    D = reorder(D, ord);

    % ---- coverage against the study's point list ----------------------------
    P = qs_pointlist('load', C.paths.pointlist);
    have = ismember(P.id, D.id);
    D.meta = struct('case_id',caseId, 'study_id',studyId, ...
                    'n_points_in_list',numel(P.id), 'n_merged',n, ...
                    'n_missing',sum(~have), 'n_duplicates_dropped',ndup, ...
                    'git',C.meta.git_commit, 'merged_at',C.meta.resolved_at);
    D.missing_ids = P.id(~have);

    if exist(C.paths.merged_dir,'dir')~=7, mkdir(C.paths.merged_dir); end
    save(fullfile(C.paths.merged_dir,'points.mat'), 'D', qs_matver());
    write_csv(fullfile(C.paths.merged_dir,'points.csv'), D);

    if ~quiet
        qs_report(D);
        fprintf('  wrote %s\n', fullfile(C.paths.merged_dir,'points.mat'));
    end
end

% ======================================================================
function D = reorder(D, ord)
    f = fieldnames(D);
    for k = 1:numel(f)
        v = D.(f{k});
        if (isnumeric(v)||islogical(v)||iscell(v)) && numel(v)==numel(ord)
            D.(f{k}) = v(ord);
        end
    end
end

function write_csv(fname, D)
    fid = fopen(fname,'w');
    cols = {'id','pc_Pa','dT_K','t_end','level','label','label_name','Amp_wh', ...
            'Freq','pfrac','MeanPeak','StdPeak','MeanRMSE','StdRMSE', ...
            'nonstationary','divergent','near_threshold','trans_flag','trans_time', ...
            'trans_from','trans_to','ok','wall_s'};
    fprintf(fid,'%s\n', strjoin(cols,','));
    for k = 1:numel(D.id)
        vals = cell(1,numel(cols));
        for c = 1:numel(cols)
            v = D.(cols{c});
            if iscell(v)
                vals{c} = v{k};
            elseif islogical(v)
                vals{c} = sprintf('%d', v(k));
            else
                vals{c} = sprintf('%.6g', v(k));
            end
        end
        fprintf(fid,'%s\n', strjoin(vals,','));
    end
    fclose(fid);
end
