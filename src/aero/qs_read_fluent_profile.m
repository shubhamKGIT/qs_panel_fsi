function S = qs_read_fluent_profile(fname)
%QS_READ_FLUENT_PROFILE  Parse a Fluent ASCII profile (.prof) export.
%
%   S = QS_READ_FLUENT_PROFILE(fname) reads
%
%       ((name point N)
%       (field1
%       val1
%       ...
%       valN
%       )
%       (field2
%       ...
%       )
%       )
%
%   and returns a struct with one N x 1 column per field. Field names are
%   sanitised to valid identifiers ("mach-number" -> "mach_number").
%
%   A .csv with a header row is also accepted, so a profile can be converted
%   outside MATLAB if that is more convenient -- the column names become the
%   struct fields and must include at least x, y and the quantities the
%   aero build needs.
%
%   Reading is vectorised (one str2double on the whole block) rather than
%   value-by-value, which matters for the 18,743-point edge export.

    assert(exist(fname,'file')==2, 'qs_read_fluent_profile: not found:\n  %s', fname);

    [~,~,ext] = fileparts(fname);
    if strcmpi(ext,'.csv')
        S = read_csv(fname);
        return
    end

    txt   = fileread(fname);
    lines = regexp(txt, '\r?\n', 'split');

    S = struct();
    curName = '';
    buf = {};
    for k = 1:numel(lines)
        ln = strtrim(lines{k});
        if isempty(ln)
            continue
        elseif numel(ln)>=2 && strcmp(ln(1:2),'((')
            continue                                    % header
        elseif ln(1)=='('
            S = flush(S, curName, buf);
            curName = sanitize(strtrim(ln(2:end)));
            buf = {};
        elseif strcmp(ln,')')
            continue
        else
            buf{end+1} = ln; %#ok<AGROW>
        end
    end
    S = flush(S, curName, buf);

    f = fieldnames(S);
    assert(~isempty(f), 'qs_read_fluent_profile: no fields parsed from %s', fname);
    n = numel(S.(f{1}));
    for k = 2:numel(f)
        assert(numel(S.(f{k}))==n, ...
            'qs_read_fluent_profile: field "%s" has %d values but "%s" has %d', ...
            f{k}, numel(S.(f{k})), f{1}, n);
    end
end

% ======================================================================
function S = flush(S, name, buf)
    if isempty(name) || isempty(buf), return; end
    S.(name) = str2double(buf(:));
end

function s = sanitize(s)
    s = regexprep(s, '[^A-Za-z0-9_]', '_');
    if isempty(s) || ~isletter(s(1)), s = ['f_' s]; end
end

function S = read_csv(fname)
    fid = fopen(fname,'r');
    hdr = strtrim(fgetl(fid));
    fclose(fid);
    names = strtrim(strsplit(hdr, ','));
    for k = 1:numel(names), names{k} = sanitize(names{k}); end
    Mx = dlmread(fname, ',', 1, 0);
    assert(size(Mx,2)==numel(names), ...
        'qs_read_fluent_profile: %s has %d header names but %d columns', ...
        fname, numel(names), size(Mx,2));
    S = struct();
    for k = 1:numel(names)
        S.(names{k}) = Mx(:,k);
    end
end
