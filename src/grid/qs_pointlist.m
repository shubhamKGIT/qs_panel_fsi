function varargout = qs_pointlist(op, varargin)
%QS_POINTLIST  The point list: the framework's replacement for a fixed grid.
%
%   A point list is a struct of equal-length columns, one row per run:
%
%     P.id        {n x 1} cellstr   deterministic ID from (pc,dT), see qs_id
%     P.pc_Pa     [n x 1]           cavity / back pressure     [Pa]
%     P.dT_K      [n x 1]           temperature rise           [K]
%     P.t_end     [n x 1]           run length for THIS point  [s]
%     P.t_trans   [n x 1]           transient discarded        [s]
%     P.level     [n x 1]           0 = coarse, 1.. = refinement pass
%     P.ic_tag    {n x 1} cellstr   'flat' | 'flat_neg' | 'from:<parent id>'
%     P.save_hist [n x 1] logical   force saving the full centre history
%     P.origin    {n x 1} cellstr   'coarse' | 'refine' | 'manual' | 'long'
%
%   Why a list and not an (npc x ndT) array:
%     * refinement APPENDS rows; nothing is renumbered, so old results stay valid
%     * per-point t_end is just a column, so "edge cases run longer" needs no
%       special mechanism
%     * merge is a join on ID, and a rectangular grid is rebuilt only for plotting
%
%   OPERATIONS
%     P   = qs_pointlist('empty')
%     P   = qs_pointlist('make', pc_Pa, dT_K, opts)   opts: t_end,t_trans,level,ic_tag,origin,save_hist
%     P   = qs_pointlist('append', P, Q)              appends Q, dropping IDs already in P
%     tf  = qs_pointlist('has', P, id)                logical, per element of cellstr id
%     idx = qs_pointlist('find', P, id)               row index (0 if absent)
%     S   = qs_pointlist('slice', P, k, n)            rows k:n:end (round-robin, 1-based)
%     R   = qs_pointlist('rows', P, idx)              subset of rows
%     n   = qs_pointlist('count', P)
%     qs_pointlist('save', P, fname)
%     P   = qs_pointlist('load', fname)
%     qs_pointlist('tocsv', P, fname)

    switch lower(op)
        case 'empty',  varargout{1} = pl_empty();
        case 'make',   varargout{1} = pl_make(varargin{:});
        case 'append', varargout{1} = pl_append(varargin{:});
        case 'has',    varargout{1} = pl_has(varargin{:});
        case 'find',   varargout{1} = pl_find(varargin{:});
        case 'slice',  varargout{1} = pl_slice(varargin{:});
        case 'rows',   varargout{1} = pl_rows(varargin{:});
        case 'count',  varargout{1} = numel(varargin{1}.id);
        case 'save',   pl_save(varargin{:});
        case 'load',   varargout{1} = pl_load(varargin{:});
        case 'tocsv',  pl_tocsv(varargin{:});
        otherwise, error('qs_pointlist: unknown operation "%s"', op);
    end
end

% ======================================================================
function P = pl_empty()
    P = struct('id',{{}}, 'pc_Pa',[], 'dT_K',[], 't_end',[], 't_trans',[], ...
               'level',[], 'ic_tag',{{}}, 'save_hist',logical([]), 'origin',{{}});
end

function P = pl_make(pc_Pa, dT_K, opts)
    if nargin < 3, opts = struct(); end
    pc = pc_Pa(:);  dT = dT_K(:);
    n  = numel(pc);
    assert(numel(dT)==n, 'qs_pointlist:make  pc_Pa and dT_K must be the same length');

    g = @(f,d) getdef(opts,f,d);
    P.pc_Pa     = pc;
    P.dT_K      = dT;
    P.t_end     = expand(g('t_end',   3.0), n);
    P.t_trans   = expand(g('t_trans', 1.5), n);
    P.level     = expand(g('level',   0),   n);
    P.ic_tag    = expandc(g('ic_tag', 'flat'),   n);
    P.save_hist = logical(expand(g('save_hist', false), n));
    P.origin    = expandc(g('origin', 'coarse'), n);
    % ID depends on how the point is run, not only where it sits -- see qs_id
    P.id        = cellify(qs_id(pc, dT, P.t_end, P.ic_tag), n);

    % Guard the one mistake that silently wastes a whole sweep.
    bad = P.t_trans >= P.t_end;
    assert(~any(bad), ['qs_pointlist:make  t_trans >= t_end for %d point(s) ' ...
        '(first: t_end=%g, t_trans=%g). No samples would be left for statistics.'], ...
        sum(bad), P.t_end(find(bad,1)), P.t_trans(find(bad,1)));
end

function P = pl_append(P, Q)
    if isempty(Q.id), return; end
    keep = ~pl_has(P, Q.id);
    if ~any(keep), return; end
    Q = pl_rows(Q, find(keep));
    f = fieldnames(P);
    for k = 1:numel(f)
        P.(f{k}) = [P.(f{k}); Q.(f{k})];
    end
end

function tf = pl_has(P, id)
    if ischar(id), id = {id}; end
    if isempty(P.id), tf = false(numel(id),1); return; end
    tf = ismember(id(:), P.id(:));
end

function idx = pl_find(P, id)
    if ischar(id), id = {id}; end
    idx = zeros(numel(id),1);
    if isempty(P.id), return; end
    [tf, loc] = ismember(id(:), P.id(:));
    idx(tf) = loc(tf);
end

function S = pl_slice(P, k, n)
    assert(k>=1 && k<=n, 'qs_pointlist:slice  need 1 <= k <= n (got k=%d, n=%d)', k, n);
    S = pl_rows(P, (k:n:numel(P.id))');
end

function R = pl_rows(P, idx)
    idx = idx(:);
    f = fieldnames(P);
    R = struct();
    for k = 1:numel(f)
        v = P.(f{k});
        if isempty(v), R.(f{k}) = v; else, R.(f{k}) = v(idx); end
    end
end

function pl_save(P, fname)
    d = fileparts(fname);
    if ~isempty(d) && exist(d,'dir')~=7, mkdir(d); end
    save(fname, 'P', '-v7');
end

function P = pl_load(fname)
    assert(exist(fname,'file')==2, 'qs_pointlist:load  not found: %s', fname);
    S = load(fname);
    P = S.P;
end

function pl_tocsv(P, fname)
%PL_TOCSV  Plain-text view of the point list, readable without MATLAB.
    d = fileparts(fname);
    if ~isempty(d) && exist(d,'dir')~=7, mkdir(d); end
    fid = fopen(fname,'w');
    fprintf(fid,'id,pc_Pa,dT_K,t_end,t_trans,level,ic_tag,save_hist,origin\n');
    for k = 1:numel(P.id)
        fprintf(fid,'%s,%.1f,%.4f,%.4f,%.4f,%d,%s,%d,%s\n', ...
            P.id{k}, P.pc_Pa(k), P.dT_K(k), P.t_end(k), P.t_trans(k), ...
            P.level(k), P.ic_tag{k}, P.save_hist(k), P.origin{k});
    end
    fclose(fid);
end

% ---------------------------------------------------------------- helpers
function v = getdef(S, f, d)
%GETDEF  Option lookup with a default.
%   The isscalar guard is load-bearing: struct('ic_tag', someCellArray)
%   silently builds a STRUCT ARRAY, and S.(f) on one of those returns a
%   comma-separated list that blows up in confusing places. Catch it here
%   with a message that says what actually went wrong.
    v = d;
    if ~isstruct(S), return; end
    assert(isscalar(S), ...
        ['qs_pointlist: the options argument is a %dx%d struct array, not a ' ...
         'single struct. A cell-array value passed to struct() does that -- ' ...
         'wrap it in braces, e.g. struct(''ic_tag'', {myCellstr}).'], size(S,1), size(S,2));
    if isfield(S,f) && ~isempty(S.(f)), v = S.(f); end
end

function v = expand(v, n)
    if isscalar(v), v = repmat(v, n, 1); else, v = v(:); end
    assert(numel(v)==n, 'qs_pointlist: column length mismatch (%d vs %d)', numel(v), n);
end

function c = expandc(c, n)
    if ischar(c), c = repmat({c}, n, 1); else, c = c(:); end
    assert(numel(c)==n, 'qs_pointlist: column length mismatch (%d vs %d)', numel(c), n);
end

function c = cellify(c, n)
    if ischar(c), c = {c}; end
    c = c(:);
    assert(numel(c)==n, 'qs_pointlist: id length mismatch');
end
