function id = qs_id(pc_Pa, dT_K, t_end, ic_tag)
%QS_ID  Deterministic run identifier.
%
%   id = QS_ID(pc_Pa, dT_K, t_end, ic_tag)
%       e.g.  'pc051927_dT0012p760_t003p00_flat'
%
%   Scalar inputs -> char row vector. Vector inputs -> cellstr (n x 1).
%   t_end and ic_tag default to 3.0 s and 'flat' if omitted.
%
%   WHY ALL FOUR FIELDS ARE IN THE KEY
%   A "run" is defined by its operating point AND how it was run. Two runs at
%   the same (p_c, dT) but different lengths are different experiments, and
%   near a regime boundary they can legitimately disagree -- that disagreement
%   is the finding, not a collision to be hidden. Keying on all four means:
%     * re-running a boundary point for 10 s ADDS a row rather than silently
%       overwriting the 3 s answer, so both are on the record
%     * a refinement point that lands on an existing value with the same
%       settings is skipped instead of recomputed
%     * results merge by ID, so appending points never renumbers anything
%     * a partial file is resumed by checking which IDs it already holds
%   qs_rasterize then picks, per (p_c, dT), the longest run available -- the
%   most trustworthy answer -- when it builds a map.
%
%   Resolution: 1 Pa in p_c, 1 mK in dT, 10 ms in t_end. All far finer than
%   the corresponding measurement uncertainties, so two points that differ by
%   less than one ID step are not physically distinguishable.

    if nargin < 3 || isempty(t_end),  t_end  = 3.0;    end
    if nargin < 4 || isempty(ic_tag), ic_tag = 'flat'; end

    pc = round(pc_Pa(:));
    dT = dT_K(:);
    n  = numel(pc);
    assert(numel(dT)==n, 'qs_id: pc_Pa and dT_K must have the same length');

    te = t_end(:);
    if isscalar(te), te = repmat(te, n, 1); end
    assert(numel(te)==n, 'qs_id: t_end must be scalar or the same length as pc_Pa');

    if ischar(ic_tag), ic = repmat({ic_tag}, n, 1); else, ic = ic_tag(:); end
    assert(numel(ic)==n, 'qs_id: ic_tag must be a single string or the same length as pc_Pa');

    out = cell(n,1);
    for k = 1:n
        s = sprintf('pc%06d_dT%08.3f_t%06.2f_%s', pc(k), dT(k), te(k), tagcode(ic{k}));
        s = strrep(s, '.', 'p');
        s = strrep(s, '-', 'm');
        out{k} = s;
    end

    if n == 1, id = out{1}; else, id = out; end
end

% ----------------------------------------------------------------------
function c = tagcode(tag)
%TAGCODE  Short, filename-safe form of an initial-condition tag.
    c = regexprep(tag, '[^A-Za-z0-9_]', '');
    if numel(c) > 16, c = c(1:16); end
    if isempty(c), c = 'ic'; end
end
