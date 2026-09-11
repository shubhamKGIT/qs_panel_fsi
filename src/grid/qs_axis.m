function v = qs_axis(spec, name)
%QS_AXIS  Build one parameter axis from a JSON grid specification.
%
%   Two accepted forms (both may carry an "extra" list of explicit values):
%
%   1) single uniform span
%        {"lo": 48.5, "hi": 55.0, "d": 0.25}
%
%   2) piecewise segments  -- this is what reproduces the existing
%      coarse/fine/coarse grids exactly
%        {"segments": [ {"lo":48.5,"hi":49.0,"d":0.5},
%                       {"lo":49.5,"hi":53.5,"d":0.1},
%                       {"lo":54.0,"hi":55.0,"d":0.5} ],
%         "extra": [51.927] }
%
%   Values are rounded to 1e-3 of the axis unit and uniqued, matching the
%   round(...,3) the original drivers used, so a study file that names the
%   old segments reproduces the old grid value-for-value.
%
%   Returns a sorted column vector in the axis' own unit (kPa or K).

    if nargin < 2, name = 'axis'; end
    assert(isstruct(spec), 'qs_axis: %s spec must be a JSON object', name);

    v = [];
    if isfield(spec,'segments') && ~isempty(spec.segments)
        S = spec.segments;
        if isstruct(S)                     % jsondecode gives a struct array
            for k = 1:numel(S), v = [v, seg(S(k), name)]; end %#ok<AGROW>
        elseif iscell(S)                   % ... or a cell array if shapes differ
            for k = 1:numel(S), v = [v, seg(S{k}, name)]; end %#ok<AGROW>
        else
            error('qs_axis: %s.segments has an unexpected type', name);
        end
    elseif isfield(spec,'lo')
        v = seg(spec, name);
    elseif ~isfield(spec,'extra')
        error('qs_axis: %s needs either "lo/hi/d", "segments", or "extra"', name);
    end

    if isfield(spec,'extra') && ~isempty(spec.extra)
        v = [v, spec.extra(:)'];
    end

    assert(~isempty(v), 'qs_axis: %s produced no values', name);
    v = unique(round(v(:)*1e3)/1e3);
end

% ----------------------------------------------------------------------
function s = seg(S, name)
    for f = {'lo','hi','d'}
        assert(isfield(S,f{1}), 'qs_axis: %s segment is missing "%s"', name, f{1});
    end
    assert(S.d > 0,      'qs_axis: %s segment step d must be > 0', name);
    assert(S.hi >= S.lo, 'qs_axis: %s segment needs hi >= lo', name);

    % Plain colon semantics: when the step does not divide the span exactly,
    % the last value falls SHORT of hi and hi itself is not included. This
    % matches what the original drivers got from lo:d:hi -- e.g.
    % 11.5:0.15:15.5 ends at 15.40, giving 27 values, not 28.
    % Set "include_hi": true on the segment, or list it under "extra", to add
    % the endpoint deliberately (which is what the case-3 driver did for
    % 74.5 kPa and 16.1 K).
    s = S.lo : S.d : S.hi;
    if isfield(S,'include_hi') && S.include_hi && abs(s(end) - S.hi) > 1e-9
        s = [s, S.hi];
    end
end
