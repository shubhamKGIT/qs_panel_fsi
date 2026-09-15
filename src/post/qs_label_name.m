function s = qs_label_name(code)
%QS_LABEL_NAME  Human-readable name for a regime code (scalar or array).
%
%   0 unknown | 1 static | 2 LCO | 3 broadband | 4 divergent
%
%   THESE FIVE ARE THE COMPLETE SET. Nonstationarity is a FLAG on a run, not a
%   regime: a run can be static and nonstationary at once, which is the whole
%   point of keeping them separate. (An earlier version of this file listed a
%   code 5 "nonstationary" that nothing ever assigned; it is gone, so the label
%   legend and the code now agree.)
%
%   qs_codes() prints this legend together with every other coded value in the
%   framework, and derives the label part from this function so the two cannot
%   drift apart.
    names = {'unknown','static','LCO','broadband','divergent'};
    if isscalar(code)
        s = pick(names, code);
    else
        s = cell(size(code));
        for k = 1:numel(code), s{k} = pick(names, code(k)); end
    end
end

function s = pick(names, c)
    if ~isfinite(c) || c < 0 || c > numel(names)-1
        s = 'unknown';
    else
        s = names{c+1};
    end
end
