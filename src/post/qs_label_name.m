function s = qs_label_name(code)
%QS_LABEL_NAME  Human-readable name for a regime code (scalar or array).
%   0 unknown | 1 static | 2 LCO | 3 broadband | 4 divergent | 5 nonstationary
    names = {'unknown','static','LCO','broadband','divergent','nonstationary'};
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
