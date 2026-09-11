function code = qs_label_window(amp_wh, pfrac, C)
%QS_LABEL_WINDOW  Regime label for a single analysis window.
%
%   Codes (kept numeric so a whole label sequence costs one int8 per window):
%       0  unknown        spectrum could not be formed
%       1  static         amplitude below tolerance -- settled, no oscillation
%       2  LCO            oscillating, one dominant tone
%       3  broadband      oscillating, no dominant tone (chaotic-looking)
%   Codes 4 (divergent) and 5 (nonstationary) are decided over the WHOLE
%   record, not per window, and are assigned in qs_classify.
%
%   Thresholds come from C.classify:
%       amp_tol_wh        below this in w/h the window is "static"
%       peak_power_frac   above this peak power fraction the window is "LCO"

    ampTol = C.classify.amp_tol_wh;
    pTol   = C.classify.peak_power_frac;

    if ~isfinite(amp_wh)
        code = 0;
    elseif amp_wh < ampTol
        code = 1;                       % static
    elseif isfinite(pfrac) && pfrac >= pTol
        code = 2;                       % limit cycle
    else
        code = 3;                       % broadband / chaotic
    end
end
