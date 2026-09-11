function W = qs_windows(tvec, wc, h, C)
%QS_WINDOWS  Sliding-window signature of a centre-deflection history.
%
%   W = QS_WINDOWS(tvec, wc, h, C)
%
%   This is the machinery that makes a mid-run regime change visible. Instead
%   of one (transient, statistics) split, the whole record -- transient
%   included -- is cut into overlapping windows and each window is reduced to
%   a few numbers and a label. A point whose label changes part way through
%   and stays changed is a TRANSITION (see qs_detect_transition).
%
%   Why it matters here specifically: a calibrated case-1 point was observed
%   to snap from oscillation to static at ~7 s. Any regime label taken from a
%   3 s run near a boundary is therefore provisional until the windows say
%   the label had settled.
%
%   Config used: C.windows.T_win_s, C.windows.overlap, C.classify.*
%
%   OUTPUT (struct of equal-length columns, one row per window)
%     .t0, .t1, .tmid   window start / end / centre time [s]
%     .wmean            mean centre deflection           [m]   (which well)
%     .amp_wh           peak-to-peak / (2h)              [-]
%     .fdom             dominant frequency               [Hz]
%     .pfrac            peak power fraction              [-]
%     .label            window label code (see qs_label_name)
%     .n                number of windows
%   Storage: a 20 s record at T_win = 0.5 s and 50% overlap gives 79 windows,
%   i.e. a few hundred bytes -- cheap enough to keep for EVERY swept point,
%   which is what lets the sweep report transitions without keeping histories.

    Fs   = 1/C.duration.dt_out;
    Twin = C.windows.T_win_s;
    ov   = C.windows.overlap;
    assert(ov >= 0 && ov < 1, 'qs_windows: overlap must be in [0,1)');

    tvec = tvec(:);  wc = wc(:);
    T    = tvec(end) - tvec(1);

    nw   = max(1, round(Twin*Fs));            % samples per window
    step = max(1, round(nw*(1-ov)));          % samples between window starts
    if nw > numel(wc)                         % record shorter than one window
        nw = numel(wc);  step = nw;
    end
    starts = 1:step:(numel(wc)-nw+1);
    if isempty(starts), starts = 1; end
    n = numel(starts);

    W.n      = n;
    W.t0     = zeros(n,1);  W.t1    = zeros(n,1);  W.tmid  = zeros(n,1);
    W.wmean  = zeros(n,1);  W.amp_wh= zeros(n,1);
    W.fdom   = nan(n,1);    W.pfrac = nan(n,1);
    W.label  = zeros(n,1);

    bf = 0.05;
    if isfield(C.classify,'peak_band_frac'), bf = C.classify.peak_band_frac; end

    for k = 1:n
        a = starts(k);  b = a + nw - 1;
        seg = wc(a:b);
        W.t0(k)   = tvec(a);
        W.t1(k)   = tvec(b);
        W.tmid(k) = 0.5*(tvec(a)+tvec(b));
        W.wmean(k)  = mean(seg);
        W.amp_wh(k) = (max(seg)-min(seg)) / (2*h);
        S = qs_spectrum(seg, Fs, bf);
        W.fdom(k)  = S.fpk;
        W.pfrac(k) = S.pfrac;
        W.label(k) = qs_label_window(W.amp_wh(k), W.pfrac(k), C);
    end

    W.T_win_s  = Twin;
    W.overlap  = ov;
    W.T_record = T;
end
