function R = t05_windows_classify()
%T05  Windows, classification and transition detection, on synthetic signals
%     whose answers are known by construction.
%
%   This is the test for the genuinely NEW machinery -- the part with no
%   equivalent in the old sweep, so nothing to regress against. Each case
%   below is a signal where the correct answer is obvious, and the test is
%   that the code returns it.

    R = struct('pass',true,'skip',false,'msg',{{}});
    function chk(cond, fmt, varargin)
        if ~cond
            R.pass = false;
            R.msg{end+1} = sprintf(fmt, varargin{:});
            fprintf('    FAIL  %s\n', R.msg{end});
        end
    end

    C = qs_config('unheated_c1_NoSBLI_Periodic','regression_unheated_c1');
    C.duration.dt_out = 1e-4;
    C.windows = struct('T_win_s',0.5,'overlap',0.5,'settle_windows',1,'min_persist_windows',2);
    C.classify.amp_tol_wh = 0.10;
    C.classify.peak_power_frac = 0.35;
    C.classify.stationarity_tol = 0.15;
    C.classify.growth_tol = 3.0;
    C.classify.near_threshold_band = 0.5;

    h  = 6.35e-4;
    dt = C.duration.dt_out;
    t  = (0:dt:8).';
    f0 = 263;

    % ---------- 1. a settled flat panel is static -------------------------
    wc = 0.01*h*sin(2*pi*f0*t);                     % amplitude 0.005 w/h << 0.1
    [K,W] = run_case(t, wc, h, C);
    chk(K.label==1, 'a tiny oscillation was labelled "%s", expected static', K.label_name);
    chk(~K.flags.nonstationary, 'a steady signal was flagged nonstationary');
    fprintf('    static           -> %s\n', K.label_name);

    % ---------- 1b. a dead panel is NOT "nonstationary" -------------------
    % A still panel's residual decays and its first/last amplitude ratio is
    % roundoff. Without an amplitude floor on the drift test this fires on
    % every static cell in a map, and each one then gets promoted to a 10 s
    % and then a 20 s re-run. This is the guard against that.
    wc = 1e-8*h*exp(-t/1.3).*sin(2*pi*f0*t) + 1e-9*h*sin(2*pi*37*t);
    [K,~] = run_case(t, wc, h, C);
    chk(K.label==1, 'a dead panel was labelled "%s", expected static', K.label_name);
    chk(~K.flags.nonstationary, ...
        ['a panel at rest (%.1e w/h) was flagged NONSTATIONARY. The drift test is ' ...
         'being applied to numerical residue; check classify.stationarity_amp_floor.'], K.amp_last);
    chk(~K.flags.drift_tested, 'the drift test should not even have been applied at %.1e w/h', K.amp_ref);
    fprintf('    dead panel       -> %s, nonstationary %d (drift test applied: %d)\n', ...
        K.label_name, K.flags.nonstationary, K.flags.drift_tested);

    % ---------- 1c. a borderline amplitude IS still judged on drift -------
    % 0.07 w/h is below the static tolerance but far above residue, so a run
    % that is genuinely still growing there must be caught.
    wc = h*(0.02 + 0.05*t/max(t)).*sin(2*pi*f0*t);
    [K,~] = run_case(t, wc, h, C);
    chk(K.flags.drift_tested, 'a %.3f w/h record was not drift-tested; the floor is too high', K.amp_ref);
    chk(K.flags.nonstationary, 'a steadily growing borderline record was not flagged nonstationary');
    fprintf('    borderline drift -> %s, amp %.4f -> %.4f w/h, nonstationary %d\n', ...
        K.label_name, K.amp_first, K.amp_last, K.flags.nonstationary);

    % ---------- 2. a clean tone is an LCO, at the right frequency ---------
    wc = 2*h*sin(2*pi*f0*t);                        % amplitude 2 w/h
    [K,W] = run_case(t, wc, h, C);
    chk(K.label==2, 'a pure tone was labelled "%s", expected LCO', K.label_name);
    fest = median(W.fdom(2:end));
    chk(abs(fest-f0) < 5, 'dominant frequency came out %.1f Hz, expected %.0f Hz', fest, f0);
    chk(abs(K.amp_last - 2.0) < 0.05, 'amplitude came out %.3f w/h, expected 2.0', K.amp_last);
    fprintf('    clean tone       -> %s, f = %.1f Hz, amp = %.3f w/h\n', K.label_name, fest, K.amp_last);

    % ---------- 3. many incommensurate tones are broadband ----------------
    rng_reset();
    wc = zeros(size(t));
    fs = f0*[0.61 0.83 1.00 1.27 1.53 1.81 2.13 2.44 2.91 3.37];
    for k = 1:numel(fs)
        wc = wc + (2*h/numel(fs))*sin(2*pi*fs(k)*t + 1.7*k);
    end
    [K,~] = run_case(t, wc, h, C);
    chk(K.label==3, 'a multi-tone signal was labelled "%s", expected broadband', K.label_name);
    fprintf('    multi-tone       -> %s\n', K.label_name);

    % ---------- 4. a mid-run transition, at a known time ------------------
    % LCO for the first 4 s, then decaying to static. This is the case-1
    % behaviour in miniature: a 3 s window would have called it an LCO.
    tt = 4.0;
    env = ones(size(t));
    env(t>tt) = exp(-(t(t>tt)-tt)/0.15);
    wc = 2*h*env.*sin(2*pi*f0*t);
    [K,W] = run_case(t, wc, h, C);
    T = qs_detect_transition(W, C);
    chk(T.flag, 'a clear LCO->static transition at %.1f s was NOT detected', tt);
    if T.flag
        chk(abs(T.time - tt) < 1.0, 'transition time came out %.2f s, expected ~%.1f s', T.time, tt);
        chk(T.from==2 && T.to==1, 'transition reported %s -> %s, expected LCO -> static', T.from_name, T.to_name);
        fprintf('    transition       -> %s -> %s at %.2f s (true %.1f s)\n', T.from_name, T.to_name, T.time, tt);
    end
    chk(K.flags.nonstationary, 'a record containing a transition was not flagged nonstationary');
    chk(K.label==1, 'the final label after the transition is "%s", expected static', K.label_name);

    % a 3 s record of the same signal must NOT report a transition --
    % this is the false-positive guard
    m3 = t <= 3.0;
    [K3,W3] = run_case(t(m3), wc(m3), h, C);
    T3 = qs_detect_transition(W3, C);
    chk(~T3.flag, 'a transition was invented inside the first 3 s, before it happens');
    chk(K3.label==2, 'the first 3 s should still read as LCO, got "%s"', K3.label_name);
    fprintf('    same signal, 3 s -> %s, no transition (correct)\n', K3.label_name);

    % ---------- 5. growth is divergent ------------------------------------
    wc = 0.05*h*exp(t/1.2).*sin(2*pi*f0*t);
    [K,~] = run_case(t, wc, h, C);
    chk(K.flags.divergent, 'an exponentially growing signal was not flagged divergent (label "%s")', K.label_name);
    fprintf('    growing          -> %s (divergent flag %d)\n', K.label_name, K.flags.divergent);

    % ---------- 6. near-threshold detection --------------------------------
    % peak-to-peak / 2h = 0.12, just above the 0.10 tolerance and inside the
    % +/-50% band, so the static/LCO call here could flip with a longer run
    wc = 0.12*h*sin(2*pi*f0*t);
    [K,~] = run_case(t, wc, h, C);
    chk(K.flags.near_threshold, ...
        'amplitude %.3f w/h against a tolerance of %.2f was not flagged near_threshold', ...
        K.amp_last, C.classify.amp_tol_wh);
    fprintf('    near threshold   -> amp %.3f w/h, flagged %d\n', K.amp_last, K.flags.near_threshold);
end

% ======================================================================
function [K, W] = run_case(t, wc, h, C)
    W = qs_windows(t, wc, h, C);
    m = struct('Amp_wh',(max(wc)-min(wc))/(2*h), 'wc',wc);
    K = qs_classify(m, W, C);
end

function rng_reset()
    try, rng(0); catch, end
end
