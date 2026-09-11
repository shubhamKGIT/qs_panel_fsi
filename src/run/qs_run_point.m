function R = qs_run_point(M, C, pt, want_history)
%QS_RUN_POINT  Integrate ONE operating point and reduce it to a result record.
%
%   R = QS_RUN_POINT(M, C, pt)
%   R = QS_RUN_POINT(M, C, pt, want_history)
%
%     M   model from qs_build_model (built once, reused for every point)
%     C   resolved config from qs_config
%     pt  one row of a point list, e.g. qs_pointlist('rows', P, k)
%     want_history  force keeping the full centre history (default: pt.save_hist)
%
%   This single function is the ONLY place a coupled run is launched. The
%   sweep calls it in a loop; the long 20 s run calls it once. That is what
%   guarantees a 20 s nominal run and a sweep cell are the same computation
%   with a different t_end, so one can be plotted on top of the other.
%
%   OUTPUT (struct R) -- one row of the results table
%     .id .pc_Pa .dT_K .t_end .t_trans .level .ic_tag
%     .Amp_wh .MeanPeak .StdPeak .MeanRMSE .StdRMSE .Freq .pfrac .signflip
%     .label .label_name .label_start
%     .nonstationary .divergent .near_threshold .short_record
%     .trans_flag .trans_time .trans_from .trans_to .n_changes
%     .wall_s .ok .errmsg
%     .W        the window signature table (always kept: a few hundred bytes)
%     .wc       full centre history, single precision (ONLY if requested/flagged)
%     .tvec     its time vector (only alongside .wc)
%     .yend     final state [q;qdot] (only if C.store.save_final_state)

    if nargin < 4, want_history = []; end

    id    = pt.id{1};
    pc    = pt.pc_Pa(1);
    dT    = pt.dT_K(1);
    tend  = pt.t_end(1);
    ttr   = pt.t_trans(1);
    ictag = pt.ic_tag{1};

    R = blank_record(pt);
    t0 = tic;

    try
        % ---- per-point operating condition -------------------------------
        % aero is a local copy inside this function, so overwriting pback
        % cannot leak into the next point.
        aero = M.aero;
        aero.pback = pc;

        tvec = 0 : C.duration.dt_out : tend;

        cfg = struct();
        cfg.dT             = dT;
        cfg.tspan          = tvec;
        cfg.integrator     = solver_handle(C.solver.integrator);
        cfg.odeopts        = odeset('RelTol', C.solver.RelTol, 'AbsTol', C.solver.AbsTol);
        cfg.reconstruct    = false;
        cfg.staticPressure = [];
        cfg.pressureModel  = @(s) pm_quasisteady_ept(s, aero);
        cfg.y0             = qs_initial_state(M, C, ictag);

        % ---- integrate ----------------------------------------------------
        out = rc19_solve_structural(M.rom, cfg);

        % ---- reduce --------------------------------------------------------
        m = qs_metrics(M, C, tvec, out.q, ttr);
        W = qs_windows(tvec, m.wc, M.h, C);
        K = qs_classify(m, W, C);
        T = qs_detect_transition(W, C);

        R.Amp_wh   = m.Amp_wh;    R.MeanPeak = m.MeanPeak;  R.StdPeak = m.StdPeak;
        R.MeanRMSE = m.MeanRMSE;  R.StdRMSE  = m.StdRMSE;
        R.Freq     = m.Freq;      R.pfrac    = m.pfrac;     R.signflip = m.signflip;

        R.label       = K.label;        R.label_name  = K.label_name;
        R.label_start = K.label_start;
        R.nonstationary  = K.flags.nonstationary;
        R.divergent      = K.flags.divergent;
        R.near_threshold = K.flags.near_threshold;
        R.short_record   = K.flags.short_record;
        R.amp_first = K.amp_first;  R.amp_last = K.amp_last;

        R.trans_flag = T.flag;   R.trans_time = T.time;
        R.trans_t_lo = T.time_lo; R.trans_t_hi = T.time_hi;
        R.trans_from = T.from;   R.trans_to   = T.to;
        R.n_changes  = T.n_changes;

        R.W  = compact_windows(W);
        R.ok = true;

        % ---- keep the full history only when it earns its place -----------
        if isempty(want_history)
            want_history = pt.save_hist(1) || history_wanted(R, C);
        end
        if want_history
            R.wc   = single(m.wc(:).');
            R.tvec = single(tvec(:).');
        end

        if isfield(C.store,'save_final_state') && C.store.save_final_state
            R.yend = [out.q(end,:).'; out.qdot(end,:).'];
        end

    catch ME
        R.ok     = false;
        R.errmsg = sprintf('%s | %s', ME.identifier, ME.message);
        warning('qs_run_point:failed', 'point %s (p_c=%.1f kPa, dT=%.3f K) FAILED: %s', ...
            id, pc/1e3, dT, ME.message);
    end

    R.wall_s = toc(t0);
end

% ======================================================================
function tf = history_wanted(R, C)
%HISTORY_WANTED  Storage policy: which points keep their full centre history.
%   A 20 s record at dt_out = 1e-4 is 800 kB in single precision, so keeping
%   every point of a 2000-point sweep would be ~1.6 GB. The window signature
%   (a few hundred bytes) is always kept; the full trace is kept only for the
%   points that actually show something.
    tf = false;
    if ~isfield(C.store,'save_history_if'), return; end
    want = C.store.save_history_if;
    if ischar(want), want = {want}; end
    for k = 1:numel(want)
        switch lower(want{k})
            case 'transition',     tf = tf || R.trans_flag;
            case 'nonstationary',  tf = tf || R.nonstationary;
            case 'near_threshold', tf = tf || R.near_threshold;
            case 'divergent',      tf = tf || R.divergent;
            case 'all',            tf = true;
            case 'none'            % explicit no-op
            otherwise
                warning('qs_run_point:badPolicy','unknown store.save_history_if entry "%s"', want{k});
        end
    end
end

function Wc = compact_windows(W)
%COMPACT_WINDOWS  Shrink the window table for storage.
    Wc.tmid   = single(W.tmid(:).');
    Wc.amp_wh = single(W.amp_wh(:).');
    Wc.wmean  = single(W.wmean(:).');
    Wc.fdom   = single(W.fdom(:).');
    Wc.pfrac  = single(W.pfrac(:).');
    Wc.label  = int8(W.label(:).');
    Wc.T_win_s = W.T_win_s;
    Wc.overlap = W.overlap;
end

function h = solver_handle(nameStr)
    switch lower(nameStr)
        case 'ode15s', h = @ode15s;
        case 'ode45',  h = @ode45;
        case 'ode23s', h = @ode23s;
        case 'ode113', h = @ode113;
        otherwise, error('qs_run_point: unsupported solver.integrator "%s"', nameStr);
    end
end

function R = blank_record(pt)
    R.id      = pt.id{1};
    R.pc_Pa   = pt.pc_Pa(1);
    R.dT_K    = pt.dT_K(1);
    R.t_end   = pt.t_end(1);
    R.t_trans = pt.t_trans(1);
    R.level   = pt.level(1);
    R.ic_tag  = pt.ic_tag{1};

    R.Amp_wh = NaN; R.MeanPeak = NaN; R.StdPeak = NaN;
    R.MeanRMSE = NaN; R.StdRMSE = NaN; R.Freq = NaN; R.pfrac = NaN; R.signflip = 1;
    R.label = 0; R.label_name = 'unknown'; R.label_start = 0;
    R.nonstationary = false; R.divergent = false;
    R.near_threshold = false; R.short_record = false;
    R.amp_first = NaN; R.amp_last = NaN;
    R.trans_flag = false; R.trans_time = NaN;
    R.trans_t_lo = NaN; R.trans_t_hi = NaN;
    R.trans_from = 0; R.trans_to = 0; R.n_changes = 0;
    R.W = []; R.wc = []; R.tvec = []; R.yend = [];
    R.wall_s = NaN; R.ok = false; R.errmsg = '';
end
