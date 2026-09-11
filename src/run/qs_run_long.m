function qs_run_long(caseId, studyId, which)
%QS_RUN_LONG  Long single-point run(s) through the same code path as the sweep.
%
%   qs_run_long(caseId, studyId)        run every point in the long study
%   qs_run_long(caseId, studyId, k)     run only point k of that list
%
%   A "long" study is just a point list with a large t_end and save_hist on.
%   The integration, the metrics, the window signature and the classification
%   are computed by exactly the same functions the sweep uses (qs_run_point ->
%   qs_metrics / qs_windows / qs_classify), so a 20 s nominal run and a 3 s
%   sweep cell are the same measurement at different durations and can be
%   plotted on the same map.
%
%   What this adds on top of qs_run_point:
%     * modal history q, qdot saved in full             (postproc_data.mat)
%     * deflection at the DIC measurement locations     (deformation_at_DIC.mat)
%     * optional instantaneous field streaming          (fields/)
%     * the DIC comparison figures
%
%   Field streaming is controlled by the study's "fields" block:
%       "fields": {"save": true, "mode": "per_second", "precision": "single",
%                  "stride": 1}
%   Mode 'per_second' writes one file per simulated second, which is what
%   makes a 20 s run recoverable if it dies part way through.

    if nargin < 3, which = []; end
    C = qs_config(caseId, studyId);
    assert(exist(C.paths.pointlist,'file')==2, ...
        'qs_run_long: run qs_init_study(''%s'',''%s'') first', caseId, studyId);

    P = qs_pointlist('load', C.paths.pointlist);
    if isempty(which), which = 1:numel(P.id); end

    M = qs_build_model(C);
    qs_preflight(M, C);

    for k = which(:).'
        pt = qs_pointlist('rows', P, k);
        outdir = fullfile(C.paths.results_dir, pt.id{1});
        if exist(outdir,'dir')~=7, mkdir(outdir); end

        fprintf('=== long run %d/%d : p_c = %.2f kPa, dT = %.3f K, %.1f s ===\n', ...
            k, numel(P.id), pt.pc_Pa(1)/1e3, pt.dT_K(1), pt.t_end(1));

        % The heavy lifting: identical to a sweep point, with the history kept.
        t0 = tic;
        [R, extra] = run_with_state(M, C, pt);
        fprintf('  integrated in %.1f min -> %s\n', toc(t0)/60, R.label_name);

        if ~R.ok
            warning('qs_run_long:failed','point %s failed: %s', R.id, R.errmsg);
            continue
        end

        report_point(R);

        % ---- modal history (saved FIRST: never lose an expensive run to a
        %      plotting error) ------------------------------------------------
        tvec = extra.tvec;  q = extra.q;  qdot = extra.qdot; %#ok<NASGU>
        metrics = R;  metrics.W = R.W; %#ok<NASGU>
        save(fullfile(outdir,'postproc_data.mat'), ...
             'tvec','q','qdot','metrics',qs_matver());

        % ---- deflection at the DIC locations -------------------------------
        if M.dic.present
            w_dic = q * M.dic.Phi_dic.';        %#ok<NASGU>  Nt x Ndic [m]
            X_dic = M.dic.X;  Y_dic = M.dic.Y;  %#ok<NASGU>
            save(fullfile(outdir,'deformation_at_DIC.mat'), ...
                 'tvec','w_dic','X_dic','Y_dic',qs_matver());
        end

        % ---- instantaneous fields ------------------------------------------
        Fb = getdef(C.study,'fields', struct('save',false));
        if getdef(Fb,'save',false)
            fdir = fullfile(outdir,'fields');
            if exist(fdir,'dir')~=7, mkdir(fdir); end
            qs_write_fields(fdir, M, C, tvec, q, qdot, Fb);
            fprintf('  fields -> %s\n', fdir);
        end

        % ---- figures (non-fatal) --------------------------------------------
        try
            qs_plot_history(C, R, outdir);
        catch ME
            warning('qs_run_long:figures','figures failed (%s); data is saved.', ME.message);
        end
        fprintf('  done -> %s\n\n', outdir);
    end
end

% ======================================================================
function [R, extra] = run_with_state(M, C, pt)
%RUN_WITH_STATE  qs_run_point, but keeping the full modal history.
%   qs_run_point deliberately throws the modal history away (2000 points x
%   200k steps x 18 modes would be unmanageable). A long run wants it, so the
%   integration is repeated here with the same settings and the same
%   reduction functions -- not a different model.
    extra = struct('tvec',[],'q',[],'qdot',[]);
    aero = M.aero;  aero.pback = pt.pc_Pa(1);
    tvec = 0 : C.duration.dt_out : pt.t_end(1);

    cfg = struct('dT', pt.dT_K(1), 'tspan', tvec, 'integrator', @ode15s, ...
                 'odeopts', odeset('RelTol',C.solver.RelTol,'AbsTol',C.solver.AbsTol), ...
                 'reconstruct', false, 'staticPressure', []);
    cfg.pressureModel = @(s) pm_quasisteady_ept(s, aero);
    cfg.y0 = qs_initial_state(M, C, pt.ic_tag{1});

    R = [];
    try
        out = rc19_solve_structural(M.rom, cfg);
        extra.tvec = tvec(:);  extra.q = out.q;  extra.qdot = out.qdot;

        m = qs_metrics(M, C, tvec, out.q, pt.t_trans(1));
        W = qs_windows(tvec, m.wc, M.h, C);
        K = qs_classify(m, W, C);
        T = qs_detect_transition(W, C);

        R = pack(pt, m, K, T, W, tvec);
        R.ok = true;  R.errmsg = '';
    catch ME
        R = pack(pt, [], [], [], [], []);
        R.ok = false;
        R.errmsg = sprintf('%s | %s', ME.identifier, ME.message);
    end
end

function R = pack(pt, m, K, T, W, tvec)
    R = struct('id',pt.id{1}, 'pc_Pa',pt.pc_Pa(1), 'dT_K',pt.dT_K(1), ...
               't_end',pt.t_end(1), 't_trans',pt.t_trans(1), 'ic_tag',pt.ic_tag{1}, ...
               'Amp_wh',NaN,'MeanPeak',NaN,'StdPeak',NaN,'MeanRMSE',NaN,'StdRMSE',NaN, ...
               'Freq',NaN,'pfrac',NaN,'signflip',1,'label',0,'label_name','unknown', ...
               'label_start',0,'nonstationary',false,'divergent',false, ...
               'near_threshold',false,'short_record',false,'amp_first',NaN,'amp_last',NaN, ...
               'trans_flag',false,'trans_time',NaN,'trans_t_lo',NaN,'trans_t_hi',NaN, ...
               'trans_from',0,'trans_to',0,'n_changes',0, ...
               'W',[],'wc',[],'tvec',[],'ok',false,'errmsg','');
    if isempty(m), return; end
    R.Amp_wh=m.Amp_wh; R.MeanPeak=m.MeanPeak; R.StdPeak=m.StdPeak;
    R.MeanRMSE=m.MeanRMSE; R.StdRMSE=m.StdRMSE; R.Freq=m.Freq; R.pfrac=m.pfrac;
    R.signflip=m.signflip;
    R.label=K.label; R.label_name=K.label_name; R.label_start=K.label_start;
    R.nonstationary=K.flags.nonstationary; R.divergent=K.flags.divergent;
    R.near_threshold=K.flags.near_threshold; R.short_record=K.flags.short_record;
    R.amp_first=K.amp_first; R.amp_last=K.amp_last;
    R.trans_flag=T.flag; R.trans_time=T.time; R.trans_t_lo=T.time_lo; R.trans_t_hi=T.time_hi;
    R.trans_from=T.from; R.trans_to=T.to;
    R.n_changes=T.n_changes;
    R.W = W;  R.wc = m.wc(:);  R.tvec = tvec(:);
end

function report_point(R)
    fprintf('  amplitude %.3f w/h   f = %.0f Hz   peak-power fraction %.2f\n', ...
        R.Amp_wh, R.Freq, R.pfrac);
    fprintf('  mean peak %.3f mm   std peak %.3f mm   meanRMSE %.3f   stdRMSE %.3f\n', ...
        R.MeanPeak, R.StdPeak, R.MeanRMSE, R.StdRMSE);
    if R.trans_flag
        fprintf('  *** REGIME CHANGE: %s -> %s at t = %.3f s ***\n', ...
            qs_label_name(R.trans_from), qs_label_name(R.trans_to), R.trans_time);
    elseif R.nonstationary
        fprintf('  NOTE: nonstationary -- amplitude still drifting at t_end. Run longer.\n');
    end
end

function v = getdef(S, f, d)
    if isstruct(S) && isfield(S,f) && ~isempty(S.(f)), v = S.(f); else, v = d; end
end
