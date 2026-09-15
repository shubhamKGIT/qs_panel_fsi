function A = t03b_attribute(pc_kPa, dT_K, refpath, olddir)
%T03B_ATTRIBUTE  Settle WHERE a regression difference comes from, for one point.
%
%   A = t03b_attribute()                 the default point: 53.20 kPa, 14.050 K
%   A = t03b_attribute(pc_kPa, dT_K)     any point on the original grid
%   A = t03b_attribute(pc, dT, refpath, olddir)
%
%   t03 compares SUMMARY NUMBERS. When they disagree, a summary cannot tell you
%   whether the cause is a setting that drifted, or two runs that started
%   identically and separated later because the point sits near a boundary.
%   This test answers that, by comparing the actual TIME HISTORIES.
%
%   It runs the same operating point three ways:
%     (1) the centre-deflection history STORED in the original map (W_center)
%     (2) the ORIGINAL driver code, re-run now, in this MATLAB    [if available]
%     (3) the new framework
%
%   and then reports WHEN each pair separates. The interpretation is simple:
%
%     separate at t = 0, or within a few steps
%         -> a setting differs. Something in case.json or defaults.json is not
%            what the old driver hardcoded. This is a real bug: hunt it.
%
%     identical for a while, then diverge
%         -> the two runs are the same computation and the point is sensitive.
%            The integrator followed a slightly different path (different MATLAB
%            release, different machine) and the trajectories parted. Look at
%            whether the point had even settled: if the window labels are still
%            changing at t_end, neither answer is converged and the point is a
%            poor regression target in the first place.
%
%     (2) agrees with (3) but neither agrees with (1)
%         -> the stored map was produced by a different MATLAB or machine. The
%            framework is reproducing the old CODE correctly; it is the old
%            NUMBERS that cannot be reproduced here.

    if nargin < 1 || isempty(pc_kPa), pc_kPa = 53.20;  end
    if nargin < 2 || isempty(dT_K),   dT_K   = 14.050; end
    if nargin < 3, refpath = ''; end
    if nargin < 4, olddir  = ''; end

    root = qs_root();
    if isempty(refpath), refpath = find_first({ ...
        fullfile(root,'..','Unheated_FSI_Analysis','case1_run','stability_sweep','case1_NoSBLI_Periodic','detailed_stability_map_data.mat'), ...
        fullfile(root,'tests','ref','detailed_stability_map_data.mat'), ...
        '/Volumes/my_80Gb_box/FTSI/Unheated_FSI_Analysis/case1_run/stability_sweep/case1_NoSBLI_Periodic/detailed_stability_map_data.mat'}); end
    % The pristine launch templates moved out of Unheated_FSI_Analysis into
    % FTSI/Legacy_script_template/ (2026-09-15). Both locations are listed so
    % this test works either side of that move; the new one is tried first.
    if isempty(olddir),  olddir  = find_first_dir({ ...
        fullfile(root,'..','Legacy_script_template','stability_sweep','case1_NoSBLI_Periodic'), ...
        fullfile(root,'..','Unheated_FSI_Analysis','stability_sweep','case1_NoSBLI_Periodic'), ...
        fullfile(root,'..','Unheated_FSI_Analysis','case1_run','stability_sweep','case1_NoSBLI_Periodic'), ...
        '/Volumes/my_80Gb_box/FTSI/Legacy_script_template/stability_sweep/case1_NoSBLI_Periodic', ...
        '/Volumes/my_80Gb_box/FTSI/Unheated_FSI_Analysis/stability_sweep/case1_NoSBLI_Periodic'}); end

    assert(~isempty(refpath), 't03b: original map not found -- pass its path as the 3rd argument');
    fprintf('\n================ t03b: attributing one point ================\n');
    fprintf('  point     : p_c = %.2f kPa,  dT = %.3f K\n', pc_kPa, dT_K);
    fprintf('  reference : %s\n', refpath);
    fprintf('  old code  : %s\n\n', tern(isempty(olddir),'(not found -- step 2 skipped)', olddir));

    C = qs_config('unheated_c1_NoSBLI_Periodic','regression_unheated_c1');
    h = C.case.panel.h_m;

    % ---------- (1) the stored history --------------------------------------
    mf  = matfile(refpath);
    pcv = mf.pc_vec;  dTv = mf.dT_vec;
    ip  = nearest(pcv, pc_kPa*1e3, 1);      % 1 Pa
    jd  = nearest(dTv, dT_K,      1e-6);
    assert(~isempty(ip) && ~isempty(jd), ...
        't03b: (%.2f kPa, %.3f K) is not a point on the original grid', pc_kPa, dT_K);

    wc_ref = double(squeeze(mf.W_center(ip, jd, :)));
    tvec   = double(mf.tvec).';
    stored = struct('Amp_wh', mf.Amp_wh(ip,jd), 'Freq', mf.Freq(ip,jd), ...
                    'MeanPeak', mf.MeanPeak(ip,jd), 'StdPeak', mf.StdPeak(ip,jd), ...
                    'MeanRMSE', mf.MeanRMSE(ip,jd), 'StdRMSE', mf.StdRMSE(ip,jd));
    pc_Pa = pcv(ip);  dT = dTv(jd);
    fprintf('  grid cell : index (%d,%d) -> p_c = %.1f Pa, dT = %.4f K\n\n', ip, jd, pc_Pa, dT);

    % ---------- (3) the framework ---------------------------------------------
    % Double precision for the history: this whole test is about WHERE two runs
    % separate, and the default single-precision storage would put a floor of
    % about 1e-7 w/h on that, hiding the difference between "identical" and
    % "identical to storage precision".
    C.store.history_precision = 'double';
    fprintf('  [framework] running ... ');
    M  = qs_build_model(C, true);
    pt = qs_pointlist('make', pc_Pa, dT, struct( ...
            't_end', C.duration.base_t_end, 't_trans', C.duration.base_t_transient, ...
            'ic_tag', 'flat'));
    t0 = tic;
    Rn = qs_run_point(M, C, pt, true);
    fprintf('%.1f min\n', toc(t0)/60);
    wc_new = double(Rn.wc(:));

    % ---------- (2) the original driver code, re-run now ----------------------
    wc_old = [];
    if ~isempty(olddir)
        fprintf('  [old code ] running ... ');
        t0 = tic;
        wc_old = run_old_driver_point(olddir, pc_Pa, dT, C);
        fprintf('%.1f min\n', toc(t0)/60);
    end
    fprintf('\n');

    % ---------- scalar comparison ----------------------------------------------
    fprintf('  --- summary numbers -------------------------------------------\n');
    fprintf('  %-10s %16s %16s\n', '', 'stored map', 'framework now');
    nm = {'Amp_wh','Freq','MeanPeak','StdPeak','MeanRMSE','StdRMSE'};
    got = {Rn.Amp_wh, Rn.Freq, Rn.MeanPeak, Rn.StdPeak, Rn.MeanRMSE, Rn.StdRMSE};
    for k = 1:numel(nm)
        fprintf('  %-10s %16.8g %16.8g\n', nm{k}, stored.(nm{k}), got{k});
    end
    fprintf('\n');

    % ---------- where do the histories separate? --------------------------------
    A = struct('pc_Pa',pc_Pa,'dT_K',dT,'stored',stored,'new',Rn);
    fprintf('  --- where the time histories separate -------------------------\n');
    A.sep_ref_new = separation(tvec, wc_ref, wc_new, h, 'stored map', 'framework');
    if ~isempty(wc_old)
        A.sep_ref_old = separation(tvec, wc_ref, wc_old, h, 'stored map', 'old code now');
        A.sep_old_new = separation(tvec, wc_old, wc_new, h, 'old code now', 'framework');
    end
    fprintf('\n');

    % ---------- had the point even settled? --------------------------------------
    fprintf('  --- had this point settled by t_end? ---------------------------\n');
    A.W_ref = report_windows(tvec, wc_ref, h, C, 'stored map');
    A.W_new = report_windows(tvec, wc_new, h, C, 'framework  ');
    if ~isempty(wc_old)
        A.W_old = report_windows(tvec, wc_old, h, C, 'old code   ');
    end

    % ---------- verdict ------------------------------------------------------------
    fprintf('\n  --- reading ----------------------------------------------------\n');
    % The verdict keys on MAXIMUM difference over the record, and on the 1e-6
    % crossing. Those are robust. The 1e-12 and 1e-9 crossings are reported for
    % information but are not used to decide anything: the stored history is
    % single precision, so nothing finer than about 1e-7 w/h can be resolved
    % against it and those two lines can read 0.0000 s for runs that are in
    % fact identical.
    SAME = 1e-5;                          % max diff below this = the same run
    if ~isempty(wc_old) && A.sep_old_new.max < SAME && A.sep_ref_new.max > 1e-2
        fprintf(['  THE FRAMEWORK REPRODUCES THE OLD CODE.\n' ...
                 '  Old code re-run now vs framework: max difference %.2e w/h -- the\n' ...
                 '  same run to storage precision. Both then disagree with the stored\n' ...
                 '  map in the SAME way (max %.2e w/h, first exceeding 1e-6 at the same\n' ...
                 '  time). So the stored numbers came from a different MATLAB or\n' ...
                 '  machine, not from different code, and this point amplifies that\n' ...
                 '  difference because it is sensitive -- see the window labels above.\n'], ...
                 A.sep_old_new.max, A.sep_ref_new.max);
    elseif A.sep_ref_new.t_1em6 <= 5*C.duration.dt_out
        fprintf(['  The framework and the stored run differ from the FIRST STEPS.\n' ...
                 '  That is a settings difference, not sensitivity. Compare case.json\n' ...
                 '  against the old driver line by line.\n']);
    elseif ~isempty(wc_old) && A.sep_old_new.max > SAME
        fprintf(['  The old code re-run now and the framework DIFFER from each other\n' ...
                 '  (max %.2e w/h). Something in the plumbing is not equivalent to the\n' ...
                 '  old driver. This is the case to investigate.\n'], A.sep_old_new.max);
    else
        fprintf(['  The runs start together and separate later. That is trajectory\n' ...
                 '  sensitivity, not a settings difference. Check the window labels\n' ...
                 '  above: if they are still changing at t_end, neither run has\n' ...
                 '  converged and this point should not be used as a regression target.\n']);
    end
    fprintf('===============================================================\n\n');
end

% ======================================================================
function S = separation(t, a, b, h, na, nb)
%SEPARATION  First time two centre histories differ by more than a threshold,
%   measured in panel thicknesses so the numbers mean something physical.
    n = min([numel(a), numel(b), numel(t)]);
    dd = abs(a(1:n) - b(1:n)) / h;
    S = struct();
    for thr = [1e-12 1e-9 1e-6 1e-3 1e-1]
        k = find(dd > thr, 1, 'first');
        fld = sprintf('t_%s', strrep(sprintf('1em%d', round(-log10(thr))), '-', 'm'));
        if isempty(k), S.(fld) = Inf; else, S.(fld) = t(k); end
    end
    S.max = max(dd);
    fprintf('  %s vs %s :\n', na, nb);
    fprintf('     exceeds 1e-12 w/h at t = %s\n', fmt(S.t_1em12));
    fprintf('     exceeds 1e-9  w/h at t = %s   (below storage precision -- ignore)\n', fmt(S.t_1em9));
    fprintf('     exceeds 1e-6  w/h at t = %s\n', fmt(S.t_1em6));
    fprintf('     exceeds 1e-3  w/h at t = %s\n', fmt(S.t_1em3));
    fprintf('     MAX difference over the record = %.3e w/h   <-- the number to read\n', S.max);
end

function W = report_windows(t, wc, h, C, tag)
    W = qs_windows(t(:), wc(:), h, C);
    m = struct('Amp_wh',(max(wc)-min(wc))/(2*h), 'wc', wc);
    K = qs_classify(m, W, C);
    T = qs_detect_transition(W, C);
    lab = qs_label_name(W.label);
    if ~iscell(lab), lab = {lab}; end
    fprintf('  %s : %-11s  amp first half %.4g -> last half %.4g w/h\n', ...
        tag, K.label_name, K.amp_first, K.amp_last);
    fprintf('  %s   window labels: %s\n', repmat(' ',1,numel(tag)), strjoin(lab, ' '));
    if K.flags.nonstationary
        fprintf('  %s   NOT STATIONARY -- this run had not settled by t_end\n', repmat(' ',1,numel(tag)));
    end
    if T.flag
        fprintf('  %s   regime change %s -> %s at t = %.2f s\n', ...
            repmat(' ',1,numel(tag)), T.from_name, T.to_name, T.time);
    end
end

function wc = run_old_driver_point(olddir, pc_Pa, dT, C)
%RUN_OLD_DRIVER_POINT  The ORIGINAL driver's per-point body, verbatim.
%   The old folder is put FIRST on the path so its own copies of the solver and
%   its own rc19_aero_data.mat are used -- cm_setup_aero defaults the aero file
%   to whatever sits next to itself, which is how the old layout worked. The
%   path is restored afterwards whatever happens.
    saved = path();
    cleanup = onCleanup(@() path(saved));
    addpath(olddir, '-begin');

    L = 0.254;
    rom  = rc19_setup_rom(fullfile(olddir,'RC19_ROM_v1.mat'));
    opts = struct('pinf',50287,'Minf',1.92,'gamma',1.4,'Lp',L,'pback',50429.27,'zsign',-1);
    [~, aero] = evalc('cm_setup_aero(rom, opts)');
    nm = rom.num_modes;
    cidx = round(rom.nx/2) + (round(rom.ny/2)-1)*rom.nx;

    t_end = 3.0;  dt_out = 1e-4;
    tvec  = 0:dt_out:t_end;

    cfg = struct('tspan',tvec,'integrator',@ode15s, ...
                 'odeopts',odeset('RelTol',1e-7,'AbsTol',1e-9), ...
                 'reconstruct',false,'staticPressure',[]);
    cfg.y0 = zeros(2*nm,1);  cfg.y0(1:nm) = 1e-4;

    aero.pback = pc_Pa;  cfg.dT = dT;
    cfg.pressureModel = @(s) pm_quasisteady_ept(s, aero);
    out = rc19_solve_structural(rom, cfg);
    wc  = out.q * rom.Modes(:,cidx);
    wc  = double(wc(:));
end

% ---------------------------------------------------------------- helpers
function i = nearest(v, x, tol)
    [dmin, i] = min(abs(v(:) - x));
    if dmin > tol, i = []; end
end
function s = fmt(t)
    if isinf(t), s = 'never (identical to this level)'; else, s = sprintf('%.4f s', t); end
end
function f = find_first(c)
    f = '';
    for k = 1:numel(c), if exist(c{k},'file')==2, f = c{k}; return; end, end
end
function f = find_first_dir(c)
    f = '';
    for k = 1:numel(c), if exist(c{k},'dir')==7, f = c{k}; return; end, end
end
function s = tern(c,a,b)
    if c, s = a; else, s = b; end
end
