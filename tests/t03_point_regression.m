function R = t03_point_regression(npoints, refpath)
%T03  Reproduce points from the ORIGINAL case-1 sweep through the new code.
%
%   R = t03_point_regression(npoints, refpath)
%     npoints  how many grid points to re-run (default 4, ~3 min each)
%     refpath  path to the old detailed_stability_map_data.mat
%              (default: searched next to the framework)
%
%   THIS IS THE TEST THAT MATTERS. The six core files moved into src/core and
%   src/aero byte-for-byte, so if the configuration plumbing did not change a
%   number anywhere, every metric must come back identical to the old map.
%
%   TWO THINGS MAKE A NAIVE COMPARISON USELESS, AND BOTH ARE HANDLED HERE.
%
%   1. Absolute floors. Most of the map is a panel sitting still, where the
%      amplitude is 1e-10 panel thicknesses and the peak deflection is 1e-8 mm
%      -- numerical dust. Asking two such numbers to agree to one part in a
%      million is asking roundoff to be reproducible, which it is not. Each
%      metric therefore has a floor below which both values are reported as
%      negligible rather than compared. The floors are set where the quantity
%      stops being physically meaningful, not where it becomes convenient.
%
%   2. Frequency of a record with no oscillation. The dominant frequency is the
%      argmax of a spectrum. With no oscillation in the record that is the
%      argmax of noise and can land anywhere. It is compared only when the
%      motion is large enough for a frequency to mean something.
%
%   A point that genuinely differs is NOT necessarily a bug -- it may sit near a
%   regime boundary where two runs that started identically separated later. Use
%   t03b_attribute(pc_kPa, dT_K) on any failing point: it compares the time
%   HISTORIES, which is what distinguishes a settings difference from
%   trajectory sensitivity.

    if nargin < 1 || isempty(npoints), npoints = 4; end
    if nargin < 2, refpath = ''; end

    R = struct('pass',true,'skip',false,'msg',{{}});
    function chk(cond, fmt, varargin)
        if ~cond
            R.pass = false;
            R.msg{end+1} = sprintf(fmt, varargin{:});
            fprintf('    FAIL  %s\n', R.msg{end});
        end
    end

    % ---------- comparison policy -------------------------------------------
    TOLREL = 1e-6;
    FLOOR = struct( ...
        'Amp_wh',   1e-3,  ...   % w/h. A thousandth of a panel thickness is not motion.
        'MeanPeak', 1e-4,  ...   % mm.  0.1 micron; DIC noise is larger than this.
        'StdPeak',  1e-4,  ...   % mm.
        'MeanRMSE', 1e-9,  ...   % ratio, O(1); effectively no floor
        'StdRMSE',  1e-9,  ...
        'Freq',     0);          % handled by the amplitude gate below
    FREQ_NEEDS_AMP = 1e-2;       % w/h below which a "dominant frequency" is noise

    % ---------- locate the reference map ------------------------------------
    if isempty(refpath)
        root = qs_root();
        cands = { ...
            fullfile(root,'..','Unheated_FSI_Analysis','case1_run','stability_sweep','case1_NoSBLI_Periodic','detailed_stability_map_data.mat'), ...
            fullfile(root,'tests','ref','detailed_stability_map_data.mat'), ...
            '/Volumes/my_80Gb_box/FTSI/Unheated_FSI_Analysis/case1_run/stability_sweep/case1_NoSBLI_Periodic/detailed_stability_map_data.mat'};
        for k = 1:numel(cands)
            if exist(cands{k},'file')==2, refpath = cands{k}; break; end
        end
    end
    if isempty(refpath) || exist(refpath,'file')~=2
        R.skip = true;
        R.msg{end+1} = 'reference map not found; pass its path: t03_point_regression(4, ''/path/to/detailed_stability_map_data.mat'')';
        fprintf('    SKIP: %s\n', R.msg{end});
        return
    end
    fprintf('    reference: %s\n', refpath);

    Ref = load(refpath, 'pc_vec','dT_vec','Amp_wh','Freq','MeanRMSE','StdRMSE', ...
                        'MeanPeak','StdPeak','t_end','t_transient','dt_out','tvec');

    C = qs_config('unheated_c1_NoSBLI_Periodic','regression_unheated_c1');
    chk(abs(C.duration.base_t_end       - Ref.t_end)       < 1e-12, ...
        'study t_end %.4f does not match the reference run %.4f', C.duration.base_t_end, Ref.t_end);
    chk(abs(C.duration.base_t_transient - Ref.t_transient) < 1e-12, ...
        'study t_transient %.4f does not match the reference %.4f', C.duration.base_t_transient, Ref.t_transient);
    chk(abs(C.duration.dt_out           - Ref.dt_out)      < 1e-20, ...
        'study dt_out %g does not match the reference %g', C.duration.dt_out, Ref.dt_out);

    % ---------- choose points ------------------------------------------------
    %  One point from the static region, to exercise the static-equilibrium path
    %  where the MEAN deflection is well determined even though the amplitude is
    %  numerical dust. The rest spread across the points that actually move,
    %  because those are the ones where every metric carries information.
    valid  = find(isfinite(Ref.Amp_wh) & isfinite(Ref.MeanPeak));
    moving = valid(Ref.Amp_wh(valid) >= FLOOR.Amp_wh);
    static = valid(Ref.Amp_wh(valid) <  FLOOR.Amp_wh);
    assert(~isempty(valid), 't03: the reference map has no finite results');

    sel = [];
    if ~isempty(static)
        [~, o] = sort(Ref.MeanPeak(static));
        sel = static(o(max(1,round(numel(o)/2))));          % a typical static point
    end

    %  A moving point is only a valid regression target if the ORIGINAL run had
    %  settled. Parts of this map are sensitive: two runs that start identically
    %  separate within a tenth of a second and end several panel thicknesses
    %  apart, so their summary numbers cannot be reproduced on another machine —
    %  that is the physics, not the code. Those points are screened out here by
    %  replaying the stored centre history through the same window analysis the
    %  framework uses, and keeping only the runs that were stationary at t_end.
    nmov = max(0, npoints - numel(sel));
    nstat_ok = 0;  ncand = 0;
    if ~isempty(moving) && nmov > 0
        [~, o]   = sort(Ref.Amp_wh(moving));
        cand     = moving(o(unique(round(linspace(1, numel(o), min(40, numel(o)))))));
        ncand    = numel(cand);
        mf       = matfile(refpath);
        tv       = double(Ref.tvec(:));
        keep     = false(size(cand));
        for k = 1:ncand
            [ip, jd] = ind2sub([numel(Ref.pc_vec), numel(Ref.dT_vec)], cand(k));
            try
                wcr = double(squeeze(mf.W_center(ip, jd, :)));
            catch
                keep(k) = true;  continue          % no stored history: allow it
            end
            if any(~isfinite(wcr)), continue; end
            Wr = qs_windows(tv, wcr, C.case.panel.h_m, C);
            Kr = qs_classify(struct('Amp_wh',(max(wcr)-min(wcr))/(2*C.case.panel.h_m),'wc',wcr), Wr, C);
            keep(k) = ~Kr.flags.nonstationary;
        end
        ok = cand(keep);
        nstat_ok = numel(ok);
        if isempty(ok)
            fprintf('    NOTE: none of the %d moving candidates had settled by t_end in the\n', ncand);
            fprintf('          original map. Falling back to the largest-amplitude points;\n');
            fprintf('          expect differences, and use t03b_attribute to confirm why.\n');
            ok = cand;
        end
        [~, o2] = sort(Ref.Amp_wh(ok));
        take = unique(round(linspace(1, numel(o2), min(nmov, numel(o2)))));
        sel  = [sel; ok(o2(take))];
    end
    sel = sel(1:min(numel(sel), npoints));
    nm_sel = sum(Ref.Amp_wh(sel) >= FLOOR.Amp_wh);
    fprintf('    %d point(s): %d static, %d moving-and-settled\n', numel(sel), numel(sel)-nm_sel, nm_sel);
    fprintf('    (map: %d moving of %d valid; %d of %d screened candidates had settled)\n', ...
        numel(moving), numel(valid), nstat_ok, ncand);

    M   = qs_build_model(C, true);
    npc = numel(Ref.pc_vec);
    nMatch = 0;  nNegl = 0;  nDiff = 0;  badPoints = {};

    for a = 1:numel(sel)
        [ip, jd] = ind2sub([npc, numel(Ref.dT_vec)], sel(a));
        pc = Ref.pc_vec(ip);  dT = Ref.dT_vec(jd);

        pt = qs_pointlist('make', pc, dT, struct( ...
                't_end', C.duration.base_t_end, 't_trans', C.duration.base_t_transient, ...
                'ic_tag', 'flat'));

        fprintf('    [%d/%d] p_c = %.2f kPa, dT = %.3f K  (map amp %.4g w/h) ... ', ...
            a, numel(sel), pc/1e3, dT, Ref.Amp_wh(ip,jd));
        t0 = tic;
        Rp = qs_run_point(M, C, pt, false);
        fprintf('%.1f min\n', toc(t0)/60);

        if ~Rp.ok
            chk(false, 'point p_c=%.2f dT=%.3f failed to run: %s', pc/1e3, dT, Rp.errmsg);
            continue
        end

        amp_ok = (Ref.Amp_wh(ip,jd) >= FREQ_NEEDS_AMP) && (Rp.Amp_wh >= FREQ_NEEDS_AMP);
        pointBad = false;

        pointBad = cmp('Amp_wh',   Rp.Amp_wh,   Ref.Amp_wh(ip,jd),   true)    || pointBad;
        pointBad = cmp('MeanPeak', Rp.MeanPeak, Ref.MeanPeak(ip,jd), true)    || pointBad;
        pointBad = cmp('StdPeak',  Rp.StdPeak,  Ref.StdPeak(ip,jd),  true)    || pointBad;
        pointBad = cmp('MeanRMSE', Rp.MeanRMSE, Ref.MeanRMSE(ip,jd), true)    || pointBad;
        pointBad = cmp('StdRMSE',  Rp.StdRMSE,  Ref.StdRMSE(ip,jd),  true)    || pointBad;
        pointBad = cmp('Freq',     Rp.Freq,     Ref.Freq(ip,jd),     amp_ok)  || pointBad;

        if pointBad
            badPoints{end+1} = sprintf('t03b_attribute(%.2f, %.3f)', pc/1e3, dT); %#ok<AGROW>
        end
    end

    fprintf('    ---------------------------------------------------------\n');
    fprintf('    %d matched to %g relative, %d below the floor (not compared), %d DIFFER\n', ...
        nMatch, TOLREL, nNegl, nDiff);
    if ~isempty(badPoints)
        fprintf('    To find out WHY a point differs -- settings, or trajectory\n');
        fprintf('    sensitivity near a boundary -- run:\n');
        for k = 1:numel(badPoints)
            fprintf('        %s\n', badPoints{k});
        end
        R.msg{end+1} = sprintf('diagnose with: %s', strjoin(badPoints, ' ; '));
    end

    % ------------------------------------------------------------------
    function bad = cmp(name, got, want, compare)
        bad = false;
        if ~isfinite(want)
            return                                  % the old map had no value here
        end
        if ~compare
            fprintf('          %-9s %14.8g   (not compared: motion too small for a frequency)\n', name, got);
            nNegl = nNegl + 1;
            return
        end
        fl = FLOOR.(name);
        if abs(want) < fl && abs(got) < fl
            fprintf('          %-9s %14.8g   (both below %g -- negligible, not compared)\n', name, got, fl);
            nNegl = nNegl + 1;
            return
        end
        dv = abs(got - want) / max(abs(want), eps);
        if dv < TOLREL
            fprintf('          %-9s %14.8g   (match, %.1e)\n', name, got, dv);
            nMatch = nMatch + 1;
        else
            nDiff = nDiff + 1;
            bad = true;
            chk(false, '%s at p_c=%.2f kPa dT=%.3f K: new %.10g vs old %.10g (%.2e relative)', ...
                name, pc/1e3, dT, got, want, dv);
        end
    end
end
