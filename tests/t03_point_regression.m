function R = t03_point_regression(npoints, refpath)
%T03  Reproduce points from the ORIGINAL case-1 sweep through the new code.
%
%   R = t03_point_regression(npoints, refpath)
%     npoints  how many grid points to re-run (default 3, ~3 min each)
%     refpath  path to the old detailed_stability_map_data.mat
%              (default: searched next to the framework)
%
%   THIS IS THE TEST THAT MATTERS. The six core files moved into src/core and
%   src/aero byte-for-byte, so if the configuration plumbing did not change a
%   number anywhere, every metric must come back identical to the old map --
%   not close, identical to round-off. Anything else means a value that used
%   to be hardcoded is now being read from the wrong place.
%
%   If it fails, the report names the metric and both values, which is enough
%   to find the setting that drifted.

    if nargin < 1 || isempty(npoints), npoints = 3; end
    if nargin < 2, refpath = ''; end

    R = struct('pass',true,'skip',false,'msg',{{}});
    function chk(cond, fmt, varargin)
        if ~cond
            R.pass = false;
            R.msg{end+1} = sprintf(fmt, varargin{:});
            fprintf('    FAIL  %s\n', R.msg{end});
        end
    end

    % ---------- locate the reference map ----------------------------------
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
        R.msg{end+1} = 'reference map not found; pass its path: t03_point_regression(3, ''/path/to/detailed_stability_map_data.mat'')';
        fprintf('    SKIP: %s\n', R.msg{end});
        return
    end
    fprintf('    reference: %s\n', refpath);

    % partial load -- the file is ~145 MB and W_center is not needed
    Ref = load(refpath, 'pc_vec','dT_vec','Amp_wh','Freq','MeanRMSE','StdRMSE', ...
                        'MeanPeak','StdPeak','t_end','t_transient','dt_out');

    C = qs_config('unheated_c1_NoSBLI_Periodic','regression_unheated_c1');
    chk(abs(C.duration.base_t_end       - Ref.t_end)       < 1e-12, ...
        'study t_end %.4f does not match the reference run %.4f', C.duration.base_t_end, Ref.t_end);
    chk(abs(C.duration.base_t_transient - Ref.t_transient) < 1e-12, ...
        'study t_transient %.4f does not match the reference %.4f', C.duration.base_t_transient, Ref.t_transient);
    chk(abs(C.duration.dt_out           - Ref.dt_out)      < 1e-20, ...
        'study dt_out %g does not match the reference %g', C.duration.dt_out, Ref.dt_out);

    % ---------- choose points spread across the response ------------------
    valid = find(isfinite(Ref.Amp_wh) & isfinite(Ref.Freq));
    assert(~isempty(valid), 't03: the reference map has no finite results');
    [~, ord] = sort(Ref.Amp_wh(valid));
    take = round(linspace(1, numel(ord), min(npoints, numel(ord))));
    sel  = valid(ord(take));

    M = qs_build_model(C, true);
    npc = numel(Ref.pc_vec);

    tolrel = 1e-6;
    for a = 1:numel(sel)
        [ip, jd] = ind2sub([npc, numel(Ref.dT_vec)], sel(a));
        pc = Ref.pc_vec(ip);  dT = Ref.dT_vec(jd);

        pt = qs_pointlist('make', pc, dT, struct( ...
                't_end', C.duration.base_t_end, 't_trans', C.duration.base_t_transient, ...
                'ic_tag', 'flat'));

        fprintf('    [%d/%d] p_c = %.2f kPa, dT = %.3f K ... ', a, numel(sel), pc/1e3, dT);
        t0 = tic;
        Rp = qs_run_point(M, C, pt, false);
        fprintf('%.1f min\n', toc(t0)/60);

        if ~Rp.ok
            chk(false, 'point p_c=%.2f dT=%.3f failed to run: %s', pc/1e3, dT, Rp.errmsg);
            continue
        end

        cmp('Amp_wh',   Rp.Amp_wh,   Ref.Amp_wh(ip,jd));
        cmp('Freq',     Rp.Freq,     Ref.Freq(ip,jd));
        cmp('MeanPeak', Rp.MeanPeak, Ref.MeanPeak(ip,jd));
        cmp('StdPeak',  Rp.StdPeak,  Ref.StdPeak(ip,jd));
        cmp('MeanRMSE', Rp.MeanRMSE, Ref.MeanRMSE(ip,jd));
        cmp('StdRMSE',  Rp.StdRMSE,  Ref.StdRMSE(ip,jd));
    end

    if R.pass
        fprintf('    every metric reproduces the original map to %g relative\n', tolrel);
    end

    function cmp(name, got, want)
        if ~isfinite(want)
            return                       % the old map had no value here
        end
        d = abs(got - want) / max(abs(want), eps);
        if d < tolrel
            fprintf('          %-9s %14.8g  (match, %.1e)\n', name, got, d);
        else
            chk(false, '%s at p_c=%.2f kPa dT=%.3f K: new %.10g vs old %.10g (%.2e relative)', ...
                name, pc/1e3, dT, got, want, d);
        end
    end
end
