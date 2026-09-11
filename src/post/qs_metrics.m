function m = qs_metrics(M, C, tvec, q, t_trans)
%QS_METRICS  Scalar response metrics for one completed run.
%
%   m = QS_METRICS(M, C, tvec, q, t_trans)
%     M       model from qs_build_model
%     tvec    output time vector [s]
%     q       modal displacement history (Nt x nm)
%     t_trans transient to discard before computing statistics [s]
%
%   Every formula here is byte-for-byte what run_detailed_stability_map.m
%   computed inline. It is factored out for one reason: the sweep and the
%   long 20 s run must produce the SAME numbers from the same history, so
%   a long run can be dropped onto a stability map without recalibrating
%   anything. Changing a definition here changes it in both places at once.
%
%   OUTPUT (struct m)
%     .Amp_wh    peak-to-peak centre deflection / (2h)     [-]
%     .MeanPeak  max |mean deflection| over the panel      [mm]
%     .StdPeak   max std of deflection over the panel      [mm]
%     .MeanRMSE  normalised RMSE of mean field vs DIC      [-]  (NaN if no DIC)
%     .StdRMSE   normalised RMSE of std field  vs DIC      [-]  (NaN if no DIC)
%     .Freq      dominant centre-deflection frequency      [Hz]
%     .pfrac     peak power fraction (tone vs broadband)   [-]
%     .wc        centre deflection history, full record    [m] (Nt x 1)
%     .signflip  +1/-1 sign applied to match the DIC mean field

    Fs   = 1/C.duration.dt_out;
    mask = tvec(:) >= t_trans;
    assert(any(mask), 'qs_metrics: transient cut leaves no samples (t_trans=%g)', t_trans);

    % ---- centre deflection history --------------------------------------
    wc = q * M.rom.Modes(:, M.cidx);
    m.wc = wc(:);

    % ---- modal statistics over the post-transient window -----------------
    Qs   = q(mask,:);
    qbar = mean(Qs,1).';
    Cq   = cov(Qs);

    meanf = M.Phi * qbar;
    stdf  = sqrt(max(sum((M.Phi*Cq).*M.Phi, 2), 0));
    m.MeanPeak = max(abs(meanf)) * 1e3;      % mm
    m.StdPeak  = max(stdf)       * 1e3;      % mm

    % ---- amplitude in panel thicknesses ----------------------------------
    segc = wc(mask);
    m.Amp_wh = (max(segc) - min(segc)) / (2*M.h);

    % ---- spectrum of the centre history ----------------------------------
    bf = 0.05;
    if isfield(C.classify,'peak_band_frac'), bf = C.classify.peak_band_frac; end
    S = qs_spectrum(segc, Fs, bf);
    m.Freq  = S.fpk;
    m.pfrac = S.pfrac;

    % ---- comparison against the DIC measurement --------------------------
    m.MeanRMSE = NaN;  m.StdRMSE = NaN;  m.signflip = 1;
    if M.dic.present
        mdic = reshape(M.dic.Phi_dic*qbar, size(M.dic.X));
        sdic = reshape(sqrt(max(sum((M.dic.Phi_dic*Cq).*M.dic.Phi_dic,2),0)), size(M.dic.X));
        % The ROM sign convention is arbitrary relative to the DIC's; take
        % whichever orientation is closer, and record which was used.
        if norm(-mdic(:)-M.dic.mean(:)) < norm(mdic(:)-M.dic.mean(:))
            mdic = -mdic;  m.signflip = -1;
        end
        m.MeanRMSE = norm(mdic(:)-M.dic.mean(:)) / M.dic.nrm_mean;
        m.StdRMSE  = norm(sdic(:)-M.dic.std(:))  / M.dic.nrm_std;
    end
end
