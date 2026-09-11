function G = qs_plot_map(caseId, studyId)
%QS_PLOT_MAP  Regime map for a merged study.
%
%   G = QS_PLOT_MAP(caseId, studyId)
%
%   Three panels:
%     1. regime map    -- one colour per label, with transition points marked
%                         by a ring: fill = regime it ended in, ring = regime
%                         it started in. A ringed point is one whose 3 s answer
%                         would have been wrong.
%     2. amplitude     -- peak-to-peak w/h, log colour scale
%     3. frequency     -- dominant centre frequency, oscillating points only
%
%   The case nominal and its +/-1 sigma measurement box are drawn on every
%   panel, because the only question the map finally has to answer is which
%   regime the experiment's own operating point falls in, and how far the
%   boundary is from it in units of the measurement error.
%
%   Cells with no point are left blank -- after refinement the lattice is
%   deliberately ragged, and filling the holes by interpolation would invent
%   boundary detail that was never computed.

    C = qs_config(caseId, studyId);
    S = load(fullfile(C.paths.merged_dir,'points.mat'));
    D = S.D;
    G = qs_rasterize(D);

    pcv = G.pc_kPa;  dTv = G.dT_K;
    fh = figure('Position',[60 60 1500 470],'Color','w');

    % ---------------- 1. regime map ---------------------------------------
    subplot(1,3,1);
    L = G.label;  L(~G.filled) = NaN;
    imagesc(dTv, pcv, L, 'AlphaData', ~isnan(L));
    set(gca,'YDir','normal');  colormap(gca, regime_colors());
    caxis([0.5 4.5]); %#ok<CAXIS>
    cb = colorbar('Ticks',1:4,'TickLabels',{'static','LCO','broadband','divergent'});
    set(get(cb,'Label'),'String','regime at t_{end}');
    xlabel('\DeltaT [K]');  ylabel('p_c [kPa]');
    title(sprintf('Regime map  (%d points, %d refined)', sum(G.filled(:)), sum(G.level(G.filled)>0)));
    hold on;
    mark_transitions(D);
    mark_nominal(C);

    % ---------------- 2. amplitude -----------------------------------------
    subplot(1,3,2);
    A = G.Amp_wh;  A(~G.filled) = NaN;
    imagesc(dTv, pcv, log10(max(A,1e-4)), 'AlphaData', ~isnan(A));
    set(gca,'YDir','normal');  colormap(gca, parula_safe());
    cb = colorbar;  set(get(cb,'Label'),'String','log_{10} (w_{pp}/2h)');
    xlabel('\DeltaT [K]');  ylabel('p_c [kPa]');
    title('Peak-to-peak amplitude');
    hold on;  mark_nominal(C);

    % ---------------- 3. frequency ------------------------------------------
    subplot(1,3,3);
    F = G.Freq;  F(~G.filled | G.label==1) = NaN;      % blank the static region
    imagesc(dTv, pcv, F, 'AlphaData', ~isnan(F));
    set(gca,'YDir','normal');  colormap(gca, parula_safe());
    cb = colorbar;  set(get(cb,'Label'),'String','f [Hz]');
    xlabel('\DeltaT [K]');  ylabel('p_c [kPa]');
    title('Dominant frequency (oscillating points)');
    hold on;  mark_nominal(C);

    try
        sgtitle(sprintf('%s / %s   [%s]', caseId, studyId, C.meta.git_commit), ...
                'Interpreter','none');
    catch
    end

    if exist(C.paths.fig_dir,'dir')~=7, mkdir(C.paths.fig_dir); end
    qs_savefig(fh, fullfile(C.paths.fig_dir,'fig_regime_map'), struct('G',G));
    fprintf('qs_plot_map: wrote %s\n', fullfile(C.paths.fig_dir,'fig_regime_map.png'));
end

% ======================================================================
function mark_transitions(D)
%MARK_TRANSITIONS  Ring every point whose regime changed mid-run.
    idx = find(D.trans_flag & D.ok);
    if isempty(idx), return; end
    cols = regime_colors();
    for k = idx(:).'
        cf = clampcol(cols, D.trans_from(k));
        ct = clampcol(cols, D.trans_to(k));
        plot(D.dT_K(k), D.pc_Pa(k)/1e3, 'o', 'MarkerSize', 9, ...
             'MarkerFaceColor', ct, 'MarkerEdgeColor', cf, 'LineWidth', 2);
    end
end

function mark_nominal(C)
%MARK_NOMINAL  The measured operating point and its +/-1 sigma box.
    pc = C.case.pc_nom_Pa/1e3;  dT = C.case.dT_nom_K;
    plot(dT, pc, 'kp', 'MarkerSize', 14, 'MarkerFaceColor','w', 'LineWidth',1.5);
    spc = 0;  sdT = 0;
    if isfield(C.case,'pc_sigma_Pa'), spc = C.case.pc_sigma_Pa/1e3; end
    if isfield(C.case,'dT_sigma_K'),  sdT = C.case.dT_sigma_K; end
    if spc > 0 && sdT > 0
        plot([dT-sdT dT+sdT dT+sdT dT-sdT dT-sdT], ...
             [pc-spc pc-spc pc+spc pc+spc pc-spc], 'k--', 'LineWidth',1);
    end
end

function c = regime_colors()
%REGIME_COLORS  1 static, 2 LCO, 3 broadband, 4 divergent.
    c = [0.85 0.87 0.90;    % static      - pale grey
         0.20 0.45 0.80;    % LCO         - blue
         0.90 0.55 0.15;    % broadband   - orange
         0.75 0.15 0.20];   % divergent   - red
end

function c = clampcol(cols, code)
    if code >= 1 && code <= size(cols,1), c = cols(code,:); else, c = [0 0 0]; end
end

function m = parula_safe()
    try, m = parula(256); catch, m = jet(256); end
end
