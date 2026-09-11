function qs_plot_history(C, R, outdir)
%QS_PLOT_HISTORY  Centre-deflection history with the window signature under it.
%
%   qs_plot_history(C, R, outdir)
%     R must carry .wc, .tvec and .W (a long run always does; a swept point
%     does when the storage policy kept its history).
%
%   Three stacked panels sharing a time axis:
%     1. centre deflection, with the transient cut and any transition marked
%     2. window amplitude (w/h) and the static tolerance line
%     3. window regime label and dominant frequency
%
%   Reading them together is the point: a transition shows up as the moment
%   the amplitude track crosses the tolerance and the label strip changes
%   colour, which is far easier to judge than staring at the raw trace.

    if nargin < 3 || isempty(outdir), outdir = C.paths.fig_dir; end
    if exist(outdir,'dir')~=7, mkdir(outdir); end
    assert(~isempty(R.wc) && ~isempty(R.W), 'qs_plot_history: record has no stored history');

    t  = double(R.tvec(:));
    wc = double(R.wc(:))*1e3;          % mm
    W  = R.W;

    fh = figure('Position',[60 60 1100 780],'Color','w');

    % ---- 1. trace ---------------------------------------------------------
    ax1 = subplot(3,1,1);
    plot(t, wc, 'LineWidth', 0.5);  grid on;  hold on;
    yl = ylim;
    plot([R.t_trans R.t_trans], yl, 'r--', 'LineWidth', 1);
    if R.trans_flag
        plot([R.trans_time R.trans_time], yl, 'k-', 'LineWidth', 2);
        text(R.trans_time, yl(2), sprintf('  %s \\rightarrow %s', ...
            qs_label_name(R.trans_from), qs_label_name(R.trans_to)), ...
            'VerticalAlignment','top','FontWeight','bold');
    end
    ylabel('centre w [mm]');
    title(sprintf('%s   p_c = %.2f kPa,  \\DeltaT = %.3f K,  %.0f s   ->  %s', ...
        R.id, R.pc_Pa/1e3, R.dT_K, R.t_end, R.label_name), 'Interpreter','tex');

    % ---- 2. window amplitude ----------------------------------------------
    ax2 = subplot(3,1,2);
    semilogy(double(W.tmid), max(double(W.amp_wh),1e-4), '-o', 'MarkerSize',3); grid on; hold on;
    yl = ylim;
    plot([t(1) t(end)], [C.classify.amp_tol_wh C.classify.amp_tol_wh], 'r--','LineWidth',1);
    text(t(1), C.classify.amp_tol_wh, ' static tolerance', 'VerticalAlignment','bottom','Color','r');
    if R.trans_flag, plot([R.trans_time R.trans_time], yl, 'k-','LineWidth',2); end
    ylabel('window w_{pp}/2h');

    % ---- 3. label strip + frequency ---------------------------------------
    ax3 = subplot(3,1,3);
    lab = double(W.label);
    cols = [0.85 0.87 0.90; 0.20 0.45 0.80; 0.90 0.55 0.15; 0.75 0.15 0.20];
    hold on;
    for k = 1:numel(lab)
        c = [0 0 0];
        if lab(k) >= 1 && lab(k) <= 4, c = cols(lab(k),:); end
        plot(double(W.tmid(k)), 1, 's', 'MarkerSize',10, 'MarkerFaceColor',c, 'MarkerEdgeColor','none');
    end
    yyaxis right
    plot(double(W.tmid), double(W.fdom), '.-');
    ylabel('window f_{dom} [Hz]');
    yyaxis left
    set(gca,'YTick',[]);  ylim([0.5 1.5]);  grid on;
    xlabel('t [s]');
    ylabel('regime');

    try, linkaxes([ax1 ax2 ax3],'x'); catch, end
    xlim(ax1,[t(1) t(end)]);

    base = fullfile(outdir, sprintf('fig_history_%s', R.id));
    qs_savefig(fh, base, struct('tvec',R.tvec,'wc',R.wc,'W',W,'record',rmfield_safe(R,{'wc','tvec'})));
    fprintf('  figure -> %s.png\n', base);
end

function S = rmfield_safe(S, flds)
    for k = 1:numel(flds)
        if isfield(S,flds{k}), S = rmfield(S,flds{k}); end
    end
end
