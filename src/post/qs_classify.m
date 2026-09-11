function K = qs_classify(m, W, C)
%QS_CLASSIFY  Overall regime label and flags for one completed run.
%
%   K = QS_CLASSIFY(m, W, C)
%     m  metrics  from qs_metrics
%     W  windows  from qs_windows
%
%   The primary label is taken from the FINAL window, not from statistics over
%   the whole post-transient record. That is the deliberate choice: with a
%   possible mid-run transition, the attractor the run ended on is the answer,
%   and the statistics over a record that contains a transition are a blend of
%   two regimes and describe neither.
%
%   OUTPUT (struct K)
%     .label        final regime code (see qs_label_name)
%     .label_name   its name
%     .label_start  regime code of the first window after settle_windows
%     .flags.nonstationary   window labels changed, or amplitude drifted
%     .flags.divergent       amplitude growing at the end of the record
%     .flags.near_threshold  amplitude sits within a band of amp_tol_wh
%     .flags.short_record    fewer windows than the transition test needs
%     .amp_first, .amp_last  mean amplitude over the first / last half [w/h]
%     .growth_ratio          amp_last / amp_first
%
%   The three flags are what the refinement stage reads to decide which points
%   deserve a longer run: near_threshold and nonstationary are exactly the
%   points whose 3 s label cannot be trusted.

    settle = 0;
    if isfield(C.windows,'settle_windows'), settle = C.windows.settle_windows; end

    use = (settle+1) : W.n;
    if isempty(use), use = W.n; end
    L   = W.label(use);
    amp = W.amp_wh(use);
    t   = W.tmid(use);

    K.label       = L(end);
    K.label_name  = qs_label_name(K.label);
    K.label_start = L(1);

    % ---- amplitude drift across the record -------------------------------
    half = max(1, floor(numel(amp)/2));
    K.amp_first = mean(amp(1:half));
    K.amp_last  = mean(amp(end-half+1:end));
    denom = max([K.amp_first, K.amp_last, eps]);
    drift = abs(K.amp_last - K.amp_first) / denom;
    K.growth_ratio = K.amp_last / max(K.amp_first, eps);

    % ---- flags ------------------------------------------------------------
    K.flags.nonstationary = (drift > C.classify.stationarity_tol) || (numel(unique(L)) > 1);

    % Divergent: amplitude growing exponentially through the used windows AND
    % ending well above the static tolerance. A log-linear fit is used rather
    % than a first/last ratio so a single noisy window cannot trigger it.
    K.flags.divergent = false;
    if numel(amp) >= 4 && all(amp > 0)
        p = polyfit(t(:), log(amp(:)), 1);
        span = t(end) - t(1);
        if p(1)*span > log(C.classify.growth_tol) && K.amp_last > C.classify.amp_tol_wh
            K.flags.divergent = true;
            K.label = 4;
            K.label_name = qs_label_name(4);
        end
    end

    % Near threshold: the amplitude is close enough to the static tolerance
    % that the static/LCO call could flip with a longer run. THIS is the flag
    % that catches a slow transition hiding under a "static" label.
    band = 0.5;
    if isfield(C.classify,'near_threshold_band'), band = C.classify.near_threshold_band; end
    K.flags.near_threshold = abs(K.amp_last - C.classify.amp_tol_wh) <= band*C.classify.amp_tol_wh;

    minp = 2;
    if isfield(C.windows,'min_persist_windows'), minp = C.windows.min_persist_windows; end
    K.flags.short_record = (numel(use) < (minp + 2));
end
