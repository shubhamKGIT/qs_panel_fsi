function T = qs_detect_transition(W, C)
%QS_DETECT_TRANSITION  Find a mid-run regime change in a window label sequence.
%
%   T = QS_DETECT_TRANSITION(W, C)   with W from qs_windows.
%
%   Rule: walk back from the end of the record and find the last window whose
%   label differs from the final label. The transition time is the centre of
%   the window immediately after it. The final label must then persist for at
%   least C.windows.min_persist_windows windows, otherwise the change is too
%   close to the end of the run to distinguish from a fluctuation and no
%   transition is reported (but flags.nonstationary will still be set by
%   qs_classify, which is what promotes the point to a longer run).
%
%   Windows before C.windows.settle_windows are ignored, so the initial
%   start-up from the flat initial condition is not reported as a transition.
%
%   TIMING RESOLUTION. The reported time is the centre of the first window
%   that carries the new label, so it is biased LATE by up to about half a
%   window plus however long the physical change took to move the amplitude
%   across the classification threshold. time_lo / time_hi bracket it with
%   that window's own start and end. Do not quote a transition time to better
%   than the window length; shorten C.windows.T_win_s if you need to.
%
%   OUTPUT (struct T)
%     .flag      true if a persistent regime change was found
%     .time      transition time estimate [s]  (NaN when flag is false)
%     .time_lo   start of the first window with the new label [s]
%     .time_hi   end   of that window [s]
%     .from,.to  regime codes before and after
%     .from_name,.to_name
%     .n_changes number of label changes in the used part of the sequence
%     .persist   how many windows the final label lasted

    settle = 0;
    if isfield(C.windows,'settle_windows'), settle = C.windows.settle_windows; end
    minp = 2;
    if isfield(C.windows,'min_persist_windows'), minp = C.windows.min_persist_windows; end

    use = (settle+1) : W.n;
    if isempty(use), use = W.n; end
    L  = W.label(use);
    t  = W.tmid(use);
    t0 = W.t0(use);
    t1 = W.t1(use);

    T = struct('flag',false, 'time',NaN, 'time_lo',NaN, 'time_hi',NaN, ...
               'from',L(1), 'to',L(end), ...
               'from_name',qs_label_name(L(1)), 'to_name',qs_label_name(L(end)), ...
               'n_changes',sum(diff(L(:))~=0), 'persist',0);

    k = find(L ~= L(end), 1, 'last');
    if isempty(k)
        T.persist = numel(L);
        return                         % one label throughout: no transition
    end

    idx = k + 1;                       % first window of the final label
    T.persist = numel(L) - idx + 1;
    if T.persist < minp
        return                         % too close to the end to trust
    end

    T.flag      = true;
    T.time      = t(idx);
    T.time_lo   = t0(idx);
    T.time_hi   = t1(idx);
    T.from      = L(k);
    T.to        = L(end);
    T.from_name = qs_label_name(T.from);
    T.to_name   = qs_label_name(T.to);
end
