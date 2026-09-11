function qs_report(D)
%QS_REPORT  Print the summary of a merged study: what the sweep found and,
%   more importantly, which of its answers should not be trusted yet.

    n = numel(D.id);
    fprintf('\n--- %s / %s -------------------------------------\n', ...
        D.meta.case_id, D.meta.study_id);
    fprintf('  points merged      : %d of %d in the list', D.meta.n_merged, D.meta.n_points_in_list);
    if D.meta.n_missing > 0
        fprintf('   (%d MISSING)', D.meta.n_missing);
    end
    fprintf('\n');
    if D.meta.n_duplicates_dropped > 0
        fprintf('  duplicates dropped : %d (re-run slices)\n', D.meta.n_duplicates_dropped);
    end
    nfail = sum(~D.ok);
    if nfail > 0
        fprintf('  FAILED points      : %d   (see errmsg in points.csv)\n', nfail);
    end

    fprintf('  regime census:\n');
    for c = 1:4
        m = sum(D.label==c & D.ok);
        if m > 0
            fprintf('    %-12s %5d  (%4.1f%%)\n', qs_label_name(c), m, 100*m/max(n,1));
        end
    end
    u = sum(D.label==0 & D.ok);
    if u>0, fprintf('    %-12s %5d\n','unknown',u); end

    fprintf('  reliability flags:\n');
    fprintf('    transitions    %5d   <- regime changed mid-run\n', sum(D.trans_flag));
    fprintf('    nonstationary  %5d   <- had not settled by t_end\n', sum(D.nonstationary));
    fprintf('    near threshold %5d   <- amplitude close to the static tolerance\n', sum(D.near_threshold));
    susp = sum(D.trans_flag | D.nonstationary | D.near_threshold);
    fprintf('    -> %d point(s) (%.1f%%) carry a provisional label\n', susp, 100*susp/max(n,1));

    if any(D.trans_flag)
        fprintf('  transitions found:\n');
        idx = find(D.trans_flag);
        [~, ord] = sort(D.trans_time(idx));
        idx = idx(ord);
        for k = 1:min(10,numel(idx))
            j = idx(k);
            fprintf('    p_c=%7.2f kPa  dT=%6.3f K   %s -> %s  at %5.2f s  (t_end %.1f s)\n', ...
                D.pc_Pa(j)/1e3, D.dT_K(j), qs_label_name(D.trans_from(j)), ...
                qs_label_name(D.trans_to(j)), D.trans_time(j), D.t_end(j));
        end
        if numel(idx) > 10
            fprintf('    ... and %d more (see points.csv)\n', numel(idx)-10);
        end
    end

    % best calibration against the DIC, among oscillating points
    if any(isfinite(D.MeanRMSE))
        score = 0.5*D.MeanRMSE + 0.5*D.StdRMSE;
        score(D.label~=2 & D.label~=3) = Inf;      % only LCO / broadband points
        if all(~isfinite(score)), score = 0.5*D.MeanRMSE + 0.5*D.StdRMSE; end
        [b, j] = min(score);
        if isfinite(b)
            fprintf('  best match to DIC  : p_c=%.2f kPa, dT=%.3f K  (%s)\n', ...
                D.pc_Pa(j)/1e3, D.dT_K(j), qs_label_name(D.label(j)));
            fprintf('    meanRMSE %.3f  stdRMSE %.3f  f=%.0f Hz  amp=%.3f w/h\n', ...
                D.MeanRMSE(j), D.StdRMSE(j), D.Freq(j), D.Amp_wh(j));
            if D.trans_flag(j) || D.nonstationary(j)
                fprintf('    WARNING: this point carries a provisional label -- re-run it longer.\n');
            end
        end
    end
    fprintf('-----------------------------------------------------------\n\n');
end
