function S = qs_codes(section)
%QS_CODES  Print every coded value in the framework, and what each one means.
%
%   qs_codes           print the whole legend
%   qs_codes('labels') print one section: labels | flags | metrics | points |
%                      config | files | units
%   S = qs_codes(...)  also return it as a struct
%
%   Numbers stored in a results table are worthless if you cannot look up what
%   they stand for six months later. This is that lookup, and it lives in the
%   code rather than only in a document so it is always the version you are
%   actually running.
%
%   The LABEL section is derived by calling qs_label_name, so it cannot drift
%   away from what the classifier actually assigns. Everything else is
%   documentation and is kept in step by hand -- if you add a new coded value,
%   add it here and to docs/CODES.md in the same commit.

    if nargin < 1 || isempty(section), section = 'all'; end
    section = lower(section);
    show = @(s) any(strcmp(section, {'all', s}));
    Sout = struct();

    % ================= LABELS ==========================================
    codes = 0:4;
    names = cell(numel(codes),1);
    for k = 1:numel(codes), names{k} = qs_label_name(codes(k)); end
    Sout.labels = struct('code', num2cell(codes(:)), 'name', names);

    if show('labels')
        hdr('REGIME LABELS   (fields: label, label_start, trans_from, trans_to)');
        p('Assigned per analysis window. A run''s label is its FINAL window,');
        p('because with a mid-run transition the whole-record statistics blend');
        p('two regimes and describe neither.');
        blank();
        rowc('code','name','test');
        line();
        rowc('0','unknown','spectrum could not be formed (record too short)');
        rowc('1','static','window amplitude < classify.amp_tol_wh');
        rowc('2','LCO','above tolerance AND pfrac >= classify.peak_power_frac');
        rowc('3','broadband','above tolerance, no dominant tone');
        rowc('4','divergent','whole record: amplitude grows by > classify.growth_tol');
        blank();
        p('These five are the complete set. There is no label for');
        p('"nonstationary" -- that is a flag, see below.');
    end

    % ================= FLAGS ===========================================
    if show('flags')
        hdr('FLAGS   (logical, one per run, independent of the label)');
        row('name','means');
        line();
        row('nonstationary','amplitude drifted > classify.stationarity_tol between the');
        row('','first and last half, OR the window labels were not all equal.');
        row('','The drift test is skipped below classify.stationarity_amp_floor:');
        row('','a still panel''s residual is ~1e-8 w/h and its first/last ratio');
        row('','is roundoff, which would flag every static cell in the map.');
        row('drift_tested','whether that drift test was applied at all. A clean');
        row('','nonstationary = 0 with drift_tested = 0 means "never asked",');
        row('','not "settled".');
        row('near_threshold','|amp_last - amp_tol_wh| within classify.near_threshold_band');
        row('','of the tolerance: the static/LCO call could flip with more time.');
        row('divergent','same test that sets label 4.');
        row('short_record','fewer windows than the transition test needs.');
        row('trans_flag','a regime change was found that persisted to the end.');
        row('ok','the integration completed; false means see errmsg.');
    end

    % ================= METRICS =========================================
    if show('metrics')
        hdr('METRICS   (one per run; nothing is averaged across the sweep)');
        p('mean and std are over TIME within a single run; the "peak" is then');
        p('over SPACE, across the panel. Amplitudes are the CENTRE NODE only.');
        blank();
        row('name','definition','units');
        line();
        row('Amp_wh','(max-min)/(2h) of the centre node over the WHOLE','w/h');
        row('','post-transient record. Half-range, so for a sinusoid','');
        row('','it is the single-sided amplitude in panel thicknesses.','');
        row('amp_last','mean window amplitude over the LAST half. This is','w/h');
        row('','the one the label describes; prefer it on a map.','');
        row('amp_first','the same over the first half','w/h');
        row('amp_ref','max(amp_first, amp_last); gates the drift test','w/h');
        row('growth_ratio','amp_last / amp_first','-');
        row('Freq','dominant frequency of the centre node','Hz');
        row('pfrac','power within +/-5% of the peak, over total power','-');
        row('MeanPeak','max |time-mean deflection| over the panel','mm');
        row('StdPeak','max time-std of deflection over the panel','mm');
        row('MeanRMSE','||model mean - DIC mean|| / ||DIC mean||','-');
        row('StdRMSE','the same for the std field','-');
        row('signflip','+1 or -1: orientation used to match the DIC sign','-');
        row('trans_time','transition time estimate (centre of the first window','s');
        row('','carrying the new label; biased LATE by up to half a window)','');
        row('trans_t_lo','that window''s start -- the earliest consistent time','s');
        row('trans_t_hi','that window''s end','s');
        row('wall_s','wall-clock time for this point','s');
        blank();
        p('Window table (kept for EVERY point): tmid, amp_wh, wmean, fdom,');
        p('pfrac, label. A few hundred bytes, which is what makes transitions');
        p('visible without storing full histories.');
    end

    % ================= POINT LIST ======================================
    if show('points')
        hdr('POINT LIST COLUMNS');
        row('name','meaning');
        line();
        row('id','pc<Pa>_dT<K>_t<s>_<ic>, e.g. pc051927_dT0012p760_t003p00_flat');
        row('level','0 = original coarse grid; 1, 2, ... = refinement pass number');
        row('origin','coarse | refine | manual | long');
        row('ic_tag','flat | flat_neg | rest | from:<point id>');
        row('','flat    = +seed on every mode (what every old sweep used)');
        row('','flat_neg= -seed; the cheapest basin probe');
        row('','rest    = exactly zero; linear checks only');
        row('save_hist','force keeping the full centre history for this point');
        row('t_end','run length for THIS point, in seconds');
        row('t_trans','transient discarded before statistics, in seconds');
    end

    % ================= CONFIG ENUMS ====================================
    if show('config')
        hdr('CONFIGURATION VALUES THAT ARE NAMES, NOT NUMBERS');
        row('setting','allowed values');
        line();
        row('study.mode','sweep | long');
        row('solver.integrator','ode15s | ode45 | ode23s | ode113');
        row('ic.mode','flat | flat_neg | rest');
        row('aero.c3_mode','vandyke2 (c3 = 0, the published Eq.8)');
        row('','classical_cubic ((gamma+1)/12)');
        row('store.save_history_if','any of: transition | nonstationary |');
        row('','near_threshold | divergent | all | none');
        row('store.history_precision','single (default) | double');
        row('','single puts a ~1e-7 w/h floor on run-to-run comparison');
        row('duration_policy.promote_if','any of: boundary | nonstationary |');
        row('','near_threshold | transition | divergent');
        row('fields.mode','per_second (default) | per_step | single');
        row('fields.precision','single | double');
        row('case.net_load_sign_convention','+1 or -1; absent = the sign gate stays off');
        blank();
        p('Preflight returns F.checks with one entry per gate:');
        p('   T0 | load_sign | buckling | box | mach     each with a .pass');
    end

    % ================= FILES ===========================================
    if show('files')
        hdr('FILE AND FOLDER NAMES');
        row('pattern','what it is');
        line();
        row('pointlist.mat / .csv','the study''s list of runs');
        row('manifest.json','resolved settings + git commit that produced them');
        row('points/partial_<k>_of_<n>.mat','one worker''s results, checkpointed');
        row('merged/points.mat / .csv','all workers joined by run id');
        row('histories/<id>_hist.mat','full centre history, for flagged points');
        row('histories/<id>_state.mat','final [q; qdot], for continuing a run');
        row('figures/fig_regime_map.*','the three-panel map');
        row('<id>/postproc_data.mat','long runs: full modal history');
        row('<id>/fields/fields_sec_NNN.mat','long runs: one file per simulated second');
    end

    % ================= UNITS ===========================================
    if show('units')
        hdr('UNITS -- WHERE kPa BECOMES Pa');
        row('quantity','JSON files','inside the code');
        line();
        row('cavity pressure','kPa (grid.pc_kPa)','Pa (pc_Pa, pc_nom_Pa)');
        row('temperature rise','K','K');
        row('amplitude','-','w/h (panel thicknesses)');
        row('deflection peaks','-','mm');
        row('pressures in results','-','Pa');
        row('time','s','s');
        blank();
        p('Rule of thumb: anything ending _Pa or _K or _m is SI; a JSON grid');
        p('block is in kPa because that is what you read off a gauge.');
    end

    if show('all')
        blank();
        p('Full write-up: docs/CODES.md, and the design memo in docs/.');
        line();
    end

    if nargout > 0, S = Sout; end
end

% ---------------------------------------------------------------- printing
function hdr(t)
    fprintf('\n%s\n%s\n', t, repmat('=', 1, min(72, numel(t))));
end
function p(t),     fprintf('  %s\n', t); end
function blank(),  fprintf('\n'); end
function line(),   fprintf('  %s\n', repmat('-', 1, 72)); end
function rowc(a, b, c)
    fprintf('  %-6s %-13s %s\n', a, b, c);
end
function row(a, b, c)
    if nargin < 3
        fprintf('  %-26s %s\n', a, b);
    else
        fprintf('  %-14s %-48s %s\n', a, b, c);
    end
end
