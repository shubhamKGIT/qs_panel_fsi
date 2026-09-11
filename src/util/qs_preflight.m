function F = qs_preflight(M, C)
%QS_PREFLIGHT  Sanity checks that must pass before spending compute on a sweep.
%
%   F = QS_PREFLIGHT(M, C)
%
%   Three checks, each of which has already cost this project a full sweep:
%
%   1. AERO TEMPERATURE. The implied stagnation temperature at the BL edge,
%         T0 = (a_edge^2 / (gamma*R)) * (1 + 0.2*M_edge^2)
%      must match the case's operating T0 and be near-constant across the
%      panel. A cold-flow RANS behind a heated case makes a_edge ~15% too
%      low, which makes the piston-theory damping term Zdot/a_l ~15% too
%      large, which suppresses limit cycles and shrinks the flutter region.
%      The map still looks plausible; it is just wrong.
%
%   2. STATIC LOAD SIGN. mean(p_l) - p_c is the mean net load. Its sign must
%      agree with the direction of the measured DIC mean deflection. If it
%      does not, the panel is being pushed the wrong way and every mean-field
%      comparison in the sweep is meaningless.
%
%   3. BUCKLING THRESHOLD. dT_cr from the smallest eigenvalue of
%      Ke*v = dT*Kt*v, reported against the study's dT range. Sweeping far
%      below dT_cr spends time on a flat unbuckled panel; sweeping without
%      knowing where it sits makes the map hard to read.
%
%   Set preflight.strict = true in the study to make a failure an error
%   instead of a warning.

    R_gas = 287.058;
    g     = C.aero.gamma;

    fprintf('\n--- preflight: %s / %s ---------------------------\n', ...
        C.case.case_id, C.study.study_id);

    F = struct('pass',true, 'checks',struct());

    % ---------- 1. implied stagnation temperature --------------------------
    T_edge = (M.aero.al.^2) / (g*R_gas);
    T0     = T_edge .* (1 + 0.5*(g-1)*M.aero.Ml.^2);
    T0min  = min(T0);  T0max = max(T0);  T0mean = mean(T0);

    T0exp = NaN;
    if isfield(C.case,'T0_K'), T0exp = C.case.T0_K; end
    tol = 0.05;
    if isfield(C.preflight,'T0_tol_frac'), tol = C.preflight.T0_tol_frac; end

    fprintf('  aero T0 implied   : %.1f K  (spread %.1f .. %.1f K)\n', T0mean, T0min, T0max);
    okT = true;
    if isfinite(T0exp)
        rel = abs(T0mean - T0exp)/T0exp;
        fprintf('  case T0 expected  : %.1f K   -> %.1f%% difference\n', T0exp, 100*rel);
        if rel > tol
            okT = false;
            note(C, 'AERO TEMPERATURE MISMATCH: the aero file was computed at T0 = %.1f K but this case ran at %.1f K.\n    Damping scales as 1/a_edge, so it is off by a factor %.3f. Limit cycles will be %s.', ...
                T0mean, T0exp, sqrt(T0exp/T0mean), ternary(T0mean<T0exp,'suppressed','exaggerated'));
        end
    else
        fprintf('  case T0 expected  : (not set in case.json -- add "T0_K" to enable this check)\n');
    end
    if (T0max - T0min)/T0mean > 0.02
        okT = false;
        note(C, 'implied T0 varies by %.1f%% across the panel; an adiabatic BL-edge solution should be near-constant. Check the edge profile extraction.', ...
            100*(T0max-T0min)/T0mean);
    end
    F.checks.T0 = struct('pass',okT,'mean',T0mean,'min',T0min,'max',T0max,'expected',T0exp);
    F.pass = F.pass && okT;

    % ---------- 2. mean net load sign --------------------------------------
    pl_mean = mean(M.aero.pl);
    pc      = C.case.pc_nom_Pa;
    dp      = pl_mean - pc;
    fprintf('  mean p_l          : %.2f kPa   p_c nominal: %.2f kPa\n', pl_mean/1e3, pc/1e3);
    fprintf('  mean net load     : %+.3f kPa\n', dp/1e3);

    okS = true;
    if M.dic.present
        dic_mean_signed = mean(M.dic.mean(:));
        fprintf('  DIC mean deflect. : %+.4f mm (panel-average)\n', dic_mean_signed*1e3);
        if ~isfield(C.case,'net_load_sign_convention')
            fprintf('  (set case.net_load_sign_convention = +1 or -1 to enable the sign check)\n');
        elseif abs(dic_mean_signed) > 1e-6 && abs(dp) > 50
            % A positive net load (p_l > p_c) pushes the panel away from the
            % flow. Whether that is +w or -w depends on the ROM's sign
            % convention, which is recorded per case.
            expect = C.case.net_load_sign_convention;
            if sign(dp)*expect ~= sign(dic_mean_signed)
                okS = false;
                note(C, ['STATIC LOAD SIGN: mean net load is %+.2f kPa but the measured mean deflection is %+.3f mm.\n' ...
                    '    They disagree. Either p_l, p_c or net_load_sign_convention is wrong --\n' ...
                    '    fix this before sweeping, every mean-field RMSE depends on it.'], dp/1e3, dic_mean_signed*1e3);
            end
        end
    else
        fprintf('  DIC mean deflect. : (no reference file -- sign check skipped)\n');
    end
    F.checks.load_sign = struct('pass',okS,'pl_mean',pl_mean,'pc',pc,'dp',dp);
    F.pass = F.pass && okS;

    % ---------- 3. buckling threshold and first frequency -------------------
    dT_cr = NaN;
    try
        ev = eig(M.rom.Ke, M.rom.Kt);
        ev = ev(isfinite(ev) & real(ev) > 0);
        if ~isempty(ev), dT_cr = min(real(ev)); end
    catch
    end
    f1 = NaN;
    try
        w2 = eig(M.rom.Ke, M.rom.M);
        w2 = sort(real(w2(isfinite(w2) & real(w2)>0)));
        if ~isempty(w2), f1 = sqrt(w2(1))/(2*pi); end
    catch
    end
    fprintf('  dT_cr (buckling)  : %.2f K\n', dT_cr);
    fprintf('  f1 (unloaded ROM) : %.1f Hz', f1);
    if isfield(C.case,'f1_installed_Hz') && isfinite(C.case.f1_installed_Hz)
        fprintf('   (measured installed: %.1f Hz)', C.case.f1_installed_Hz);
    end
    fprintf('\n');

    okB = true;
    if isfield(C.study,'grid') && isfinite(dT_cr)
        dTv = qs_axis(C.study.grid.dT_K,'grid.dT_K');
        below = sum(dTv < dT_cr);
        fprintf('  study dT range    : %.2f .. %.2f K  (%d of %d values below dT_cr)\n', ...
            min(dTv), max(dTv), below, numel(dTv));
        if below > 0.5*numel(dTv)
            okB = false;
            note(C, 'more than half the dT grid sits below the buckling threshold %.2f K; most of the sweep will be a flat unbuckled panel.', dT_cr);
        end
    end
    F.checks.buckling = struct('pass',okB,'dT_cr',dT_cr,'f1',f1);
    F.pass = F.pass && okB;

    % ---------- 4. does the study's box contain the operating point? ---------
    %  The sweep exists to answer "which regime is the EXPERIMENT in, and how
    %  far away is the boundary". If the measured operating point is not inside
    %  the swept box, the sweep cannot answer that question no matter how well
    %  it runs. This is the cheapest possible check and it catches the most
    %  expensive mistake: pointing a study written for one case at another case
    %  whose cavity pressure sits somewhere else entirely.
    okX = true;
    F.checks.box = struct('pass',true);
    if isfield(C.study,'grid')
        pcv = qs_axis(C.study.grid.pc_kPa,'grid.pc_kPa') * 1e3;
        dTv = qs_axis(C.study.grid.dT_K,  'grid.dT_K');
        pcn = C.case.pc_nom_Pa;   dTn = C.case.dT_nom_K;
        spc = getfielddef(C.case,'pc_sigma_Pa', NaN);
        sdT = getfielddef(C.case,'dT_sigma_K',  NaN);

        fprintf('  study p_c range   : %.2f .. %.2f kPa   (nominal %.2f kPa)\n', ...
            min(pcv)/1e3, max(pcv)/1e3, pcn/1e3);
        fprintf('  study dT  range   : %.2f .. %.2f K     (nominal %.3f K)\n', ...
            min(dTv), max(dTv), dTn);

        [okpc, mpc] = inside('p_c', pcn, pcv, spc, 'kPa', 1e3);
        [okdT, mdT] = inside('dT',  dTn, dTv, sdT, 'K',   1);
        okX = okpc && okdT;
        if ~okX
            note(C, ['STUDY BOX DOES NOT CONTAIN THE OPERATING POINT.\n' ...
                     '    %s%s\n' ...
                     '    This sweep cannot say which regime the experiment is in.\n' ...
                     '    Copy the study file, move the range onto the nominal, and re-init.'], mpc, mdT);
        else
            % inside the box, but is there room to see the boundary either side?
            margin_pc = min(pcn-min(pcv), max(pcv)-pcn);
            margin_dT = min(dTn-min(dTv), max(dTv)-dTn);
            if isfinite(spc) && spc > 0 && margin_pc < 3*spc
                note(C, 'the nominal p_c sits only %.1f sigma from the edge of the swept range; widen it or the boundary may fall outside.', margin_pc/spc);
                okX = false;
            end
            if isfinite(sdT) && sdT > 0 && margin_dT < 1*sdT
                note(C, 'the nominal dT sits only %.1f sigma from the edge of the swept range.', margin_dT/sdT);
                okX = false;
            end
        end
        F.checks.box = struct('pass',okX, 'pc_range',[min(pcv) max(pcv)], ...
                              'dT_range',[min(dTv) max(dTv)], 'pc_nom',pcn, 'dT_nom',dTn);
    end
    F.pass = F.pass && okX;

    % ---------- 5. is piston theory usable everywhere on the panel? ---------
    %  c1 = M/sqrt(M^2-1) is singular at M = 1. A BL-edge Mach number at or
    %  below 1 anywhere on the panel makes the coefficients complex or huge,
    %  and nothing downstream would tell you -- the run just produces nonsense.
    Mlo = min(M.aero.Ml);
    fprintf('  min BL-edge Mach  : %.3f\n', Mlo);
    okM = true;
    if ~isreal(M.aero.c1) || ~isreal(M.aero.c2)
        okM = false;
        note(C, 'the Van Dyke coefficients are COMPLEX: some BL-edge Mach numbers are subsonic. Piston theory does not apply here.');
    elseif Mlo < 1.05
        okM = false;
        note(C, 'minimum BL-edge Mach is %.3f. c1 = M/sqrt(M^2-1) is %.1f there -- piston theory is stretched thin near M = 1.', Mlo, Mlo/sqrt(Mlo^2-1));
    end
    F.checks.mach = struct('pass',okM,'Ml_min',Mlo);
    F.pass = F.pass && okM;

    % ---------- verdict ------------------------------------------------------
    if F.pass
        fprintf('  PREFLIGHT PASSED\n');
    else
        fprintf('  PREFLIGHT RAISED ISSUES (see above)\n');
        if isfield(C.preflight,'strict') && C.preflight.strict
            error('qs_preflight: strict mode -- refusing to run. Fix the issues or set preflight.strict = false.');
        end
    end
    fprintf('-----------------------------------------------------------\n\n');
end

% ======================================================================
function note(C, fmt, varargin)
    msg = sprintf(fmt, varargin{:});
    fprintf('  !! %s\n', msg);
    if ~(isfield(C.preflight,'strict') && C.preflight.strict)
        warning('qs_preflight:issue', '%s', msg);
    end
end

function s = ternary(c,a,b)
    if c, s = a; else, s = b; end
end

function v = getfielddef(S, f, d)
    if isstruct(S) && isfield(S,f) && ~isempty(S.(f)), v = S.(f); else, v = d; end
end

function [ok, msg] = inside(name, nom, vec, sigma, unit, scale)
%INSIDE  Is the nominal value within the swept range, and by how much?
    ok = (nom >= min(vec)) && (nom <= max(vec));
    msg = '';
    if ~ok
        if nom < min(vec), gap = min(vec) - nom; where = 'BELOW';
        else,              gap = nom - max(vec); where = 'ABOVE'; end
        extra = '';
        if isfinite(sigma) && sigma > 0
            extra = sprintf(' = %.0f sigma', gap/sigma);
        end
        msg = sprintf('%s nominal %.3f %s is %s the swept range by %.3f %s%s. ', ...
            name, nom/scale, unit, where, gap/scale, unit, extra);
    end
end
