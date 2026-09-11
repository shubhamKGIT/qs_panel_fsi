function R = t02_model_and_preflight()
%T02  Model assembly and the preflight gates.
%   Loads the real unheated case-1 data and checks that qs_build_model
%   reproduces what the original driver built inline, then checks that
%   preflight actually detects the known cold-flow aero problem.

    R = struct('pass',true,'skip',false,'msg',{{}});
    function chk(cond, fmt, varargin)
        if ~cond
            R.pass = false;
            R.msg{end+1} = sprintf(fmt, varargin{:});
            fprintf('    FAIL  %s\n', R.msg{end});
        end
    end

    C = qs_config('unheated_c1_NoSBLI_Periodic','regression_unheated_c1');
    M = qs_build_model(C);

    % ---------- shapes and the values the driver hardcoded ----------------
    chk(M.rom.nx==101 && M.rom.ny==51, 'ROM grid is %dx%d, expected 101x51', M.rom.nx, M.rom.ny);
    chk(M.nm == size(M.rom.Modes,1),   'mode count mismatch');
    chk(size(M.Phi,1)==M.rom.nx*M.rom.ny && size(M.Phi,2)==M.nm, 'Phi has the wrong shape');
    fprintf('    %d modes on a %dx%d grid\n', M.nm, M.rom.nx, M.rom.ny);

    cidx_old = round(M.rom.nx/2) + (round(M.rom.ny/2)-1)*M.rom.nx;
    chk(M.cidx == cidx_old, 'centre node index changed: %d vs %d', M.cidx, cidx_old);

    chk(abs(M.h - 6.35e-4) < 1e-12, 'panel thickness is %g, expected 6.35e-4', M.h);
    chk(M.aero.zsign == -1, 'zsign is %d, expected -1 (dissipative single-surface piston theory)', M.aero.zsign);
    chk(abs(M.aero.pinf - 50287) < 1e-9, 'aero pinf is %g, expected 50287', M.aero.pinf);
    chk(numel(M.aero.pl) == M.rom.nx*M.rom.ny, 'aero pl has %d nodes, expected %d', ...
        numel(M.aero.pl), M.rom.nx*M.rom.ny);
    chk(all(isfinite(M.aero.c1)) && all(isfinite(M.aero.c2)), 'Van Dyke coefficients are not finite');
    chk(all(M.aero.Ml > 1), 'some BL-edge Mach numbers are subsonic; piston theory c1 would be complex');

    % ---------- DIC basis --------------------------------------------------
    chk(M.dic.present, 'DIC reference was not loaded from %s', C.paths.ref_file);
    if M.dic.present
        chk(size(M.dic.Phi_dic,1)==numel(M.dic.X), 'Phi_dic rows (%d) != DIC points (%d)', ...
            size(M.dic.Phi_dic,1), numel(M.dic.X));
        chk(size(M.dic.Phi_dic,2)==M.nm, 'Phi_dic columns (%d) != modes (%d)', size(M.dic.Phi_dic,2), M.nm);
        chk(all(isfinite(M.dic.Phi_dic(:))), 'Phi_dic contains non-finite values (bad interpolation)');
        chk(M.dic.nrm_mean > 0 && M.dic.nrm_std > 0, 'DIC reference norms are zero');

        % rebuild the basis the way the ORIGINAL driver did and compare
        xv = M.rom.x_grid(:,1);  yv = M.rom.y_grid(1,:).';
        xq = M.dic.X(:);  yq = M.dic.Y(:) + 0.5*0.127;
        Pd = zeros(numel(xq), M.nm);
        for j = 1:M.nm
            Fj = griddedInterpolant({xv,yv}, reshape(M.rom.Modes(j,:),M.rom.nx,M.rom.ny), 'linear','nearest');
            Pd(:,j) = Fj(xq, yq);
        end
        d = max(abs(Pd(:) - M.dic.Phi_dic(:)));
        chk(d < 1e-12, 'DIC modal basis differs from the original by %g', d);
        fprintf('    DIC basis matches the original driver (max diff %.2e)\n', d);
    end

    % ---------- preflight ---------------------------------------------------
    F = qs_preflight(M, C);

    % The aero file for this case is the known COLD-FLOW RANS (implied
    % T0 ~ 291 K) while the case ran at 388 K. Preflight exists to catch
    % exactly this, so a PASS here would mean the gate is broken.
    chk(~F.checks.T0.pass, ...
        ['preflight did NOT flag the known cold-flow aero file. Implied T0 came ' ...
         'out as %.1f K against an expected %.1f K; the check is not working.'], ...
        F.checks.T0.mean, F.checks.T0.expected);
    if ~F.checks.T0.pass
        fprintf('    preflight correctly flagged the cold-flow aero: T0 %.1f K vs %.1f K expected\n', ...
            F.checks.T0.mean, F.checks.T0.expected);
    end

    chk(isfinite(F.checks.buckling.dT_cr) && F.checks.buckling.dT_cr > 0, ...
        'dT_cr did not compute: %g', F.checks.buckling.dT_cr);
    if isfinite(F.checks.buckling.dT_cr)
        % recorded value from the earlier analysis
        chk(abs(F.checks.buckling.dT_cr - 5.65) < 0.5, ...
            'dT_cr is %.2f K; the recorded value for this ROM is 5.65 K', F.checks.buckling.dT_cr);
        fprintf('    dT_cr = %.2f K (recorded: 5.65 K),  f1 = %.1f Hz\n', ...
            F.checks.buckling.dT_cr, F.checks.buckling.f1);
    end
end
