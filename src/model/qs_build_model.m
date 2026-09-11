function M = qs_build_model(C, quiet)
%QS_BUILD_MODEL  Assemble everything a run needs, once, from the config.
%
%   M = QS_BUILD_MODEL(C)   with C from qs_config.
%
%   This is the single place where case data turns into model objects. It
%   replaces the block that used to sit at the top of every copied driver:
%   ROM load, aero setup, DIC reference and the modal basis interpolated to
%   the DIC points. Nothing here depends on the operating point (p_c, dT), so
%   it is built once per process and reused for every point in the slice.
%
%   OUTPUT (struct M)
%     .rom        from rc19_setup_rom
%     .aero       from cm_setup_aero   (aero.pback is overwritten per point)
%     .nm         number of modes
%     .Phi        rom.Modes.'          (nodes x modes)
%     .cidx       panel-centre node index
%     .h          panel thickness [m]   (sets the w/h amplitude scale)
%     .dic        struct: X, Y, mean, std, Phi_dic, nrm_mean, nrm_std, present
%
%   If the reference (DIC) file is absent, M.dic.present = false and the
%   DIC-comparison metrics come back NaN. Everything else still runs, so a
%   case can be swept before its measurement file has been dropped in.

    if nargin < 2, quiet = false; end

    % ---------- structural ROM -------------------------------------------
    M.rom = rc19_setup_rom(C.paths.rom_file);
    M.nm  = M.rom.num_modes;
    M.Phi = M.rom.Modes.';
    M.cidx = round(M.rom.nx/2) + (round(M.rom.ny/2)-1)*M.rom.nx;   % panel centre
    M.h    = C.case.panel.h_m;

    % ---------- aerodynamics ---------------------------------------------
    opts = struct( ...
        'aerofile', C.paths.aero_file, ...
        'pinf',     C.case.pinf_Pa, ...
        'Minf',     C.case.Minf, ...
        'gamma',    C.aero.gamma, ...
        'Lp',       C.case.panel.L_m, ...
        'pback',    C.case.pc_nom_Pa, ...   % per-point value overrides this
        'zsign',    C.aero.zsign);
    % evalc is used only to swallow cm_setup_aero's banner in quiet mode.
    % Octave's evalc does not forward the expression's output arguments, so
    % there the banner is simply printed -- harmless, and Octave is only ever
    % used for checks here.
    if quiet && ~qs_isoctave()
        [~, M.aero] = evalc('cm_setup_aero(M.rom, opts)');
    else
        M.aero = cm_setup_aero(M.rom, opts);
    end

    % Optional: true 3rd-order classical cubic instead of Van Dyke 2nd order.
    if isfield(C.aero,'c3_mode') && strcmpi(C.aero.c3_mode,'classical_cubic')
        M.aero.c3 = ((C.aero.gamma+1)/12) * ones(size(M.aero.Ml));
        if ~quiet, fprintf('qs_build_model: c3 set to classical cubic (gamma+1)/12\n'); end
    end

    % ---------- DIC reference + modal basis at the DIC points -------------
    M.dic = build_dic(C, M);

    if ~quiet
        fprintf('qs_build_model: %d modes, grid %dx%d, DIC %s\n', ...
            M.nm, M.rom.nx, M.rom.ny, tern(M.dic.present, ...
            sprintf('%d points', numel(M.dic.X)), 'ABSENT (RMSE metrics -> NaN)'));
    end
end

% ======================================================================
function D = build_dic(C, M)
    D.present = false;
    D.X = []; D.Y = []; D.mean = []; D.std = []; D.Phi_dic = [];
    D.nrm_mean = NaN; D.nrm_std = NaN;

    if exist(C.paths.ref_file,'file')~=2
        return
    end
    R = load(C.paths.ref_file);
    if ~isfield(R,'DIC')
        warning('qs_build_model:noDIC', ...
            '%s has no DIC struct; RMSE metrics will be NaN.', C.paths.ref_file);
        return
    end

    D.X    = R.DIC.X_m;
    D.Y    = R.DIC.Y_m;
    D.mean = R.DIC.MeanDef_m;
    D.std  = R.DIC.StdDef_m;

    % The DIC grid is centred on y = 0; the ROM grid runs 0..W. The offset is
    % a case setting (dic_y_offset_m), defaulting to half the panel width --
    % the same 0.5*0.127 the original drivers hardcoded.
    yoff = 0.5 * C.case.panel.W_m;
    if isfield(C.case,'dic_y_offset_m') && ~isempty(C.case.dic_y_offset_m)
        yoff = C.case.dic_y_offset_m;
    end

    xv = M.rom.x_grid(:,1);
    yv = M.rom.y_grid(1,:).';
    xq = D.X(:);
    yq = D.Y(:) + yoff;

    D.Phi_dic = zeros(numel(xq), M.nm);
    for j = 1:M.nm
        Fj = griddedInterpolant({xv,yv}, ...
             reshape(M.rom.Modes(j,:), M.rom.nx, M.rom.ny), 'linear', 'nearest');
        D.Phi_dic(:,j) = Fj(xq, yq);
    end

    D.nrm_mean = norm(D.mean(:));
    D.nrm_std  = norm(D.std(:));
    D.present  = true;
end

function s = tern(c,a,b)
    if c, s = a; else, s = b; end
end
