function aero = cm_setup_aero(rom, opts)
%CM_SETUP_AERO  Assemble the aerodynamic data for the quasi-steady coupled model.
%
%   aero = CM_SETUP_AERO(rom, opts) LOADS the pre-interpolated structural-grid
%   aerodynamic data from a .mat file (default: rc19_aero_data.mat, located
%   next to this function) and builds the struct the EPT pressure model needs.
%   rc19_aero_data.mat is provided pre-computed (RANS surface pressure + 7 mm
%   boundary-layer-edge conditions already interpolated onto the structural grid).
%
%   The implemented fluid model (no turbulent-BL term) is
%
%       p(x,y,t) = P_rans(x,y)  +  P_PT(x,y,t)                       (Eq.1, TBL omitted)
%       p_QS = gamma * p_edge * [ c1*(vn/a_l) + c2*(vn/a_l)^2 + c3*(vn/a_l)^3 ]   (Eq.6)
%       vn   = U_l * Z' + Zdot                                       (Eq.7)
%       c1   = M_l / sqrt(M_l^2 - 1)
%       c2   = ( M_l^4*(gamma+1) - 4*(M_l^2-1) ) / ( 4*(M_l^2-1)^2 ) (Eq.8)
%       c3   = 0     (Van Dyke 2nd-order => quadratic perturbation)
%   NOTE: the published Eq.6 prints "c2" on the cubic term, but that is a
%   misprint -- Eq.8 sets c3 = 0 and the text describes a quadratic
%   perturbation (confirmed by deriving the piston-theory series). For a
%   true 3rd-order CLASSICAL cubic, set aero.c3 = (gamma+1)/12 below.
%
%   INPUTS
%     rom   : struct from rc19_setup_rom (used to verify the grid matches)
%     opts  : (optional) struct, fields:
%               .aerofile  path to rc19_aero_data.mat (default: next to this file)
%               .pinf      reference far-field pressure [Pa]   (default: from file)
%               .Minf      incoming Mach number                (default 1.92)
%               .gamma     ratio of specific heats             (default 1.4)
%               .Lp        panel length [m]                    (default: from file)
%               .pback     back / cavity reference pressure [Pa] (default = pinf)
%               .zsign     orientation of structural w vs piston Z (Z = zsign*w);
%                          zsign = -1 (default) makes single-surface piston
%                          theory dissipative (physically correct).
%
%   OUTPUT (struct aero), all nodal arrays (nx*ny x 1), node-ordered like rom.Modes:
%     .pl     mean RANS SURFACE pressure  (P_rans, fixed)            [Pa]
%     .p_edge local BL-EDGE static pressure (scales P_PT)            [Pa]
%     .Ml, .al, .Ul   local BL-edge Mach, sound speed, velocity
%     .c1,.c2,.c3     Van Dyke coefficients (c3 = 0)
%     .gamma,.pinf,.Minf,.pback,.zsign,.Lp

    if nargin < 2, opts = struct(); end
    here = fileparts(mfilename('fullpath'));
    if ~isfield(opts,'aerofile') || isempty(opts.aerofile)
        opts.aerofile = fullfile(here,'rc19_aero_data.mat');
    end
    if ~isfield(opts,'gamma'), opts.gamma = 1.4;  end
    if ~isfield(opts,'Minf'),  opts.Minf  = 1.92; end

    if ~isfile(opts.aerofile)
        error(['cm_setup_aero:missingData','Aero data file not found:\n  %s\n', ...
               'It should be provided alongside this function (rc19_aero_data.mat).'], opts.aerofile);
    end
    D = load(opts.aerofile);

    % grid-consistency check (node ordering must match the structural ROM)
    assert(D.nx==rom.nx && D.ny==rom.ny, ...
        'aero .mat grid (%dx%d) does not match ROM grid (%dx%d).', ...
        D.nx, D.ny, rom.nx, rom.ny);

    % --- local edge / surface fields (already on the structural grid) ----
    aero.pl     = D.pl_surface(:);        % P_rans surface pressure [Pa]
    aero.p_edge = D.p_edge(:);            % BL-edge static pressure  [Pa]
    aero.Ml     = D.M_edge(:);            % BL-edge Mach number      [-]
    aero.al     = D.a_edge(:);            % BL-edge sound speed      [m/s]
    aero.Ul     = aero.Ml .* aero.al;     % BL-edge velocity         [m/s]

    % --- Van Dyke second-order coefficients (Eq.8) -----------------------
    g  = opts.gamma;  M2 = aero.Ml.^2;
    aero.c1 = aero.Ml ./ sqrt(M2 - 1);
    aero.c2 = ( aero.Ml.^4*(g+1) - 4*(M2 - 1) ) ./ ( 4*(M2 - 1).^2 );
    aero.c3 = zeros(size(aero.Ml));       % cubic term omitted (see NOTE above)

    % --- carry constants -------------------------------------------------
    if ~isfield(opts,'pinf')  || isempty(opts.pinf),  opts.pinf  = D.pinf; end
    if ~isfield(opts,'Lp')    || isempty(opts.Lp),    opts.Lp    = D.Lp;   end
    if ~isfield(opts,'pback') || isempty(opts.pback), opts.pback = opts.pinf; end
    if ~isfield(opts,'zsign') || isempty(opts.zsign), opts.zsign = -1;     end
    aero.gamma = g;          aero.pinf  = opts.pinf;  aero.Minf  = opts.Minf;
    aero.pback = opts.pback; aero.zsign = opts.zsign; aero.Lp    = opts.Lp;

    % turbulent-BL fluctuating-pressure hook (P_TBL). Empty => P_TBL = 0.
    % Set aero.tblModel = @(state) my_tbl(state,...) to plug a model in later
    % (see pm_tbl_template.m). The pressure model adds its output to the load.
    if isfield(opts,'tblModel'), aero.tblModel = opts.tblModel; else, aero.tblModel = []; end

    fprintf(['cm_setup_aero: loaded %s\n  ', ...
             'P_rans = %.1f..%.1f kPa | p_edge = %.1f..%.1f kPa | ', ...
             'M_l = %.2f..%.2f | a_l = %.0f..%.0f m/s | U_l = %.0f..%.0f m/s\n'], ...
        opts.aerofile, min(aero.pl)/1e3,max(aero.pl)/1e3, ...
        min(aero.p_edge)/1e3,max(aero.p_edge)/1e3, ...
        min(aero.Ml),max(aero.Ml),min(aero.al),max(aero.al), ...
        min(aero.Ul),max(aero.Ul));
end
