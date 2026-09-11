function out = rc19_solve_structural(rom, cfg)
%RC19_SOLVE_STRUCTURAL  Time-integrate the RC-19 nonlinear structural ROM.
%
%   out = RC19_SOLVE_STRUCTURAL(rom, cfg) integrates
%
%       M*qddot + C*qdot + (Ke - dT*Kt)*q + Fnl(q) = F_aero(t,q,qdot) + dT*Ft
%
%   and returns the modal/physical response. The aerodynamic (pressure)
%   load can be either a fixed prescribed load OR supplied every step by a
%   pluggable PRESSURE MODEL, which is the hook for fluid-structure coupling.
%
%   --------------------------------------------------------------------
%   cfg fields (all optional unless noted)
%     .dT             temperature rise [K]                  (default 0)
%     .tspan          output time vector [s]                (default 0:2e-5:0.05)
%     .y0             initial state [q; qdot] (2*nm x 1)    (default zeros)
%     .integrator     ODE integrator handle @ode45/@ode15s  (default @ode45)
%     .odeopts        odeset(...) options                   (default RelTol 1e-6)
%     .reconstruct    reconstruct physical w(x,y,t)?        (default true)
%
%   LOAD SPECIFICATION  (choose ONE)
%     .staticPressure pressure load held constant in time:
%                       - scalar  -> uniform over mesh [Pa]
%                       - (nx*ny) vector -> spatial field [Pa]
%     .pressureModel  function handle  p = f(state)  returning the
%                     INSTANTANEOUS nodal pressure [Pa] (nx*ny x 1).
%                     When set, this is called every RHS evaluation and the
%                     current structural deformation is passed back to it
%                     through `state` (see COUPLING INTERFACE below). This
%                     enables two-way FSI coupling.
%     .staticPressure may also be combined with .pressureModel, in which
%     case the model output is ADDED to the steady prescribed load.
%
%   COUPLING INTERFACE  -  the struct `state` passed to pressureModel:
%     state.t        current time [s]
%     state.q        modal displacement      (nm x 1)
%     state.qdot     modal velocity          (nm x 1)
%     state.w        nodal deflection w      (nx*ny x 1) [m]
%     state.wdot     nodal velocity dw/dt    (nx*ny x 1) [m/s]
%     state.nx,.ny   grid size
%     state.x_grid,.y_grid   node coordinates (nx x ny) [m]
%   The model returns instantaneous nodal pressure (nx*ny x 1) [Pa].
%
%   --------------------------------------------------------------------
%   out fields
%     .T      time vector
%     .q      modal displacement history     (Nt x nm)
%     .qdot   modal velocity history         (Nt x nm)
%     .w      physical deflection w(x,y,t)   (nx x ny x Nt) [m]  (if reconstruct)
%     .cfg    the resolved configuration used

    % ---------- resolve defaults ------------------------------------
    if ~isfield(cfg,'dT')          || isempty(cfg.dT),          cfg.dT = 0;                 end
    if ~isfield(cfg,'tspan')       || isempty(cfg.tspan),       cfg.tspan = 0:2e-5:0.05;    end
    if ~isfield(cfg,'integrator')  || isempty(cfg.integrator),  cfg.integrator = @ode45;    end
    if ~isfield(cfg,'odeopts')     || isempty(cfg.odeopts)
        cfg.odeopts = odeset('RelTol',1e-6,'AbsTol',1e-8);
    end
    if ~isfield(cfg,'reconstruct') || isempty(cfg.reconstruct), cfg.reconstruct = true;     end
    if ~isfield(cfg,'staticPressure'), cfg.staticPressure = [];  end
    if ~isfield(cfg,'pressureModel'),  cfg.pressureModel  = [];  end

    nm = rom.num_modes;

    % ---------- constant operators (built ONCE) ---------------------
    Keff = rom.Ke - cfg.dT*rom.Kt;     % thermal pre-stress softens stiffness
    Fth  = cfg.dT*rom.Ft;              % constant thermal modal force
    Minv = rom.Minv;
    Cmat = rom.C;

    % steady (prescribed) aerodynamic load, projected to modal space once
    if isempty(cfg.staticPressure)
        Faero_const = zeros(nm,1);
    else
        Faero_const = rom.TransferMat * pressure_to_nodal(cfg.staticPressure, rom);
    end

    useModel = ~isempty(cfg.pressureModel);   % coupling on/off flag

    % ---------- initial state ---------------------------------------
    if ~isfield(cfg,'y0') || isempty(cfg.y0)
        y0 = zeros(2*nm,1);
    else
        y0 = cfg.y0(:);
    end

    % ---------- integrate -------------------------------------------
    [T,Y] = cfg.integrator(@rhs, cfg.tspan, y0, cfg.odeopts);

    % ---------- pack results ----------------------------------------
    out.T    = T;
    out.q    = Y(:,1:nm);
    out.qdot = Y(:,nm+1:2*nm);
    out.cfg  = cfg;
    if cfg.reconstruct
        out.w = rc19_reconstruct(out.q, rom);
    end

    % =================================================================
    %  Nested RHS  -  shares the precomputed operators by closure
    %  (no large data copied per call => efficient).
    % =================================================================
    function dy = rhs(t,y)
        q    = y(1:nm);
        qdot = y(nm+1:2*nm);

        % --- aerodynamic / pressure load -----------------------------
        Faero = Faero_const;
        if useModel
            % Reconstruct the instantaneous physical deformation and pass
            % it back to the pressure model  (the FSI coupling step).
            state.t     = t;
            state.q     = q;       state.qdot = qdot;
            state.w     = rom.Modes.' * q;        % nodal deflection [m]
            state.wdot  = rom.Modes.' * qdot;     % nodal velocity [m/s]
            state.nx    = rom.nx;  state.ny = rom.ny;
            state.x_grid = rom.x_grid;  state.y_grid = rom.y_grid;

            p     = cfg.pressureModel(state);     % instantaneous nodal pressure
            Faero = Faero + rom.TransferMat * p(:);
        end

        % --- nonlinear restoring force -------------------------------
        Fnl = rc19_nonlinear_force(q, rom);

        % --- state derivative ----------------------------------------
        dy            = zeros(2*nm,1);
        dy(1:nm)      = qdot;
        dy(nm+1:2*nm) = Minv*(Faero + Fth - Cmat*qdot - Keff*q - Fnl);
    end
end

% ---------------------------------------------------------------------
function pv = pressure_to_nodal(P, rom)
%PRESSURE_TO_NODAL  Expand a scalar (uniform) or accept a nodal field.
    if isscalar(P)
        pv = P*ones(rom.nx*rom.ny,1);
    else
        pv = P(:);
        assert(numel(pv)==rom.nx*rom.ny, ...
            'staticPressure must be scalar or length nx*ny = %d', rom.nx*rom.ny);
    end
end
