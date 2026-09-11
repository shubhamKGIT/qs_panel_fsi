function [pnet, p_total, p_qs, p_tbl] = pm_quasisteady_ept(state, aero)
%PM_QUASISTEADY_EPT  Mean + quasi-steady EPT (+ optional TBL) pressure model.
%   This is the pluggable pressure-model HOOK for rc19_solve_structural. It
%   implements the fluid model
%
%       p(x,y,t) = P_rans(x,y) + P_PT(x,y,t) + P_TBL(x,y,t)
%       P_rans   = aero.pl                      (fixed RANS SURFACE pressure)
%       P_PT     = gamma*p_edge*[ c1*(vn/a_l) + c2*(vn/a_l)^2 + c3*(vn/a_l)^3 ]
%       vn       = U_l*Z' + Zdot                ( Z = aero.zsign * w )
%       P_TBL    = aero.tblModel(state)         (turbulent-BL load; HOOK, default 0)
%       Van Dyke 2nd-order => c3 = 0 (quadratic perturbation). The "c2" on the
%       cubic term in the paper's Eq.6 is a misprint (should be c3 = 0).
%
%   P_PT is built ENTIRELY from local 7 mm boundary-layer-edge RANS
%   quantities: p_edge, M_l (-> c1,c2), a_l, and U_l = M_l*a_l.
%
%   using the LOCAL boundary-layer-edge conditions pre-computed in `aero`
%   (see cm_setup_aero). The instantaneous structural deformation is read
%   from `state`, closing the fluid-structure coupling loop.
%
%   OUTPUTS
%     pnet    : NET nodal load delivered to the structure  [Pa] (nx*ny x 1)
%               = (p_l + p_QS) - aero.pback
%               This is what rc19_solve_structural projects with TransferMat.
%     p_total : TOTAL instantaneous surface pressure P_rans+P_PT+P_TBL [Pa]
%               (use this for post-processing / writing out)
%     p_qs    : quasi-steady perturbation P_PT alone           [Pa]
%     p_tbl   : turbulent-BL fluctuating load P_TBL alone      [Pa] (0 if no hook)
%
%   NOTE on signs / orientation:
%     Z = aero.zsign*w relates the structural deflection convention to the
%     piston-theory deflection (positive toward the flow). aero.zsign = -1
%     makes the model dissipative (physically correct for single-surface
%     piston theory). The mean load (p_l - pback) uses the same TransferMat
%     projection as the validated static load, so its sign is consistent.

    nx = state.nx;  ny = state.ny;

    % structural deflection / velocity in the piston-theory convention
    Z  = reshape(aero.zsign*state.w,    nx, ny);   % [m]
    Zt = reshape(aero.zsign*state.wdot, nx, ny);   % [m/s]

    % streamwise slope Z' = dZ/dx  (x varies along dim 1; flow in +x)
    dx          = state.x_grid(2,1) - state.x_grid(1,1);
    dZdx        = zeros(nx,ny);
    dZdx(2:end-1,:) = (Z(3:end,:) - Z(1:end-2,:)) / (2*dx);
    dZdx(1,:)       = (Z(2,:)     - Z(1,:))       / dx;
    dZdx(end,:)     = (Z(end,:)   - Z(end-1,:))   / dx;
    dZdx = dZdx(:);

    % piston (downwash) velocity  vn = U_l*Z' + Zdot   (Eq.3)
    vn = aero.Ul .* dZdx + Zt(:);

    % quasi-steady enriched piston-theory pressure  (Eq.6 of Brouwer 2023)
    % NOTE 1: P_PT is built ENTIRELY from local 7 mm BL-edge RANS quantities,
    %         so the leading pressure here is the EDGE pressure p_edge
    %         (NOT the RANS surface pressure, which is the fixed P_rans term).
    % NOTE 2: this is Van Dyke SECOND-order theory => it is a QUADRATIC
    %         perturbation and c3 = 0 (the cubic term vanishes). The published
    %         Eq.6 prints "c2" on the cubic term, but that is a MISPRINT:
    %         Eq.8 defines c3 = 0, the text calls it a "quadratic perturbation",
    %         and a 2nd-order theory has no cubic term (verified by derivation).
    %         To use a true 3rd-order CLASSICAL cubic instead, set
    %         aero.c3 = (gamma+1)/12 in cm_setup_aero.
    s    = vn ./ aero.al;
    p_qs = aero.gamma .* aero.p_edge .* ( aero.c1.*s + aero.c2.*s.^2 + aero.c3.*s.^3 );   % = P_PT

    % ---- turbulent-boundary-layer fluctuating pressure  P_TBL  (HOOK) ----
    %  Pressure decomposition (Brouwer 2023, Eq.1):  p = P_rans + P_PT + P_TBL.
    %  P_TBL is the fluctuating load from the turbulent BL / unsteady SBLI and
    %  is what sustains the measured limit-cycle oscillation. It is NOT modelled
    %  yet -- this is the plug-in point. Provide aero.tblModel as a function
    %  handle  p_tbl = tblModel(state)  returning nodal pressure [Pa] (nx*ny x 1).
    %  Default (empty) => P_TBL = 0, i.e. the current quasi-steady model.
    if isfield(aero,'tblModel') && ~isempty(aero.tblModel)
        p_tbl = aero.tblModel(state);
        p_tbl = p_tbl(:);
    else
        p_tbl = zeros(nx*ny, 1);
    end

    % total surface pressure  =  P_rans + P_PT + P_TBL
    p_total = aero.pl + p_qs + p_tbl;
    pnet    = p_total - aero.pback;   % net structural load (back/cavity offset)
end
