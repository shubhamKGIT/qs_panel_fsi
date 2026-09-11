function p_tbl = pm_tbl_template(state, tbl)
%PM_TBL_TEMPLATE  Template / HOOK for the turbulent-boundary-layer load P_TBL.
%
%   p_tbl = PM_TBL_TEMPLATE(state, tbl) returns the instantaneous fluctuating
%   turbulent-boundary-layer (and unsteady-SBLI) surface pressure on the panel
%   as a nodal vector [Pa] (nx*ny x 1), node-ordered like rom.Modes columns.
%
%   This is the third term of the pressure decomposition (Brouwer 2023, Eq.1):
%       p = P_rans + P_PT + P_TBL
%   P_TBL is the load that SUSTAINS the measured limit-cycle oscillation. It is
%   NOT implemented yet -- this template just returns zeros and documents the
%   interface so a real model can be dropped in later.
%
%   WIRE IT UP (in cm_main_coupled_RC19.m or cm_compare_reference.m), e.g.:
%       tbl  = struct(...);                          % your TBL parameters
%       aero.tblModel = @(s) pm_tbl_template(s, tbl);% attach the hook
%   The EPT pressure model (pm_quasisteady_ept) will then add p_tbl to the load.
%
%   `state` provides everything a TBL model might need:
%       state.t                 current time [s]
%       state.q, state.qdot     modal displacement / velocity
%       state.w, state.wdot     nodal deflection [m] / velocity [m/s]
%       state.nx, state.ny      grid size
%       state.x_grid, state.y_grid   node coordinates (nx x ny) [m]
%
%   IMPLEMENTATION NOTES (for the future task):
%     - In the experiment p_TBL was computed for a RIGID panel and applied
%       UNCOUPLED (independent of the structural response). A first model can
%       therefore ignore state.w/state.wdot and depend only on (x, y, t).
%     - Typical choices: a band-limited random pressure field with the measured
%       wall-pressure PSD/coherence, or a time series interpolated from a
%       precomputed (x, y, t) database, or a Corcos/empirical TBL model.
%     - Return a column vector of length state.nx*state.ny in [Pa].

    % ---- PLACEHOLDER: no TBL load yet ----
    p_tbl = zeros(state.nx * state.ny, 1);

    % ---- EXAMPLE skeleton (commented) ------------------------------------
    % x = state.x_grid(:);  y = state.y_grid(:);  t = state.t;
    % p_tbl = tbl.amplitude * sin(2*pi*tbl.freq*t) .* tbl.shape(x,y);   % e.g.
end
