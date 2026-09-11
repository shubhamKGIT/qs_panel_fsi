function [Fnl, J] = rc19_nonlinear_force(q, rom)
%RC19_NONLINEAR_FORCE  Geometric nonlinear modal restoring force (and Jacobian).
%   Fnl       = RC19_NONLINEAR_FORCE(q, rom) returns the quadratic+cubic
%               modal force for modal displacement q.
%   [Fnl, J]  = RC19_NONLINEAR_FORCE(q, rom) also returns the tangent
%               stiffness  J = dFnl/dq  (used by the static Newton solver).
%
%     Fnl(i) = sum_{quad}  K2(i,j,m)   * q(j)*q(m)
%            + sum_{cube}  K3(i,j,m,n) * q(j)*q(m)*q(n)
%
%   Evaluation is fully vectorised via ACCUMARRAY: each stored coefficient
%   contributes once, scattered into the i-th modal equation. No dense
%   3-/4-D tensors are ever built, so this is fast enough to call inside
%   the ode45 right-hand side every step.

    nm = rom.num_modes;
    nl = rom.nl;

    % --- nonlinear force ---------------------------------------------
    Fq = accumarray(nl.qi, nl.kq .* q(nl.qj) .* q(nl.qm),             [nm,1]);
    Fc = accumarray(nl.ci, nl.kc .* q(nl.cj) .* q(nl.cm) .* q(nl.cn), [nm,1]);
    Fnl = Fq + Fc;

    % --- tangent stiffness J = dFnl/dq (only if requested) -----------
    if nargout > 1
        J = zeros(nm,nm);
        % quadratic:  d/dq_a [ K2 q_j q_m ]
        J = J + accumarray([nl.qi nl.qj], nl.kq .* q(nl.qm), [nm nm]);
        J = J + accumarray([nl.qi nl.qm], nl.kq .* q(nl.qj), [nm nm]);
        % cubic:      d/dq_a [ K3 q_j q_m q_n ]
        J = J + accumarray([nl.ci nl.cj], nl.kc .* q(nl.cm) .* q(nl.cn), [nm nm]);
        J = J + accumarray([nl.ci nl.cm], nl.kc .* q(nl.cj) .* q(nl.cn), [nm nm]);
        J = J + accumarray([nl.ci nl.cn], nl.kc .* q(nl.cj) .* q(nl.cm), [nm nm]);
    end
end
