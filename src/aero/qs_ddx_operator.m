function Dx = qs_ddx_operator(rom)
%QS_DDX_OPERATOR  Sparse d/dx operator on the structural grid (nx*ny x nx*ny).
%   Dx*f returns the streamwise derivative of a nodal field f, using the SAME
%   stencil as pm_quasisteady_ept (2nd-order central interior, 1st-order
%   one-sided at the x edges; x varies along dim 1, flow in +x). Built once so
%   the instantaneous P_QS field can be reconstructed for a whole chunk of
%   time steps at a time:  dZdx = Dx*Z  with Z of size (nx*ny x Nsteps).
%
%   Identical in behaviour to cm_ddx_operator; vectorised assembly so a
%   101x51 grid builds in milliseconds instead of a double loop.

    nx = rom.nx;  ny = rom.ny;
    dx = rom.x_grid(2,1) - rom.x_grid(1,1);
    n  = nx*ny;

    i = (1:nx).';
    interior = i(2:end-1);

    I = []; J = []; V = [];
    for j = 1:ny
        off = (j-1)*nx;
        % interior: central difference
        r = interior + off;
        I = [I; r; r];                                          %#ok<AGROW>
        J = [J; interior+1+off; interior-1+off];                %#ok<AGROW>
        V = [V; repmat(1/(2*dx),numel(r),1); repmat(-1/(2*dx),numel(r),1)]; %#ok<AGROW>
        % leading edge: forward difference
        I = [I; 1+off; 1+off];  J = [J; 2+off; 1+off];  V = [V; 1/dx; -1/dx];   %#ok<AGROW>
        % trailing edge: backward difference
        I = [I; nx+off; nx+off]; J = [J; nx+off; nx-1+off]; V = [V; 1/dx; -1/dx]; %#ok<AGROW>
    end
    Dx = sparse(I, J, V, n, n);
end
