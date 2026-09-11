function rom = rc19_setup_rom(matfile)
%RC19_SETUP_ROM  Load and pre-process the RC-19 structural ROM.
%   rom = RC19_SETUP_ROM(matfile) reads the reduced-order-model data file
%   and returns a struct holding every quantity the solver needs, with all
%   time-independent pre-processing done ONCE here (so the time-integration
%   hot loop stays cheap).
%
%   INPUT
%     matfile : ROM .mat file name (default 'RC19_ROM_v1.mat')
%
%   OUTPUT (struct rom)
%     .num_modes        number of modes (=18)
%     .nx, .ny          structural grid size (101 x 51)
%     .x_grid, .y_grid  2D node coordinates [m]            (nx x ny)
%     .M, .C            modal mass / damping matrices      (nm x nm)
%     .Ke               linear elastic stiffness           (nm x nm)
%     .Kt               thermal stiffness per 1 K          (nm x nm)
%     .Ft               thermal force per 1 K              (nm x 1)
%     .Minv             inverse mass matrix (precomputed)  (nm x nm)
%     .Modes            mode shapes, nm x (nx*ny)
%     .TransferMat      load projection, nm x (nx*ny)
%     .nl               struct of nonlinear index/value arrays (see below)
%
%   The geometric nonlinear stiffness is stored sparsely as a list of
%   (i,j,m,n) indices with matching coefficients. They are split here into
%   quadratic and cubic blocks so that rc19_nonlinear_force can evaluate
%   Fnl(q) by vectorised accumulation without ever forming dense tensors.

    if nargin < 1 || isempty(matfile)
        matfile = 'RC19_ROM_v1.mat';
    end

    S = load(matfile);

    rom.num_modes = size(S.Modes,1);
    rom.nx = 101;  rom.ny = 51;

    % structural grid
    rom.x_grid = reshape(S.xy_rom(1,:)', rom.nx, rom.ny);
    rom.y_grid = reshape(S.xy_rom(2,:)', rom.nx, rom.ny);

    % linear operators
    rom.M    = S.M_rom;
    rom.C    = S.C_rom;
    rom.Ke   = S.K_rom;
    rom.Kt   = S.Ktemp_rom;
    rom.Ft   = S.Ftemp_rom(:);
    rom.Minv = inv(rom.M);          % constant -> invert once

    % spatial operators
    rom.Modes       = S.Modes;          % nm x (nx*ny)
    rom.TransferMat = S.TransferMat;    % nm x (nx*ny)

    % ---- nonlinear stiffness index/value arrays --------------------
    Knl = S.Knl_rom(:);
    nq  = round(S.Nquad_rom);           % number of quadratic terms (3078)

    qd = S.Indxnl_rom(1:nq,     :);     % quadratic rows  [i j m 0]
    cu = S.Indxnl_rom(nq+1:end, :);     % cubic rows      [i j m n]

    rom.nl.qi = qd(:,1);  rom.nl.qj = qd(:,2);  rom.nl.qm = qd(:,3);
    rom.nl.kq = Knl(1:nq);
    rom.nl.ci = cu(:,1);  rom.nl.cj = cu(:,2);
    rom.nl.cm = cu(:,3);  rom.nl.cn = cu(:,4);
    rom.nl.kc = Knl(nq+1:end);
end
