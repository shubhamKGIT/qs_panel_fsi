function w = rc19_reconstruct(q_hist, rom)
%RC19_RECONSTRUCT  Modal coordinates -> physical panel deflection.
%   w = RC19_RECONSTRUCT(q_hist, rom) maps the modal displacement history
%   q_hist (Nt x nm) to the physical out-of-plane deflection field
%   w(x,y,t), returned as an (nx x ny x Nt) array [m]:
%
%       w(:,:,k) = sum_j phi_j(x,y) * q_hist(k,j)
%
%   Vectorised as a single matrix product for speed.

    Nt     = size(q_hist,1);
    Wnodes = rom.Modes.' * q_hist.';          % (nx*ny) x Nt
    w      = reshape(Wnodes, rom.nx, rom.ny, Nt);
end
