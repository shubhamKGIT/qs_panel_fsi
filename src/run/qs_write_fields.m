function qs_write_fields(fdir, M, C, tvec, q, qdot, Fb)
%QS_WRITE_FIELDS  Stream instantaneous surface pressure and deflection to disk.
%
%   qs_write_fields(fdir, M, C, tvec, q, qdot, Fb)
%
%   Reconstructs, for every saved step,
%       P_surf = P_rans + P_QS      (nodal, structural grid) [Pa]
%       W      = panel deflection   (nodal, structural grid) [m]
%   from the modal history, in chunks, so memory stays bounded regardless of
%   run length.
%
%   Fb (the study's "fields" block):
%     .mode       'per_second' (default) | 'per_step' | 'single'
%     .precision  'single' (default) | 'double'
%     .stride     save every Nth step (default 1)
%
%   'per_second' writes fields_sec_001.mat, fields_sec_002.mat, ... which is
%   what makes a 20 s run survivable: a job that dies at second 14 leaves 13
%   usable files behind.

    mode   = getdef(Fb,'mode','per_second');
    prec   = getdef(Fb,'precision','single');
    stride = getdef(Fb,'stride',1);

    nn  = M.rom.nx * M.rom.ny;
    Nt  = numel(tvec);
    sel = 1:stride:Nt;
    spc = max(1, round(1/C.duration.dt_out));      % steps per simulated second
    chunk = 20000;                                 % steps reconstructed at once

    Dx    = qs_ddx_operator(M.rom);
    cast2 = @(x) cast(x, prec);
    recon = @(Qc,Qdc) recon_surface(M.aero, Dx, M.rom, Qc, Qdc);

    switch lower(mode)
      case 'single'
        if qs_isoctave()
            error('qs_write_fields: mode "single" needs MATLAB matfile(); use "per_second" here.');
        end
        mfP = matfile(fullfile(fdir,'pressure_all.mat'),'Writable',true);
        mfW = matfile(fullfile(fdir,'deform_all.mat'),  'Writable',true);
        mfP.P_surf = zeros(nn, numel(sel), prec);
        mfW.W      = zeros(nn, numel(sel), prec);
        mfP.tvec = tvec(sel);  mfW.tvec = tvec(sel);
        for a = 1:chunk:numel(sel)
            b  = min(a+chunk-1, numel(sel));
            ks = sel(a:b);
            [Ps, Wf] = recon(q(ks,:).', qdot(ks,:).');
            mfP.P_surf(:, a:b) = cast2(Ps);
            mfW.W(:, a:b)      = cast2(Wf);
        end

      case 'per_step'
        for n = sel
            [Ps, Wf] = recon(q(n,:).', qdot(n,:).');
            P_surf = cast2(Ps);  W = cast2(Wf);  t = tvec(n); %#ok<NASGU>
            save(fullfile(fdir,sprintf('step_%08d.mat',n)),'P_surf','W','t',qs_matver());
        end

      otherwise      % per_second
        nsec = ceil(Nt/spc);
        for s = 1:nsec
            ks = (s-1)*spc+1 : min(s*spc, Nt);
            ks = ks(mod(ks-1,stride)==0);
            if isempty(ks), continue; end
            [Ps, Wf] = recon(q(ks,:).', qdot(ks,:).');
            P_surf = cast2(Ps);  W = cast2(Wf);  t = tvec(ks); %#ok<NASGU>
            save(fullfile(fdir,sprintf('fields_sec_%03d.mat',s)),'P_surf','W','t',qs_matver());
        end
    end
end

% ======================================================================
function [Psurf, W] = recon_surface(aero, Dx, rom, Qc, Qdc)
%RECON_SURFACE  P_RANS + P_QS and deflection for a chunk of modal states.
%   Same formula as pm_quasisteady_ept, vectorised across the chunk via the
%   precomputed d/dx operator instead of the per-call finite difference.
    W  = rom.Modes.' * Qc;                       % nodal deflection [m]
    Wd = rom.Modes.' * Qdc;                      % nodal velocity  [m/s]
    Z  = aero.zsign * W;
    Zt = aero.zsign * Wd;
    s  = (aero.Ul .* (Dx*Z) + Zt) ./ aero.al;
    Pqs   = aero.gamma * aero.p_edge .* (aero.c1.*s + aero.c2.*s.^2 + aero.c3.*s.^3);
    Psurf = aero.pl + Pqs;
end

function v = getdef(S, f, d)
    if isstruct(S) && isfield(S,f) && ~isempty(S.(f)), v = S.(f); else, v = d; end
end
