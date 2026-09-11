function outfile = qs_build_aero(caseId)
%QS_BUILD_AERO  Turn a case's Fluent profile exports into its aero data file.
%
%   outfile = QS_BUILD_AERO(caseId)
%
%   THIS IS THE DROP-IN STEP. Put the two exports in cases/<caseId>/aero/,
%   point case.json's "aero_build" block at them, run this once, and the case
%   is ready to sweep. Nothing else in the framework touches .prof files.
%
%   case.json block:
%     "aero_build": {
%       "surf_file": "aero/<panel surface pressure>.prof",
%       "edge_file": "aero/<BL-edge pressure/mach/soundspeed>.prof",
%       "y_shift_m": 0.0125,
%       "surf_fields": {"pressure": "pressure"},
%       "edge_fields": {"pressure": "pressure", "mach": "mach_number",
%                       "sound_speed": "sound_speed"},
%       "out": "aero/rc19_aero_data.mat"
%     }
%   The *_fields maps are optional: the field names are auto-detected from the
%   export and only need overriding when a name differs from the usual ones.
%
%   Both exports live in the CFD case's own frame; the panel-local frame used
%   everywhere else has  x_local = x_cfd,  y_local = y_cfd - y_shift_m.
%   Coverage of the panel footprint is checked before interpolating, which is
%   what catches a wrong shift before it silently corrupts the field.
%
%   After writing, the implied stagnation temperature is reported. For a
%   heated case it must come out near the case's T0_K -- that single number is
%   the check that the heated RANS, not an older cold-flow one, went in.

    C = qs_config_case_only(caseId);
    assert(isfield(C.case,'aero_build'), ...
        'qs_build_aero: cases/%s/case.json has no "aero_build" block', caseId);
    B = C.case.aero_build;

    surfFile = resolve(C.paths.case_dir, req(B,'surf_file'));
    edgeFile = resolve(C.paths.case_dir, req(B,'edge_file'));
    yShift   = getdef(B,'y_shift_m', 0.0125);
    outfile  = resolve(C.paths.case_dir, getdef(B,'out','aero/rc19_aero_data.mat'));

    % ---- target grid: the ROM's own node grid --------------------------
    S  = load(C.paths.rom_file, 'xy_rom');
    nx = 101;  ny = 51;
    if isfield(C.case,'rom_grid')
        nx = C.case.rom_grid(1);  ny = C.case.rom_grid(2);
    end
    assert(size(S.xy_rom,2) == nx*ny, ...
        'qs_build_aero: ROM has %d nodes but rom_grid says %dx%d=%d', ...
        size(S.xy_rom,2), nx, ny, nx*ny);
    x_grid = reshape(S.xy_rom(1,:)', nx, ny);
    y_grid = reshape(S.xy_rom(2,:)', nx, ny);

    Lp = C.case.panel.L_m;
    Wp = C.case.panel.W_m;

    % ---- surface pressure -----------------------------------------------
    fprintf('qs_build_aero: reading %s\n', surfFile);
    surf = qs_read_fluent_profile(surfFile);
    sf   = map_fields(surf, getdef(B,'surf_fields',struct()), {'pressure'}, 'surf_file');
    xs = surf.(sf.x);  ys = surf.(sf.y) - yShift;
    check_coverage('surf_file', xs, ys, Lp, Wp);
    pl_surface = qs_scatter_interp(xs, ys, surf.(sf.pressure), x_grid, y_grid);
    pl_surface = pl_surface(:);

    % ---- BL-edge conditions ----------------------------------------------
    fprintf('qs_build_aero: reading %s\n', edgeFile);
    edge = qs_read_fluent_profile(edgeFile);
    ef   = map_fields(edge, getdef(B,'edge_fields',struct()), ...
                      {'pressure','mach','sound_speed'}, 'edge_file');
    xe = edge.(ef.x);  ye = edge.(ef.y) - yShift;
    check_coverage('edge_file', xe, ye, Lp, Wp);
    p_edge = reshape(qs_scatter_interp(xe, ye, edge.(ef.pressure),    x_grid, y_grid), [], 1);
    M_edge = reshape(qs_scatter_interp(xe, ye, edge.(ef.mach),        x_grid, y_grid), [], 1);
    a_edge = reshape(qs_scatter_interp(xe, ye, edge.(ef.sound_speed), x_grid, y_grid), [], 1);

    % ---- report -----------------------------------------------------------
    g = 1.4;  Rg = 287.058;
    T0 = (a_edge.^2/(g*Rg)) .* (1 + 0.5*(g-1)*M_edge.^2);
    fprintf('\n  pl_surface : %8.2f .. %8.2f kPa\n', min(pl_surface)/1e3, max(pl_surface)/1e3);
    fprintf('  p_edge     : %8.2f .. %8.2f kPa\n', min(p_edge)/1e3, max(p_edge)/1e3);
    fprintf('  M_edge     : %8.3f .. %8.3f\n', min(M_edge), max(M_edge));
    fprintf('  a_edge     : %8.1f .. %8.1f m/s\n', min(a_edge), max(a_edge));
    fprintf('  implied T0 : %8.1f .. %8.1f K  (mean %.1f)\n', min(T0), max(T0), mean(T0));
    if isfield(C.case,'T0_K') && isfinite(C.case.T0_K)
        rel = abs(mean(T0)-C.case.T0_K)/C.case.T0_K;
        fprintf('  case T0    : %8.1f K   -> %.1f%% difference %s\n', C.case.T0_K, 100*rel, ...
            tern(rel<0.05,'(OK)','<<< MISMATCH: is this the heated RANS?'));
    end
    if isfield(C.case,'pc_nom_Pa')
        fprintf('  mean net load p_l - p_c = %+.3f kPa\n', (mean(pl_surface)-C.case.pc_nom_Pa)/1e3);
    end

    % ---- save ---------------------------------------------------------------
    pinf = C.case.pinf_Pa; %#ok<NASGU>
    readme = sprintf(['RC-19 aero data on the structural ROM grid (node-ordered, matching ' ...
        'rom.x_grid(:) / rom.Modes columns). pl_surface = mean steady-RANS SURFACE ' ...
        'pressure (P_rans) [Pa]; p_edge,M_edge,a_edge = BL-edge static pressure [Pa], ' ...
        'Mach [-], sound speed [m/s]. Built by qs_build_aero from %s and %s ' ...
        '(y_local = y_cfd - %.4f m). Implied T0 = %.1f K.'], ...
        surfFile, edgeFile, yShift, mean(T0)); %#ok<NASGU>
    Lp_out = Lp; %#ok<NASGU>

    d = fileparts(outfile);
    if exist(d,'dir')~=7, mkdir(d); end
    Lp = Lp_out; %#ok<NASGU>
    save(outfile, 'x_grid','y_grid','nx','ny','pl_surface','p_edge','M_edge','a_edge', ...
         'pinf','Lp','readme',qs_matver());
    fprintf('\n  saved %s\n\n', outfile);
end

% ======================================================================
function C = qs_config_case_only(caseId)
%QS_CONFIG_CASE_ONLY  Resolve just the case block, with no study attached.
    root = qs_root();
    caseDir = fullfile(root,'cases',caseId);
    assert(exist(caseDir,'dir')==7, 'qs_build_aero: no case folder %s', caseDir);
    C.case = qs_jsonread(fullfile(caseDir,'case.json'));
    C.case.case_id = caseId;
    C.paths.case_dir = caseDir;
    C.paths.rom_file = resolve(caseDir, C.case.rom_file);
end

function f = map_fields(S, overrides, needed, tag)
%MAP_FIELDS  Work out which struct field holds x, y and each needed quantity.
    have = fieldnames(S);
    f.x = pick(have, overrides, 'x', {'x','x_coordinate','xcoord'}, tag);
    f.y = pick(have, overrides, 'y', {'y','y_coordinate','ycoord'}, tag);
    cand = struct( ...
        'pressure',    {{'pressure','p','static_pressure','absolute_pressure'}}, ...
        'mach',        {{'mach_number','mach','machnumber'}}, ...
        'sound_speed', {{'sound_speed','soundspeed','c','a'}});
    for k = 1:numel(needed)
        nm = needed{k};
        f.(nm) = pick(have, overrides, nm, cand.(nm), tag);
    end
end

function name = pick(have, overrides, key, candidates, tag)
    if isfield(overrides, key) && ~isempty(overrides.(key))
        name = overrides.(key);
        assert(any(strcmp(name, have)), ...
            'qs_build_aero: %s has no field "%s" (available: %s)', tag, name, strjoin(have,', '));
        return
    end
    for k = 1:numel(candidates)
        j = find(strcmpi(candidates{k}, have), 1);
        if ~isempty(j), name = have{j}; return; end
    end
    error(['qs_build_aero: could not find the "%s" column in %s.\n' ...
           '  Available fields: %s\n' ...
           '  Add an explicit mapping in case.json, e.g. "%s_fields": {"%s": "<name>"}'], ...
           key, tag, strjoin(have,', '), strrep(tag,'_file',''), key);
end

function check_coverage(tag, x, y, Lp, W)
    pad = 0.01;
    if min(x) > pad || max(x) < Lp-pad || min(y) > pad || max(y) < W-pad
        warning('qs_build_aero:coverage', ...
            ['%s: point cloud x=[%.4f %.4f] y=[%.4f %.4f] does not comfortably cover\n' ...
             '  the panel footprint x=[0 %.3f] y=[0 %.3f]. Check y_shift_m.'], ...
            tag, min(x), max(x), min(y), max(y), Lp, W);
    end
end

function p = resolve(baseDir, p)
    if isempty(p), return; end
    if p(1)=='/' || p(1)=='~' || (numel(p)>1 && p(2)==':'), return; end
    p = fullfile(baseDir, p);
end

function v = req(S, key)
    assert(isfield(S,key), 'qs_build_aero: aero_build."%s" is required', key);
    v = S.(key);
end

function v = getdef(S, f, d)
    if isstruct(S) && isfield(S,f) && ~isempty(S.(f)), v = S.(f); else, v = d; end
end

function s = tern(c,a,b)
    if c, s = a; else, s = b; end
end
