function P = qs_make_points(C)
%QS_MAKE_POINTS  Build the level-0 point list from the study's grid block.
%
%   P = QS_MAKE_POINTS(C)  with C from qs_config.
%
%   Reads C.study.grid:
%       grid.pc_kPa   axis spec (see qs_axis), in kPa
%       grid.dT_K     axis spec, in K
%   and C.duration for the base run length.
%
%   ORDERING is deliberate: p_c varies fastest, then dT. That is the same
%   linear order the original drivers got from ind2sub([npc ndT], idx), so
%   slicing the list k:n:end assigns the same points to the same task number
%   as the old sweep did. It makes the regression test an exact comparison
%   rather than an approximate one.

    assert(isfield(C.study,'grid'), 'qs_make_points: study has no "grid" block');
    G = C.study.grid;

    pc_kPa = qs_axis(G.pc_kPa, 'grid.pc_kPa');
    dT_K   = qs_axis(G.dT_K,   'grid.dT_K');

    npc = numel(pc_kPa);  ndT = numel(dT_K);
    [PC, DT] = ndgrid(pc_kPa, dT_K);      % pc fastest -> column-major
    pc_Pa = PC(:)*1e3;
    dT    = DT(:);

    P = qs_pointlist('make', pc_Pa, dT, struct( ...
            't_end',   C.duration.base_t_end, ...
            't_trans', C.duration.base_t_transient, ...
            'level',   0, ...
            'ic_tag',  ic_tag_from(C), ...
            'origin',  'coarse'));

    % points the study explicitly wants full histories for
    if isfield(C.store,'history_always_ids') && ~isempty(C.store.history_always_ids)
        want = C.store.history_always_ids;
        if ischar(want), want = {want}; end
        P.save_hist = P.save_hist | qs_pointlist('has', P, want(:));
    end

    fprintf('qs_make_points: %d x %d = %d points  (p_c %.2f..%.2f kPa, dT %.2f..%.2f K)\n', ...
        npc, ndT, numel(P.id), min(pc_kPa), max(pc_kPa), min(dT_K), max(dT_K));
end

% ----------------------------------------------------------------------
function tag = ic_tag_from(C)
    tag = 'flat';
    if isfield(C,'ic') && isfield(C.ic,'mode') && ~isempty(C.ic.mode)
        tag = C.ic.mode;
    end
end
