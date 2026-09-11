function y0 = qs_initial_state(M, C, ic_tag)
%QS_INITIAL_STATE  Build the initial modal state for one point.
%
%   y0 = QS_INITIAL_STATE(M, C, ic_tag)   returns [q0; qdot0]  (2*nm x 1)
%
%   Supported tags:
%     'flat'              q0 = +seed on every mode, qdot0 = 0.
%                         This is what every existing sweep used: a flat panel
%                         with a uniform seed. Above dT_cr the panel is buckled
%                         and two wells exist, so a map built this way shows
%                         "the attractor reached from a flat panel", not "the
%                         attractor". That is reproducible and defensible, but
%                         it has to be stated on any figure made from it.
%     'flat_neg'          q0 = -seed. Same run from the mirrored start; pairing
%                         the two is the cheapest possible basin probe.
%     'from:<point_id>'   continue from the final state saved by an earlier
%                         point (requires that point to have been run with
%                         store.save_final_state = true).
%
%   Seed magnitude is C.ic.seed (default 1e-4), matching the original drivers.

    nm   = M.nm;
    seed = 1e-4;
    if isfield(C,'ic') && isfield(C.ic,'seed') && ~isempty(C.ic.seed)
        seed = C.ic.seed;
    end

    y0 = zeros(2*nm,1);

    if strncmp(ic_tag, 'from:', 5)
        parent = ic_tag(6:end);
        f = fullfile(C.paths.hist_dir, [parent '_state.mat']);
        assert(exist(f,'file')==2, ...
            ['qs_initial_state: ic_tag "%s" needs the final state of point %s,\n' ...
             '  expected at %s\n' ...
             '  Run that point first with store.save_final_state = true.'], ...
            ic_tag, parent, f);
        S = load(f);
        y0 = S.yend(:);
        assert(numel(y0)==2*nm, 'qs_initial_state: saved state has %d entries, expected %d', numel(y0), 2*nm);
        return
    end

    switch lower(ic_tag)
        case 'flat'
            y0(1:nm) =  seed;
        case 'flat_neg'
            y0(1:nm) = -seed;
        case 'rest'
            % exact zero: only useful for linear checks, never for a buckled panel
        otherwise
            error('qs_initial_state: unknown ic_tag "%s"', ic_tag);
    end
end
