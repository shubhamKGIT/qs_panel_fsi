function G = qs_rasterize(D, field)
%QS_RASTERIZE  Put a merged point table back onto a rectangular grid for plotting.
%
%   G = QS_RASTERIZE(D)          rasterizes the standard set of fields
%   G = QS_RASTERIZE(D, field)   rasterizes one named field
%
%   The point list is deliberately NOT a grid: refinement inserts midpoints
%   only where the boundary is, so the set of (p_c, dT) values is ragged. This
%   function builds the lattice of all unique values that DO occur and drops
%   each point into its cell. Cells with no point stay NaN -- they are holes in
%   a ragged lattice, not failures, and the plotting code shows them as gaps.
%
%   WHERE SEVERAL RUNS SHARE AN OPERATING POINT (a 3 s coarse run and a 10 s
%   re-run, say) the LONGEST run wins. That is the whole reason both are kept:
%   the short answer stays on the record in points.csv, while the map shows the
%   answer with the most evidence behind it.
%
%   OUTPUT (struct G)
%     .pc_kPa, .dT_K     axis vectors (sorted unique values)
%     .<field>           npc x ndT array for each rasterized field
%     .t_end             which run length produced each cell
%     .filled            logical, true where a point exists
%     .n_superseded      how many shorter runs were hidden by a longer one

    if nargin < 2 || isempty(field)
        fields = {'label','Amp_wh','Freq','pfrac','MeanPeak','StdPeak', ...
                  'MeanRMSE','StdRMSE','trans_time','trans_from','trans_to', ...
                  'level','amp_last'};
        logfields = {'trans_flag','nonstationary','near_threshold','divergent','ok'};
    else
        fields = {field};  logfields = {};
    end

    ok = D.ok;
    pc_kPa = unique(round(D.pc_Pa(ok))/1e3);
    dT_K   = unique(round(D.dT_K(ok)*1e3)/1e3);
    npc = numel(pc_kPa);  ndT = numel(dT_K);

    G.pc_kPa = pc_kPa;   G.dT_K = dT_K;
    G.filled = false(npc,ndT);
    G.t_end  = nan(npc,ndT);
    for f = fields,    G.(f{1}) = nan(npc,ndT);   end
    for f = logfields, G.(f{1}) = false(npc,ndT); end

    % longest run first, so a later (shorter) run never overwrites it
    [~, ord] = sort(D.t_end(:), 'descend');
    nsup = 0;

    for kk = 1:numel(ord)
        k = ord(kk);
        if ~D.ok(k), continue; end
        i = nearest_index(pc_kPa, D.pc_Pa(k)/1e3);
        j = nearest_index(dT_K,   D.dT_K(k));
        if isempty(i) || isempty(j), continue; end
        if G.filled(i,j)
            nsup = nsup + 1;
            continue                      % already taken by a longer run
        end
        G.filled(i,j) = true;
        G.t_end(i,j)  = D.t_end(k);
        for f = fields,    G.(f{1})(i,j) = D.(f{1})(k);          end
        for f = logfields, G.(f{1})(i,j) = logical(D.(f{1})(k)); end
    end

    G.n_superseded = nsup;
end

% ----------------------------------------------------------------------
function i = nearest_index(v, x)
    [d, i] = min(abs(v - x));
    if d > 1e-6, i = []; end
end
