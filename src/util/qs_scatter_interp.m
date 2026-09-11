function vq = qs_scatter_interp(x, y, v, xq, yq)
%QS_SCATTER_INTERP  Linear scattered interpolation with nearest-neighbour fill.
%
%   Matches the behaviour rc19_build_aero_data relied on:
%       scatteredInterpolant(x, y, v, 'linear', 'nearest')
%   i.e. linear inside the convex hull of the input cloud, nearest-neighbour
%   for the handful of target nodes just outside it (the exact panel edges,
%   since Fluent exports cell-centre data that does not quite reach them).
%
%   Falls back to griddata + an explicit nearest fill where
%   scatteredInterpolant is unavailable, so the aero build also runs under
%   Octave for testing.

    x = x(:);  y = y(:);  v = v(:);
    if exist('scatteredInterpolant','class') == 8 || exist('scatteredInterpolant','file') == 2
        F  = scatteredInterpolant(x, y, v, 'linear', 'nearest');
        vq = F(xq(:), yq(:));
    else
        vq = griddata(x, y, v, xq(:), yq(:), 'linear');
        bad = ~isfinite(vq);
        if any(bad)
            qx = xq(:);  qy = yq(:);
            idx = find(bad);
            for k = 1:numel(idx)
                [~, j] = min((x-qx(idx(k))).^2 + (y-qy(idx(k))).^2);
                vq(idx(k)) = v(j);
            end
        end
    end
    vq = reshape(vq, size(xq));
end
