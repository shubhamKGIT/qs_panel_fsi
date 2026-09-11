function F = griddedInterpolant(gridvecs, V, method, extrap)
%GRIDDEDINTERPOLANT  Minimal Octave stand-in for MATLAB's griddedInterpolant.
%
%   ONLY used when running the test suite under Octave -- qs_startup adds this
%   folder to the path exclusively in that case. Under MATLAB the built-in is
%   used and this file is never seen.
%
%   Supports the one call shape the framework uses:
%       F = griddedInterpolant({xv, yv}, V, 'linear', 'nearest');
%       vq = F(xq, yq);
%   'nearest' extrapolation is implemented by clamping the query point to the
%   grid and then taking the value at the nearest node, which matches MATLAB's
%   behaviour for query points just outside a rectilinear grid.

    if nargin < 3, method = 'linear'; end
    if nargin < 4, extrap = 'nearest'; end
    assert(iscell(gridvecs) && numel(gridvecs)==2, ...
        'compat griddedInterpolant: only the 2-D {xv,yv} form is supported');

    xv = gridvecs{1}(:);
    yv = gridvecs{2}(:);
    F = @(xq,yq) eval_interp(xv, yv, V, xq, yq, method, extrap);
end

function vq = eval_interp(xv, yv, V, xq, yq, method, extrap)
    xq = xq(:);  yq = yq(:);
    if strcmpi(extrap,'nearest')
        xc = min(max(xq, xv(1)), xv(end));
        yc = min(max(yq, yv(1)), yv(end));
    else
        xc = xq;  yc = yq;
    end
    % interp2 expects meshgrid orientation; V here is (numel(xv) x numel(yv))
    vq = interp2(yv.', xv, V, yc, xc, method);
    bad = ~isfinite(vq);
    if any(bad) && strcmpi(extrap,'nearest')
        for k = find(bad).'
            [~, i] = min(abs(xv - xc(k)));
            [~, j] = min(abs(yv - yc(k)));
            vq(k) = V(i,j);
        end
    end
end
