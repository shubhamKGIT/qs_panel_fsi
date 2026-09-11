function qs_savefig(fh, basename, data)
%QS_SAVEFIG  Save a figure as .fig and .png, plus the data behind it.
%   qs_savefig(fh, basename, data)
%   The .mat companion is the point: a figure that cannot be regenerated or
%   re-plotted from its own numbers is a dead end six months later.
    try, savefig(fh, [basename '.fig']); catch, end
    try
        exportgraphics(fh, [basename '.png'], 'Resolution', 150);
    catch
        try, print(fh, [basename '.png'], '-dpng', '-r150'); catch, end
    end
    if nargin > 2 && ~isempty(data)
        try, save([basename '_data.mat'], '-struct', 'data', qs_matver()); catch, end
    end
end
