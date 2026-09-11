function tf = qs_isoctave()
%QS_ISOCTAVE  True when running under GNU Octave rather than MATLAB.
%   Used only to guard MATLAB-only conveniences (matfile, exportgraphics,
%   sgtitle). The physics path is identical in both.
    tf = logical(exist('OCTAVE_VERSION','builtin'));
end
