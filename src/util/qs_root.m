function r = qs_root()
%QS_ROOT  Absolute path of the qs_framework root directory.
%   Derived from this file's own location (src/util/qs_root.m), so it is
%   correct no matter what the current working directory is.
    d = fileparts(mfilename('fullpath'));   % .../src/util
    r = fileparts(fileparts(d));            % .../
end
