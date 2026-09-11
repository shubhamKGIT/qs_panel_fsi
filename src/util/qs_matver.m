function v = qs_matver()
%QS_MATVER  The .mat format flag every save in the framework uses.
%
%   MATLAB  -> '-v7.3': HDF5-backed. Needed for anything over 2 GB and for
%   PARTIAL loading, which is how the regression test reads a few variables
%   out of a 145 MB map instead of pulling the whole thing into memory.
%
%   Octave  -> '-v7': Octave cannot write v7.3. Octave is only used for
%   syntax and plumbing checks here, never for production runs, so the 2 GB
%   limit never bites.
    if exist('OCTAVE_VERSION','builtin')
        v = '-v7';
    else
        v = '-v7.3';
    end
end
