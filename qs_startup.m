function qs_startup()
%QS_STARTUP  Put the whole qs_framework source tree on the MATLAB path.
%
%   Run this ONCE at the start of any session, from anywhere:
%       run('/path/to/qs_framework/qs_startup.m')
%   or, if you are already in the framework root:
%       qs_startup
%
%   Every entry point (qs_run_sweep, qs_run_long, qs_merge, ...) assumes this
%   has been called. The SLURM launchers do it for you.
%
%   Nothing else in the framework ever calls addpath, so there is exactly one
%   place where the path is decided.

    here = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(here,'src')));

    % The test drivers go on the path too, so run_all_tests can be called from
    % any working directory. tests/compat is NOT included here -- see below.
    if exist(fullfile(here,'tests'),'dir')==7
        addpath(fullfile(here,'tests'));
    end

    % Octave only: compatibility shims for MATLAB-only functions used by the
    % model code (griddedInterpolant). Never added under MATLAB.
    if exist('OCTAVE_VERSION','builtin')
        shim = fullfile(here,'tests','compat');
        if exist(shim,'dir')==7, addpath(shim); end
    end

    fprintf('qs_framework: path set from %s\n', here);
end
