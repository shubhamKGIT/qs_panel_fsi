function S = qs_jsonread(fname)
%QS_JSONREAD  Read a JSON file into a struct, with a useful error message.
%
%   Comment lines are NOT part of JSON. To keep notes inside a config file,
%   use a key whose name starts with an underscore (e.g. "_note"): those are
%   read like any other field and simply ignored by the framework.

    assert(exist(fname,'file')==2, 'qs_jsonread: file not found:\n  %s', fname);
    if exist('jsondecode','builtin') ~= 5 && exist('jsondecode','file') == 0
        error(['qs_jsonread: this MATLAB/Octave has no jsondecode. ' ...
               'MATLAB R2016b+ or Octave 7+ is required.']);
    end
    txt = fileread(fname);
    try
        S = jsondecode(txt);
    catch ME
        error('qs_jsonread: %s is not valid JSON.\n  %s', fname, ME.message);
    end
    assert(isstruct(S) && isscalar(S), ...
        'qs_jsonread: %s must contain a single JSON object at the top level.', fname);
end
