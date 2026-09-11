function qs_jsonwrite(fname, S)
%QS_JSONWRITE  Write a struct to a JSON file (pretty-printed where supported).
%   Used for the results manifest, so a run's exact settings can be read
%   without MATLAB.

    try
        txt = jsonencode(S, 'PrettyPrint', true);   % R2021a+
    catch
        txt = jsonencode(S);                        % older: single line
    end
    d = fileparts(fname);
    if ~isempty(d) && exist(d,'dir')~=7, mkdir(d); end
    fid = fopen(fname,'w');
    assert(fid>0, 'qs_jsonwrite: cannot open %s for writing', fname);
    fprintf(fid, '%s\n', txt);
    fclose(fid);
end
