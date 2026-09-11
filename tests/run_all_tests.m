function ok = run_all_tests(which_tests)
%RUN_ALL_TESTS  Sanity-check the framework and write a report you can send back.
%
%   run_all_tests               run every test
%   run_all_tests([1 2 5])      run only those
%
%   Tests, in increasing cost:
%     t01  config, grids, IDs, point list      no physics, < 1 s
%     t02  model build + preflight             loads the case data, ~10 s
%     t03  single-point regression vs the old map   ~3 min per point
%     t04  sweep -> merge -> refine round trip  ~2 min (short runs)
%     t05  windows / classify / transition      synthetic signals, ~5 s
%
%   START WITH  run_all_tests([1 2 5])  -- those three exercise every piece of
%   new plumbing without waiting for an integration. If they pass, run 3 and 4.
%
%   A report is written to tests/last_test_report.txt. Send that file back if
%   anything fails; it carries the failing values, not just the word FAIL.

    if nargin < 1 || isempty(which_tests), which_tests = 1:5; end

    here = fileparts(mfilename('fullpath'));
    run(fullfile(fileparts(here),'qs_startup.m'));

    tests = { @t01_config_and_points, 't01 config / grids / IDs / point list'
              @t02_model_and_preflight,'t02 model build + preflight'
              @t03_point_regression,   't03 single-point regression vs old map'
              @t04_sweep_roundtrip,    't04 sweep -> merge -> refine round trip'
              @t05_windows_classify,   't05 windows / classify / transition' };

    lines = {};
    lines{end+1} = sprintf('qs_framework test report   %s', datestr(now,'yyyy-mm-dd HH:MM:SS')); %#ok<TNOW1,DATST>
    lines{end+1} = sprintf('root    : %s', qs_root());
    lines{end+1} = sprintf('engine  : %s', version());
    lines{end+1} = repmat('=',1,72);

    nfail = 0;  nskip = 0;
    for k = which_tests(:).'
        name = tests{k,2};
        fprintf('\n>>> %s\n', name);
        t0 = tic;
        try
            R = tests{k,1}();
        catch ME
            R = struct('pass',false,'skip',false,'msg',{{sprintf('EXCEPTION %s: %s', ME.identifier, ME.message)}});
            for s = 1:min(5,numel(ME.stack))
                R.msg{end+1} = sprintf('    at %s line %d', ME.stack(s).name, ME.stack(s).line);
            end
        end
        el = toc(t0);

        if isfield(R,'skip') && R.skip
            status = 'SKIP';  nskip = nskip + 1;
        elseif R.pass
            status = 'PASS';
        else
            status = 'FAIL';  nfail = nfail + 1;
        end
        fprintf('<<< %s  (%.1f s)\n', status, el);

        lines{end+1} = sprintf('[%s] %-45s %6.1f s', status, name, el); %#ok<AGROW>
        if isfield(R,'msg')
            for m = 1:numel(R.msg)
                lines{end+1} = sprintf('       %s', R.msg{m}); %#ok<AGROW>
            end
        end
    end

    lines{end+1} = repmat('=',1,72);
    lines{end+1} = sprintf('%d test(s) run, %d failed, %d skipped', numel(which_tests), nfail, nskip);

    rep = fullfile(here,'last_test_report.txt');
    fid = fopen(rep,'w');
    for k = 1:numel(lines), fprintf(fid,'%s\n', lines{k}); end
    fclose(fid);

    fprintf('\n%s\n', repmat('=',1,72));
    for k = 1:numel(lines), fprintf('%s\n', lines{k}); end
    fprintf('\nreport written to %s\n', rep);

    ok = (nfail == 0);
end
