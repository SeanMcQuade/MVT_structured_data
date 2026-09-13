function failed = runShards(stage, day, nWorkers, varargin)
% MVT.RUNSHARDS  Run one stage across N background MATLAB processes.
%
% Purpose
%   The `make -j` equivalent, callable from inside MATLAB, for machines with no
%   make (Windows). There is no in-process parallelism anywhere in this
%   pipeline: no stage uses parfor, `opts.UseParfor` is a no-op, and a plain
%   mvt.build call is strictly sequential. Concurrency comes from running
%   several MATLAB processes over disjoint shards of the same day, which is
%   exactly what the Makefile does.
%
% Inputs
%   stage      'slim' | 'full' | 'gps' | ... (any mvt.build stage)
%   day        16, 17, or 18
%   nWorkers   number of MATLAB processes to launch (1 runs in-process)
%   Name/value:
%     'PollSeconds'  how often to check for completion (default 10)
%     ...            any other option is forwarded to mvt.build/mvt.options
%
% Outputs
%   failed  number of shards that reported an error (0 on success). Errors are
%           also raised unless an output is requested.
%
% Example
%   mvt.runShards('slim', 18, 8)        % eight concurrent MATLAB processes
%
% Notes
%   Only stages that shard by segment gain from this: 'slim' and 'full' split
%   their 24 segments via mvt.shardIndices. 'gps', 'samples', 'fields', 'macro'
%   and 'micro' ignore Shard and would each do the whole job N times, so this
%   refuses to launch more than one worker for them.
%
%   Size nWorkers by memory, not cores: each process holds a decoded segment,
%   several GB at peak.
%
%   Each worker logs to <results>/.mvt/logs/<stage>-<day>-shard<k>of<N>.log and
%   writes a status file when it exits, which is how completion and failure are
%   detected across processes.
%
% Dependencies
%   mvt.paths, mvt.options, mvt.ensureDir, mvt.build (in the workers)
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

shardedStages = {'slim', 'full'};

parser = inputParser;
parser.KeepUnmatched = true;
parser.addParameter('PollSeconds', 10);
parser.parse(varargin{:});
pollSeconds = parser.Results.PollSeconds;
passthrough = parser.Unmatched;

stage = lower(char(stage));
mvt.assertDay(day);
if ~isscalar(nWorkers) || nWorkers < 1 || nWorkers ~= fix(nWorkers)
    error('mvt:runShards:badWorkers', 'nWorkers must be a positive integer.');
end
if nWorkers > 1 && ~ismember(stage, shardedStages)
    error('mvt:runShards:notSharded', ...
        ['Stage ''%s'' does not shard; %d workers would each repeat the whole ', ...
         'job. Only %s split their segments.'], ...
        stage, nWorkers, strjoin(shardedStages, ' and '));
end

if nWorkers == 1
    mvt.build(stage, day, passthrough);
    failed = 0;
    return
end

p = mvt.paths();
logDir = fullfile(p.stateDir, 'logs');
mvt.ensureDir(logDir);
stamp = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>

matlabExe = fullfile(matlabroot, 'bin', 'matlab');
markers = cell(1, nWorkers);
logs = cell(1, nWorkers);

fprintf('[mvt] %s 2022-11-%d: launching %d MATLAB workers\n', stage, day, nWorkers);
for k = 1:nWorkers
    tag = sprintf('%s-%d-shard%dof%d', stage, day, k, nWorkers);
    logs{k} = fullfile(logDir, [tag '.log']);
    markers{k} = fullfile(logDir, [tag '.' stamp '.status']);
    inner = sprintf(['cd(''%s''); try, mvt.build(''%s'', %d, ''Shard'', [%d %d]); ' ...
        's = 0; catch err, disp(getReport(err)); s = 1; end; ' ...
        'fid = fopen(''%s'', ''w''); fprintf(fid, ''%%d'', s); fclose(fid); exit(s);'], ...
        p.scriptsDir, stage, day, k, nWorkers, markers{k});
    % -logfile rather than shell redirection: on Windows, wrapping this in
    % `cmd /c "... > ""log"" 2>&1"` nests quotes three deep and is fragile.
    % MATLAB writes the log itself, so the spawn line stays simple on both
    % platforms. `inner` contains only single-quoted MATLAB strings, so the
    % -batch argument needs no escaping.
    cmd = sprintf('"%s" -batch "%s" -logfile "%s"', matlabExe, inner, logs{k});
    if ispc
        spawn = sprintf('start "mvt %s" /B %s', tag, cmd);
    else
        spawn = sprintf('%s &', cmd);
    end
    [status, msg] = system(spawn);
    if status ~= 0
        error('mvt:runShards:spawnFailed', ...
            'Could not start worker %d of %d: %s', k, nWorkers, strtrim(msg));
    end
    fprintf('[mvt]   shard %d/%d -> %s\n', k, nWorkers, logs{k});
end

fprintf('[mvt] waiting for %d workers (polling every %ds)\n', nWorkers, pollSeconds);
started = tic;
report = mvt.progress(nWorkers, sprintf('%s 2022-11-%d workers', stage, day), ...
    'UpdateSeconds', pollSeconds);
done = false(1, nWorkers);
while ~all(done)
    pause(pollSeconds);
    for k = 1:nWorkers
        done(k) = done(k) || isfile(markers{k});
    end
    report(sum(done));
end
report();

failed = 0;
for k = 1:nWorkers
    status = str2double(strtrim(fileread(markers{k})));
    delete(markers{k});
    if ~isequal(status, 0)
        failed = failed + 1;
        fprintf(2, '[mvt] shard %d/%d FAILED; see %s\n', k, nWorkers, logs{k});
    end
end
fprintf('[mvt] %s 2022-11-%d: %d/%d shards ok in %.0f s\n', ...
    stage, day, nWorkers - failed, nWorkers, toc(started));

if failed > 0 && nargout == 0
    error('mvt:runShards:shardFailed', ...
        '%d of %d shards failed; see the logs in %s', failed, nWorkers, logDir);
end
end
