function make(varargin)
% MAKE  Build the pipeline from inside MATLAB - the Makefile equivalent.
%
% Purpose
%   One command that runs the stages in the right order, on machines without
%   `make` (Windows), and without having to remember which stages shard. Every
%   stage still performs its own staleness check, so re-running is cheap and
%   only out-of-date work is redone.
%
% Usage (command syntax works, so the quotes are optional)
%   make                       % everything out of date, all three days
%   make all                   % same
%   make data                  % gps, slim, samples, fields
%   make figures               % macro, micro
%   make slim                  % one stage, all three days
%   make slim Days 18          % one stage, one day
%   make all Days 18 Workers 6 % one day, slim/full across 6 processes
%   make status                % what is stale, and why; builds nothing
%   make config                % show the resolved paths and settings
%   make all DryRun true       % plan only, write nothing
%   make slim Force true       % rebuild regardless of timestamps
%
% Inputs
%   target   'all' (default) | 'data' | 'figures' | a single stage name
%            ('gps', 'slim', 'full', 'samples', 'fields', 'macro', 'micro',
%            'av') | 'status' | 'config'
%   Name/value:
%     'Log'      write a timestamped log under <results>/.mvt/logs
%                (default true; Log false disables it)
%     'KeepGoing' attempt every stage and report failures at the end,
%                instead of stopping at the first one (default false)
%     'Workers'  MATLAB processes for the stages that shard (default 1).
%                Only 'slim' and 'full' shard; everything else ignores it.
%                Size by memory, not cores: each worker peaks at several GB.
%     'Days'     days to build (default [16 17 18])
%     ...        any other option goes to mvt.options: Force, Clean, DryRun,
%                Verbose, SettleSeconds
%
% Notes
%   Command syntax passes every argument as text, so `make slim Days 18` gives
%   '18' rather than 18; numeric-looking values are converted here.
%
%   'av' is cross-day and runs once, after every requested day, because it
%   reads all three days' samples.
%
% Dependencies
%   mvt.options, mvt.build, mvt.runShards, mvt.status, mvt.paths
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

[target, workers, keepGoing, doLog, opts, configArgs] = parseArguments(varargin);

switch target
    case 'status'
        mvt.status(opts);
        return
    case 'config'
        configure(configArgs, opts, workers);
        return
end

stages = expandTarget(target);
days = opts.Days;   % already resolved in parseArguments: explicit > config > all
crossDay = ismember('av', stages);
perDay = stages(~strcmp(stages, 'av'));

logCleanup = startLog(target, doLog); %#ok<NASGU> closes the diary on any exit
announceEnvironment();
fprintf('[make] target ''%s'': %s\n', target, strjoin(stages, ' '));
fprintf('[make] days %s, %d worker(s) for sharded stages\n', ...
    mat2str(days), workers);
started = tic;

failures = {};
for day = days
    for iStage = 1:numel(perDay)
        stage = perDay{iStage};
        if workers > 1 && ismember(stage, {'slim', 'full'})
            thunk = @() mvt.runShards(stage, day, workers, opts);
        else
            thunk = @() mvt.build(stage, day, opts);
        end
        failures = runStage(thunk, stage, day, opts, keepGoing, failures);
    end
end

if crossDay
    failures = runStage(@() mvt.build('av', [], opts), 'av', [], opts, ...
        keepGoing, failures);
end

fprintf('[make] target ''%s'': days %s finished in %s\n', ...
    target, mat2str(days), humanTime(toc(started)));
if isempty(failures)
    fprintf('[make] all stages completed\n');
else
    fprintf(2, '[make] %d stage(s) FAILED:\n', numel(failures));
    for k = 1:numel(failures)
        fprintf(2, '         %s\n', failures{k});
    end
    if ~opts.DryRun
        error('mvt:make:stageFailed', ...
            '%d stage(s) failed; see the list above.', numel(failures));
    end
end
end

% ---------------------------------------------------------------------------
function cleanup = startLog(target, doLog)
% Tee the whole run to a timestamped file under <results>/.mvt/logs.
%
% A build takes hours and prints thousands of lines; when something goes wrong
% the command window has usually scrolled past it, or the session has been
% closed. The diary keeps the whole thing, including the stage output and the
% error report, next to the data it was building.
%
% Returns an onCleanup object: the diary closes when make returns, whether it
% finished, errored or was interrupted.
cleanup = [];
if ~doLog
    return
end
try
    p = mvt.paths();
    logDir = fullfile(p.stateDir, 'logs');
    mvt.ensureDir(logDir);
    stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    logFile = fullfile(logDir, sprintf('make-%s-%s.log', target, stamp));

    % Preserve a diary the caller already had running.
    priorState = get(0, 'Diary');
    priorFile = get(0, 'DiaryFile');
    diary(logFile);
    cleanup = onCleanup(@() closeLog(logFile, priorState, priorFile));

    fprintf('[make] logging to %s\n', logFile);
    fprintf('[make] %s | MATLAB %s | %s | data version %s\n', ...
        char(datetime('now')), version('-release'), computer(), mvt.dataVersion());
catch err
    fprintf(2, '[make] could not open a log file (%s); continuing without one\n', ...
        err.identifier);
    cleanup = [];
end
end

% ---------------------------------------------------------------------------
function closeLog(logFile, priorState, priorFile)
fprintf('[make] log written to %s\n', logFile);
diary off
if strcmpi(priorState, 'on')
    diary(priorFile);
end
end

% ---------------------------------------------------------------------------
function announceEnvironment()
% MVT_* variables silently change what a build does - most sharply MVT_DAYS,
% which narrows the run to a subset without any other outward sign. setenv
% persists for the whole MATLAB session, so one left over from an earlier
% experiment quietly shrinks every later build. Say so up front.
names = {'MVT_DAYS', 'MVT_DATA_DIR', 'MVT_RESULTS_DIR', 'MVT_FORCE', ...
    'MVT_CLEAN', 'MVT_DRYRUN', 'MVT_VERBOSE', 'MVT_SHARD', 'MVT_SETTLE_SECONDS'};
for k = 1:numel(names)
    value = strtrim(getenv(names{k}));
    if isempty(value)
        continue
    end
    if strcmp(names{k}, 'MVT_DAYS')
        fprintf('[make] environment: %s = %s (IGNORED; make builds every day)\n', ...
            names{k}, value);
    else
        fprintf('[make] environment: %s = %s\n', names{k}, value);
    end
end
end

% ---------------------------------------------------------------------------
function failures = runStage(thunk, stage, day, opts, keepGoing, failures)
% Run one stage, recording rather than raising when asked to carry on.
%
% Under DryRun a later stage usually cannot even look at its inputs, because
% the earlier stage that would have produced them wrote nothing, so a plan
% always continues. A real run stops at the first failure unless KeepGoing is
% set, in which case every stage is attempted and the failures are reported
% together at the end - one bad day should not throw away the other two.
try
    thunk();
catch err
    reason = err.identifier;
    if isempty(reason)          % MATLAB errors need not carry an identifier
        reason = 'error';
    end
    label = sprintf('%s%s (%s)', stage, dayLabel(day), reason);
    if opts.DryRun
        fprintf('[make] (dry-run) %s%s: cannot plan yet (%s)\n', stage, ...
            dayLabel(day), reason);
        return
    end
    if ~keepGoing
        fprintf(2, '[make] %s FAILED. Re-run with KeepGoing true to attempt the rest.\n', label);
        % Into the diary as well: in -batch the error text goes to stderr and
        % would otherwise be missing from the log that exists to explain it.
        disp(getReport(err, 'extended', 'hyperlinks', 'off'));
        rethrow(err);
    end
    fprintf(2, '[make] %s FAILED, continuing (KeepGoing)\n', label);
    disp(getReport(err, 'extended', 'hyperlinks', 'off'));
    failures{end+1} = label; %#ok<AGROW>
end
end

% ---------------------------------------------------------------------------
function label = dayLabel(day)
if isempty(day)
    label = '';
else
    label = sprintf(' 2022-11-%d', day);
end
end

% ---------------------------------------------------------------------------
function stages = expandTarget(target)
switch target
    case 'all'
        stages = {'gps', 'slim', 'samples', 'fields', 'macro', 'micro', 'av'};
    case 'data'
        stages = {'gps', 'slim', 'samples', 'fields'};
    case 'figures'
        stages = {'macro', 'micro'};
    case {'gps', 'slim', 'full', 'samples', 'fields', 'macro', 'micro', 'av'}
        stages = {target};
    otherwise
        error('mvt:make:unknownTarget', ...
            ['Unknown target ''%s''. Expected all, data, figures, status, ', ...
             'config, or a stage: gps slim full samples fields macro micro av.'], ...
            target);
end
end

% ---------------------------------------------------------------------------
function [target, workers, keepGoing, doLog, opts, configArgs] = parseArguments(args)
target = 'all';
if ~isempty(args) && (ischar(args{1}) || isstring(args{1})) ...
        && ~isOptionName(args{1})
    target = lower(char(args{1}));
    args(1) = [];
end

% Command syntax delivers everything as text; make the values usable.
for k = 2:2:numel(args)
    value = args{k};
    if ischar(value) || isstring(value)
        text = strtrim(char(value));
        numeric = str2num(text); %#ok<ST2NM> accepts '18' and '[16 17]'
        if ~isempty(numeric)
            args{k} = numeric;
        elseif any(strcmpi(text, {'true', 'false'}))
            args{k} = strcmpi(text, 'true');
        end
    end
end

workers = 1;
keepGoing = false;
doLog = true;
keep = true(1, numel(args));
for k = 1:2:numel(args) - 1
    switch lower(char(args{k}))
        case 'workers'
            workers = args{k + 1};
            keep(k:k + 1) = false;
        case 'keepgoing'
            keepGoing = logical(args{k + 1});
            keep(k:k + 1) = false;
        case 'log'
            doLog = logical(args{k + 1});
            keep(k:k + 1) = false;
    end
end
configArgs = args(keep);
% `make config reset` is a bare word, not a name/value pair; mvt.options would
% reject it, so keep it out of the options struct.
optionArgs = configArgs;
if numel(optionArgs) == 1 && (ischar(optionArgs{1}) || isstring(optionArgs{1})) ...
        && any(strcmpi(char(optionArgs{1}), {'reset', 'clear'}))
    optionArgs = {};
end
opts = mvt.options(optionArgs{:});

% `make all` builds every day. mvt.options honors MVT_DAYS - which the Unix
% Makefile needs, and which is why it stays there - but an ambient variable
% must not silently decide what a build covers. Precedence here is explicit
% argument, then a value set with `make config`, then all three days.
gaveDays = false;
for k = 1:2:numel(optionArgs) - 1
    if strcmpi(char(optionArgs{k}), 'Days')
        gaveDays = true;
    end
end
if ~gaveDays
    opts.Days = sessionDays();
end

if ~isscalar(workers) || workers < 1 || workers ~= fix(workers)
    error('mvt:make:badWorkers', 'Workers must be a positive integer.');
end
end

% ---------------------------------------------------------------------------
function days = sessionDays(newDays)
% Days for a build with no explicit Days argument: whatever `make config Days`
% last set in this MATLAB session, otherwise all three. Deliberately session
% scoped and in memory - a setting that outlived the session would recreate the
% invisible state that MVT_DAYS caused.
persistent configured
if nargin > 0
    configured = newDays;
end
if isempty(configured)
    days = [16 17 18];
else
    days = configured;
end
end

% ---------------------------------------------------------------------------
function tf = isOptionName(value)
tf = any(strcmpi(char(value), {'Force', 'Clean', 'DryRun', 'Verbose', ...
    'Shard', 'Days', 'SettleSeconds', 'Workers', 'KeepGoing', 'Log', 'UseParfor'}));
end

% ---------------------------------------------------------------------------
function configure(configArgs, opts, workers)
% `make config` prints the resolved settings; `make config Days 18` sets the
% days used by later builds in this session, and `make config reset` clears it.
if ~isempty(configArgs)
    if numel(configArgs) == 1 && any(strcmpi(char(configArgs{1}), {'reset', 'clear'}))
        sessionDays([]);
        fprintf('  config reset: builds cover all three days again\n');
    else
        for k = 1:2:numel(configArgs) - 1
            if strcmpi(char(configArgs{k}), 'Days')
                sessionDays(configArgs{k + 1});
                fprintf('  config set: days = %s (this MATLAB session)\n', ...
                    mat2str(configArgs{k + 1}));
            else
                fprintf(2, '  config: %s is not a persistent setting; pass it to the build\n', ...
                    char(configArgs{k}));
            end
        end
    end
    opts.Days = sessionDays();
end
p = mvt.paths();
fprintf('  repo      = %s\n', p.repoRoot);
fprintf('  data      = %s\n', p.dataDir);
fprintf('  results   = %s\n', p.resultsDir);
fprintf('  days      = %s\n', mat2str(sessionDays()));
fprintf('  workers   = %d\n', workers);
fprintf('  force     = %d\n', opts.Force);
fprintf('  dry run   = %d\n', opts.DryRun);
fprintf('  data ver. = %s\n', mvt.dataVersion());
stray = strtrim(getenv('MVT_DAYS'));
if ~isempty(stray)
    fprintf(2, ['  note: MVT_DAYS = %s is set but IGNORED by make; builds cover\n' ...
        '        the days above. Use `make all Days ...` or `make config Days ...`.\n'], stray);
end
end

% ---------------------------------------------------------------------------
function text = humanTime(seconds)
if seconds < 60
    text = sprintf('%ds', round(seconds));
elseif seconds < 3600
    text = sprintf('%dm%02ds', floor(seconds/60), round(mod(seconds, 60)));
else
    text = sprintf('%dh%02dm', floor(seconds/3600), floor(mod(seconds, 3600)/60));
end
end
