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

[target, workers, keepGoing, opts] = parseArguments(varargin);

switch target
    case 'status'
        mvt.status(opts);
        return
    case 'config'
        showConfig(opts, workers);
        return
end

stages = expandTarget(target);
days = opts.Days;
crossDay = ismember('av', stages);
perDay = stages(~strcmp(stages, 'av'));

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
function announceEnvironment()
% MVT_* variables silently change what a build does - most sharply MVT_DAYS,
% which narrows the run to a subset without any other outward sign. setenv
% persists for the whole MATLAB session, so one left over from an earlier
% experiment quietly shrinks every later build. Say so up front.
names = {'MVT_DAYS', 'MVT_DATA_DIR', 'MVT_RESULTS_DIR', 'MVT_FORCE', ...
    'MVT_CLEAN', 'MVT_DRYRUN', 'MVT_VERBOSE', 'MVT_SHARD', 'MVT_SETTLE_SECONDS'};
for k = 1:numel(names)
    value = strtrim(getenv(names{k}));
    if ~isempty(value)
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
        rethrow(err);
    end
    fprintf(2, '[make] %s FAILED, continuing (KeepGoing)\n', label);
    fprintf(2, '       %s\n', err.message);
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
function [target, workers, keepGoing, opts] = parseArguments(args)
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
keep = true(1, numel(args));
for k = 1:2:numel(args) - 1
    switch lower(char(args{k}))
        case 'workers'
            workers = args{k + 1};
            keep(k:k + 1) = false;
        case 'keepgoing'
            keepGoing = logical(args{k + 1});
            keep(k:k + 1) = false;
    end
end
opts = mvt.options(args{keep});

if ~isscalar(workers) || workers < 1 || workers ~= fix(workers)
    error('mvt:make:badWorkers', 'Workers must be a positive integer.');
end
end

% ---------------------------------------------------------------------------
function tf = isOptionName(value)
tf = any(strcmpi(char(value), {'Force', 'Clean', 'DryRun', 'Verbose', ...
    'Shard', 'Days', 'SettleSeconds', 'Workers', 'KeepGoing', 'UseParfor'}));
end

% ---------------------------------------------------------------------------
function showConfig(opts, workers)
p = mvt.paths();
fprintf('  repo      = %s\n', p.repoRoot);
fprintf('  data      = %s\n', p.dataDir);
fprintf('  results   = %s\n', p.resultsDir);
fprintf('  days      = %s\n', mat2str(opts.Days));
fprintf('  workers   = %d\n', workers);
fprintf('  force     = %d\n', opts.Force);
fprintf('  dry run   = %d\n', opts.DryRun);
fprintf('  data ver. = %s\n', mvt.dataVersion());
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
