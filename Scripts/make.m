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

[target, workers, opts] = parseArguments(varargin);

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

fprintf('[make] target ''%s'': %s\n', target, strjoin(stages, ' '));
fprintf('[make] days %s, %d worker(s) for sharded stages\n', ...
    mat2str(days), workers);
started = tic;

for day = days
    for iStage = 1:numel(perDay)
        stage = perDay{iStage};
        if workers > 1 && ismember(stage, {'slim', 'full'})
            runStage(@() mvt.runShards(stage, day, workers, opts), stage, day, opts);
        else
            runStage(@() mvt.build(stage, day, opts), stage, day, opts);
        end
    end
end

if crossDay
    runStage(@() mvt.build('av', [], opts), 'av', [], opts);
end

fprintf('[make] target ''%s'' finished in %s\n', target, humanTime(toc(started)));
end

% ---------------------------------------------------------------------------
function runStage(thunk, stage, day, opts)
% Run one stage. Under DryRun a later stage usually cannot even look at its
% inputs, because the earlier stage that would have produced them wrote
% nothing - so a plan must report that and carry on rather than abort. A real
% run still fails loudly.
try
    thunk();
catch err
    if ~opts.DryRun
        rethrow(err);
    end
    fprintf('[make] (dry-run) %s%s: cannot plan yet (%s)\n', stage, ...
        dayLabel(day), err.identifier);
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
function [target, workers, opts] = parseArguments(args)
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
keep = true(1, numel(args));
for k = 1:2:numel(args) - 1
    if strcmpi(char(args{k}), 'Workers')
        workers = args{k + 1};
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
    'Shard', 'Days', 'SettleSeconds', 'Workers', 'UseParfor'}));
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
