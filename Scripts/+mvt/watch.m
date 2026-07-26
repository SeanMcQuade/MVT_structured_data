function report = watch(varargin)
% MVT.WATCH  Live progress view of a pipeline run, from its outputs alone.
%
% Purpose
%   Answers "how far along is it?" while `make -j` is running. The Makefile
%   owns the process tree and each stage is its own MATLAB, so their stdout
%   interleaves and no single process can render a shared progress bar. This
%   watcher sidesteps that: it polls the files the stages are expected to
%   produce and reports completion from what exists on disk.
%
%   That is deliberately non-invasive. Instrumenting the stage scripts with
%   heartbeats would give finer detail, but editing a stage makes its outputs
%   stale under the mvt.isStale rule, so adding progress reporting would force
%   a rebuild of everything already computed. Watching outputs costs nothing.
%
% Inputs
%   Name/value pairs:
%     'Days'      days to watch (default [16 17 18]; also honors MVT_DAYS)
%     'Interval'  seconds between refreshes (default 2)
%     'Once'      true to print one snapshot and return (default false).
%                 Use this when logging to a file: it prints a plain block
%                 with no screen clearing.
%     'Stages'    stages to watch (default gps, slim, samples, fields, macro,
%                 micro - 'full' is omitted, see splitArgs)
%   Any other name/value pair is passed through to mvt.options.
%
% Outputs
%   report  struct array (stage, day, done, total, active, fraction), returned
%           when an output is requested. Printed otherwise.
%
% Example
%   mvt.watch                      % live, all three days, refresh every 2 s
%   mvt.watch('Days', 18)          % one day
%   mvt.watch('Once', true)        % single snapshot, log-friendly
%
% Notes
%   A target counts as "active" when one of its outputs changed within the
%   last ActiveSeconds (default 30) - that is a heuristic, not a process
%   check, so a stage that is thinking hard between writes reads as idle.
%
% Dependencies
%   mvt.options, mvt.expectedOutputs, mvt.assertDay
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

[interval, once, activeSeconds, stages, passthrough] = splitArgs(varargin);
opts = mvt.options(passthrough{:});

while true
    snap = snapshot(stages, opts, activeSeconds);
    if nargout > 0
        report = snap;      % assigned only when asked for, so `mvt.watch`
        if once              % as a statement prints the view and nothing else
            return
        end
    end
    render(snap, once);
    if once
        break
    end
    pause(interval);
end
end

% ---------------------------------------------------------------------------
function report = snapshot(stages, opts, activeSeconds)
report = struct('stage', {}, 'day', {}, 'done', {}, 'total', {}, ...
    'active', {}, 'fraction', {});
now_ = now * 86400; %#ok<TNOW1> % datenum days -> seconds, matching dir().datenum

for day = opts.Days
    for iStage = 1:numel(stages)
        stage = stages{iStage};
        try
            expected = mvt.expectedOutputs(stage, day, opts);
        catch
            continue        % a stage that cannot declare outputs is skipped
        end
        [done, total, newest] = countOutputs(expected);
        report(end+1) = struct('stage', stage, 'day', day, ...
            'done', done, 'total', total, ...
            'active', ~isempty(newest) && (now_ - newest) < activeSeconds, ...
            'fraction', done / max(total, 1)); %#ok<AGROW>
    end
end
end

% ---------------------------------------------------------------------------
function [done, total, newest] = countOutputs(expected)
% Count satisfied outputs. Each entry counts as exactly one item, whether it
% is an exact name or a glob.
%
% Counting a glob's *matches* would be wrong in both directions: the expected
% count is unknown (that is why it is a glob), so n matches would always report
% n/n = complete, and unrelated files that happen to match would inflate it.
% One item, satisfied or not, is the only honest reading - so `slim`, whose 24
% names come from the manifest, is the target that shows real partial progress.
done = 0; total = numel(expected); newest = [];
for iFile = 1:numel(expected)
    listing = dir(expected{iFile});
    if ~isempty(listing)
        listing = listing(~[listing.isdir]);
        listing = listing(~startsWith({listing.name}, '.'));
    end
    if ~isempty(listing)
        done = done + 1;
        newest = newestOf(newest, listing);
    end
end
end

% ---------------------------------------------------------------------------
function newest = newestOf(newest, listing)
if isempty(listing)
    return
end
candidate = max([listing.datenum]) * 86400;
newest = max([newest, candidate]);
end

% ---------------------------------------------------------------------------
function render(report, once)
if ~once
    clc
end
fprintf('MVT pipeline progress  (%s)\n\n', datestr(now, 'HH:MM:SS')); %#ok<TNOW1,DATST>

if isempty(report)
    fprintf('  nothing to report\n');
    return
end

totalDone = sum([report.done]);
totalWanted = sum([report.total]);
lastDay = -1;
for iRow = 1:numel(report)
    row = report(iRow);
    if row.day ~= lastDay
        fprintf('\n  2022-11-%d\n', row.day);
        lastDay = row.day;
    end
    marker = '  ';
    if row.active
        marker = '->';
    elseif row.done >= row.total && row.total > 0
        marker = 'ok';
    end
    fprintf('    %s %-8s %s %3d/%-3d\n', marker, row.stage, ...
        bar(row.fraction, 24), row.done, row.total);
end

fprintf('\n  overall %s %d/%d (%.0f%%)\n', bar(totalDone / max(totalWanted, 1), 30), ...
    totalDone, totalWanted, 100 * totalDone / max(totalWanted, 1));
if ~once
    fprintf('\n  Ctrl-C to stop watching (the pipeline keeps running).\n');
end
end

% ---------------------------------------------------------------------------
function text = bar(fraction, width)
fraction = min(1, max(0, fraction));
filled = round(fraction * width);
text = ['[' repmat('#', 1, filled) repmat('.', 1, width - filled) ']'];
end

% ---------------------------------------------------------------------------
function [interval, once, activeSeconds, stages, passthrough] = splitArgs(args)
% Pull the watcher's own options out; everything else goes to mvt.options,
% which errors on names it does not know.
interval = 2;
once = false;
activeSeconds = 30;
% 'full' is omitted: it is an optional stage (eastbound/reference data the
% paper does not use), so watching it would show 0/24 forever on a normal run.
% Pass 'Stages' to include it.
stages = {'gps', 'slim', 'samples', 'fields', 'macro', 'micro'};
passthrough = {};

iArg = 1;
while iArg <= numel(args)
    name = args{iArg};
    if (ischar(name) || isstring(name)) && iArg < numel(args)
        switch lower(char(name))
            case 'interval'
                interval = args{iArg+1}; iArg = iArg + 2; continue
            case 'once'
                once = logical(args{iArg+1}); iArg = iArg + 2; continue
            case 'activeseconds'
                activeSeconds = args{iArg+1}; iArg = iArg + 2; continue
            case 'stages'
                stages = cellstr(args{iArg+1}); iArg = iArg + 2; continue
        end
    end
    passthrough{end+1} = args{iArg}; %#ok<AGROW>
    iArg = iArg + 1;
end
end
