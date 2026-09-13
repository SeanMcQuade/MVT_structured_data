function accepted = accept(varargin)
% MVT.ACCEPT  Mark existing outputs as up to date with the current code.
%
% Purpose
%   Staleness is decided from modification times, so editing a stage marks its
%   outputs for rebuild - correct in general, but wrong when the edit provably
%   does not change results (a refactor, a comment, a new option). Rebuilding
%   the full three-day tree costs hours and ~150 GB of writes.
%
%   This updates the timestamps of outputs that already exist, so the pipeline
%   treats them as current. It never creates, deletes, or modifies content.
%
% Inputs
%   varargin  options struct and/or name/value pairs (see mvt.options);
%             Days selects the days, DryRun reports without touching
%
% Outputs
%   accepted  cellstr of files whose timestamps were updated (or would be,
%             under DryRun)
%
% Usage
%   mvt.accept                        % all stages, all days
%   mvt.accept('Days', 17)            % one day
%   mvt.accept('DryRun', true)        % list what would be accepted
%
% Warning
%   Only run this when you have established that the current code reproduces
%   the existing outputs - for example by rebuilding one segment into a
%   separate results tree (MVT_RESULTS_DIR) and comparing checksums.
%
% Dependencies
%   mvt.options, mvt.expectedOutputs, mvt.dayDir, mvt.log
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

opts = mvt.options(varargin{:});
stages = {'gps', 'slim', 'full', 'samples', 'fields', 'macro', 'micro'};
accepted = {};

for day = opts.Days
    for iStage = 1:numel(stages)
        try
            outputs = mvt.expectedOutputs(stages{iStage}, day, opts);
        catch err
            mvt.log(opts, 'skipping %s for 2022-11-%d (%s)', stages{iStage}, day, err.message);
            continue
        end
        accepted = [accepted, touchAll(outputs, opts)]; %#ok<AGROW>
    end
end

% Cross-day figures
accepted = [accepted, touchAll(mvt.expectedOutputs('av', opts.Days(1), opts), opts)];

mvt.log(opts, 'accepted %d existing output files as current', numel(accepted));
end

% ---------------------------------------------------------------------------
function touched = touchAll(specs, opts)
touched = {};
for iSpec = 1:numel(specs)
    spec = specs{iSpec};
    folder = fileparts(spec);
    listing = dir(spec);
    listing = listing(~[listing.isdir]);
    for iFile = 1:numel(listing)
        target = fullfile(folder, listing(iFile).name);
        if opts.DryRun
            mvt.log(opts, 'would accept %s', target);
        else
            touchFile(target);
        end
        touched{end+1} = target; %#ok<AGROW>
    end
end
end

% ---------------------------------------------------------------------------
function touchFile(target)
if ispc
    command = sprintf('powershell -NoProfile -Command "(Get-Item ''%s'').LastWriteTime = Get-Date"', target);
else
    command = sprintf('touch -c ''%s''', target);
end
[status, message] = system(command);
if status ~= 0
    error('mvt:accept:touchFailed', 'Could not update the timestamp of %s: %s', ...
        target, message);
end
end
