function [stale, reason] = isStale(outputs, inputs, sourceFiles, opts)
% MVT.ISSTALE  Make-style freshness test for a pipeline target.
%
% Purpose
%   Replaces the pipeline's original `if isfile(output), skip, end` guards.
%   Those guards never rebuilt after a code change, so editing a stage
%   silently produced nothing. This applies the classic makefile rule instead:
%   a target is stale when it is missing, when its inputs are newer, or when
%   the code that produced it is newer.
%
% Inputs
%   outputs      char / cellstr / string: files this target produces. Glob
%                patterns ('*') are expanded. Empty means "always stale".
%   inputs       (optional) data files the target is built from; globs allowed
%   sourceFiles  (optional) code files that define how it is built; typically
%                mvt.sources('<stage name>')
%   opts         (optional) options struct from mvt.options; only Force is read
%
% Outputs
%   stale   logical, true when the target must be (re)built
%   reason  char, human-readable explanation, suitable for logging in both
%           branches (e.g. 'source Scripts/generate_data_mvt_slim.m is newer
%           than output results/slim/2022-11-17/I-24MOTION_...json')
%
% Algorithm
%   1. opts.Force short-circuits to stale.
%   2. Expand the three file lists; a declared output list that resolves to
%      nothing is stale (nothing has been built yet).
%   3. Any missing output => stale.
%   4. Compare the OLDEST output timestamp against the NEWEST prerequisite
%      (inputs + sources). Newer prerequisite => stale. Using the oldest output
%      means a partially rebuilt target is correctly reported as stale.
%   5. Prerequisites that do not exist are ignored here; the stage itself
%      raises the more informative error when a required input is absent.
%
% Notes
%   Comparison uses a one-second tolerance because HFS+/APFS and network
%   volumes disagree about sub-second modification times, and MATLAB's `dir`
%   datenum is only reliable to the second.
%
% Dependencies
%   mvt.options (defaults only)
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

if nargin < 2, inputs = {}; end
if nargin < 3, sourceFiles = {}; end
if nargin < 4 || isempty(opts), opts = mvt.options(); end

if isfield(opts, 'Force') && opts.Force
    stale = true;
    reason = 'forced rebuild (Force = true)';
    return
end

declaredOutputs = asCellstr(outputs);
if isempty(declaredOutputs)
    stale = true;
    reason = 'no outputs declared';
    return
end

% ---- 1. missing outputs ----------------------------------------------------
outputFiles = {};
for iOut = 1:numel(declaredOutputs)
    spec = declaredOutputs{iOut};
    if contains(spec, '*')
        matched = expandGlob(spec);
        if isempty(matched)
            stale = true;
            reason = sprintf('no files match expected output %s', spec);
            return
        end
        outputFiles = [outputFiles, matched]; %#ok<AGROW>
    else
        if ~isfile(spec)
            stale = true;
            reason = sprintf('missing output %s', spec);
            return
        end
        outputFiles{end+1} = spec; %#ok<AGROW>
    end
end

% ---- 2. oldest output vs newest prerequisite -------------------------------
[oldestOutputTime, oldestOutputFile] = extremeTime(outputFiles, @min);

prereqs = [resolveExisting(inputs), resolveExisting(sourceFiles)];
if isempty(prereqs)
    stale = false;
    reason = 'up to date (no prerequisites declared)';
    return
end
[newestPrereqTime, newestPrereqFile] = extremeTime(prereqs, @max);

toleranceDays = 1 / 86400;   % one second
if newestPrereqTime > oldestOutputTime + toleranceDays
    stale = true;
    reason = sprintf('%s is newer than output %s', ...
        shortPath(newestPrereqFile), shortPath(oldestOutputFile));
else
    stale = false;
    reason = sprintf('up to date (oldest output %s)', shortPath(oldestOutputFile));
end
end

% ---------------------------------------------------------------------------
function out = asCellstr(value)
if isempty(value)
    out = {};
elseif ischar(value)
    out = {value};
elseif isstring(value)
    out = cellstr(value(:)');
elseif iscell(value)
    out = cellfun(@char, value(:)', 'UniformOutput', false);
else
    error('mvt:isStale:badFileList', ...
        'File lists must be char, string, or cellstr; got %s.', class(value));
end
out = out(~cellfun(@isempty, out));
end

% ---------------------------------------------------------------------------
function files = expandGlob(spec)
folder = fileparts(spec);
listing = dir(spec);
listing = listing(~[listing.isdir]);
files = cell(1, numel(listing));
for iFile = 1:numel(listing)
    files{iFile} = fullfile(folder, listing(iFile).name);
end
end

% ---------------------------------------------------------------------------
function files = resolveExisting(specs)
files = {};
specs = asCellstr(specs);
for iSpec = 1:numel(specs)
    spec = specs{iSpec};
    if contains(spec, '*')
        files = [files, expandGlob(spec)]; %#ok<AGROW>
    elseif isfile(spec)
        files{end+1} = spec; %#ok<AGROW>
    end
end
end

% ---------------------------------------------------------------------------
function [t, file] = extremeTime(files, pickFcn)
times = zeros(1, numel(files));
for iFile = 1:numel(files)
    info = dir(files{iFile});
    if isempty(info)
        times(iFile) = NaN;
    else
        times(iFile) = info(1).datenum;
    end
end
valid = ~isnan(times);
times = times(valid);
files = files(valid);
[t, idx] = pickFcn(times);
file = files{idx};
end

% ---------------------------------------------------------------------------
function short = shortPath(file)
p = mvt.paths();
short = file;
for root = {p.dataRoot}
    prefix = [root{1}, filesep];
    if startsWith(short, prefix)
        short = short(numel(prefix)+1:end);
        return
    end
end
end
