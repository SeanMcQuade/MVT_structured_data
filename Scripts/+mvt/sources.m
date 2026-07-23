function files = sources(stageName, opts)
% MVT.SOURCES  Code files whose modification should invalidate a stage's output.
%
% Purpose
%   Supplies the "source" half of the makefile rule in mvt.isStale: the stage's
%   own .m file plus everything it calls. Editing a fuel model must invalidate
%   the slim/full JSON; editing a plotting script must not.
%
% Inputs
%   stageName  char, function name without extension, e.g.
%              'generate_data_mvt_slim'. Must live in <repo>/Scripts.
%   opts       (optional) options struct; only Verbose is read
%
% Outputs
%   files  cellstr of absolute paths that exist on disk
%
% Algorithm
%   1. Locate <repo>/Scripts/<stageName>.m (error if absent).
%   2. Reuse the cached closure in results/.mvt/deps/<stageName>.json when it
%      was computed from the same stage-file timestamp and MATLAB release.
%   3. Otherwise ask matlab.codetools.requiredFilesAndProducts for the static
%      dependency closure, keeping only files inside the repository.
%   4. Union that with the explicit extras below, because the fuel models are
%      resolved at run time (addpath + eval of fuel_model_<class>_simplified)
%      and are therefore invisible to static analysis.
%   5. Drop this helper package (+mvt) from the closure: editing the build
%      machinery does not change scientific output, and including it would
%      invalidate every result whenever a helper is touched.
%   6. Cache the result (best effort - a read-only results tree is not fatal).
%
% Dependencies
%   mvt.paths, mvt.ensureDir
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

if nargin < 2 || isempty(opts)
    opts = mvt.options();
end

p = mvt.paths();
stageFile = fullfile(p.scriptsDir, [stageName '.m']);
if ~isfile(stageFile)
    error('mvt:sources:unknownStage', ...
        'No such stage script: %s', stageFile);
end

info = dir(stageFile);
stageStamp = info(1).datenum;
cacheFile = fullfile(p.depsDir, [stageName '.json']);

cached = readCache(cacheFile, stageStamp);
if ~isempty(cached)
    files = keepExisting(cached);
    return
end

closure = {stageFile};
try
    static = matlab.codetools.requiredFilesAndProducts(stageFile);
    if ischar(static)
        static = {static};
    end
    closure = [closure, reshape(static, 1, [])];
catch err
    mvt.log(opts, 'dependency analysis for %s fell back to the static list (%s)', ...
        stageName, err.identifier);
end

closure = [closure, extraSources(stageName, p)];
closure = keepInRepo(closure, p);
closure = dropHelperPackage(closure, p);
files = keepExisting(unique(closure, 'stable'));

writeCache(cacheFile, stageName, stageStamp, files);
end

% ---------------------------------------------------------------------------
function extras = extraSources(stageName, p)
% Run-time dependencies that static analysis cannot see.
switch stageName
    case {'generate_data_mvt_slim', 'generate_data_mvt_full'}
        % Scripts/+mvt aside, these stages addpath(Models) and dispatch to
        % fuel_model_<class>_simplified via eval, plus read the grade fit.
        listing = dir(fullfile(p.modelsDir, '*.m'));
        extras = cell(1, numel(listing) + 1);
        for iFile = 1:numel(listing)
            extras{iFile} = fullfile(p.modelsDir, listing(iFile).name);
        end
        extras{end} = fullfile(p.modelsDir, 'Eastbound_grade_fit.csv');
    otherwise
        extras = {};
end
end

% ---------------------------------------------------------------------------
function files = keepInRepo(files, p)
prefix = [p.repoRoot, filesep];
files = files(startsWith(files, prefix));
end

% ---------------------------------------------------------------------------
function files = dropHelperPackage(files, p)
helperPrefix = [fullfile(p.scriptsDir, '+mvt'), filesep];
files = files(~startsWith(files, helperPrefix));
end

% ---------------------------------------------------------------------------
function files = keepExisting(files)
files = reshape(files, 1, []);
files = files(cellfun(@isfile, files));
end

% ---------------------------------------------------------------------------
function files = readCache(cacheFile, stageStamp)
files = {};
if ~isfile(cacheFile)
    return
end
try
    cached = jsondecode(fileread(cacheFile));
catch
    return
end
if ~isfield(cached, 'stageStamp') || ~isfield(cached, 'files') || ...
        ~isfield(cached, 'release') || ~strcmp(cached.release, version('-release'))
    return
end
if abs(cached.stageStamp - stageStamp) > 1/86400
    return   % stage file changed; recompute the closure
end
files = reshape(cellstr(cached.files), 1, []);
end

% ---------------------------------------------------------------------------
function writeCache(cacheFile, stageName, stageStamp, files)
try
    mvt.ensureDir(fileparts(cacheFile));
    payload = struct( ...
        'stage', stageName, ...
        'stageStamp', stageStamp, ...
        'release', version('-release'), ...
        'files', {reshape(files, [], 1)});
    fid = fopen(cacheFile, 'w');
    if fid < 0
        return
    end
    cleanup = onCleanup(@() fclose(fid));
    fwrite(fid, jsonencode(payload), 'char');
catch
    % Caching is an optimization; never fail a run because of it.
end
end
