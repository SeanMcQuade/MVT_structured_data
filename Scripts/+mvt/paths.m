function p = paths()
% MVT.PATHS  Canonical locations of the MVT code, data, and results trees.
%
% Purpose
%   Single source of truth for every path used by the pipeline. Resolves the
%   layout from the location of this file rather than from `pwd`, so callers
%   work regardless of the current folder while still honoring the historical
%   convention that the repository is a *sibling* of `data/` and `results/`.
%
% Inputs
%   (none)
%
% Outputs
%   p  struct with fields:
%        scriptsDir   <repo>/Scripts
%        repoRoot     <repo>                     (MVT_structured_data)
%        modelsDir    <repo>/Models
%        dataRoot     parent of <repo>           (workspace root)
%        dataDir      <dataRoot>/data            (raw inputs)
%        resultsDir   <dataRoot>/results         (generated artifacts)
%        stateDir     <resultsDir>/.mvt          (build state; never released)
%        cacheDir     <stateDir>/cache           (derived caches, e.g. *_reduced.mat)
%        stampDir     <stateDir>/stamps          (make stamp files)
%        manifestDir  <stateDir>/manifests       (input -> output name maps)
%        depsDir      <stateDir>/deps            (cached dependency closures)
%
% Algorithm
%   1. Take the folder holding this file (<repo>/Scripts/+mvt).
%   2. Walk up: +mvt -> Scripts -> repo root -> workspace root.
%   3. Derive the data/results/state paths from the workspace root, honoring
%      the MVT_DATA_DIR and MVT_RESULTS_DIR overrides.
%
% Environment overrides
%   MVT_RESULTS_DIR  write results somewhere other than <dataRoot>/results.
%                    Use it to generate a second copy of the outputs for
%                    comparison against an existing run:
%                      MVT_RESULTS_DIR=/path/results_new make slim-16
%                    A relative value is taken relative to the workspace root.
%   MVT_DATA_DIR     read raw inputs from somewhere other than <dataRoot>/data
%
% Dependencies
%   None (base MATLAB).
%
% Notes
%   Nothing is created here; use mvt.ensureDir to make directories. The code
%   locations are cached, but the overrides are re-read on every call so that
%   a script can point at a different results tree mid-session.
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

persistent codePaths
if isempty(codePaths)
    here = fileparts(mfilename('fullpath'));   % <repo>/Scripts/+mvt
    scriptsDir = fileparts(here);              % <repo>/Scripts
    repoRoot = fileparts(scriptsDir);          % <repo>
    dataRoot = fileparts(repoRoot);            % workspace root
    codePaths = struct( ...
        'scriptsDir', scriptsDir, ...
        'repoRoot', repoRoot, ...
        'modelsDir', fullfile(repoRoot, 'Models'), ...
        'dataRoot', dataRoot);
end

p = codePaths;
p.dataDir = resolveOverride('MVT_DATA_DIR', p.dataRoot, 'data');
p.resultsDir = resolveOverride('MVT_RESULTS_DIR', p.dataRoot, 'results');
p.stateDir = fullfile(p.resultsDir, '.mvt');
p.cacheDir = fullfile(p.stateDir, 'cache');
p.stampDir = fullfile(p.stateDir, 'stamps');
p.manifestDir = fullfile(p.stateDir, 'manifests');
p.depsDir = fullfile(p.stateDir, 'deps');
end

% ---------------------------------------------------------------------------
function folder = resolveOverride(variableName, root, defaultName)
folder = strtrim(getenv(variableName));
if isempty(folder)
    folder = fullfile(root, defaultName);
elseif ~startsWith(folder, filesep) && ~contains(folder, ':')
    % Relative overrides are resolved against the workspace root.
    folder = fullfile(root, folder);
end
end
