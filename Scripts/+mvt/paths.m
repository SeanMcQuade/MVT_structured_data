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
%   3. Derive the data/results/state paths from the workspace root.
%   4. Cache the result in a persistent variable (paths cannot change within a
%      MATLAB session).
%
% Dependencies
%   None (base MATLAB).
%
% Notes
%   Nothing is created here; use mvt.ensureDir to make directories.
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

persistent cached
if isempty(cached)
    here = fileparts(mfilename('fullpath'));   % <repo>/Scripts/+mvt
    scriptsDir = fileparts(here);              % <repo>/Scripts
    repoRoot = fileparts(scriptsDir);          % <repo>
    dataRoot = fileparts(repoRoot);            % workspace root
    resultsDir = fullfile(dataRoot, 'results');
    stateDir = fullfile(resultsDir, '.mvt');
    cached = struct( ...
        'scriptsDir', scriptsDir, ...
        'repoRoot', repoRoot, ...
        'modelsDir', fullfile(repoRoot, 'Models'), ...
        'dataRoot', dataRoot, ...
        'dataDir', fullfile(dataRoot, 'data'), ...
        'resultsDir', resultsDir, ...
        'stateDir', stateDir, ...
        'cacheDir', fullfile(stateDir, 'cache'), ...
        'stampDir', fullfile(stateDir, 'stamps'), ...
        'manifestDir', fullfile(stateDir, 'manifests'), ...
        'depsDir', fullfile(stateDir, 'deps'));
end
p = cached;
end
