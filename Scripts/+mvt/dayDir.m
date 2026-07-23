function d = dayDir(kind, day)
% MVT.DAYDIR  Folder holding a given day's files for a given pipeline stage.
%
% Purpose
%   Centralize the '2022-11-DD' folder naming used across the pipeline so that
%   stages, the Makefile, and the tests all agree on one layout.
%
% Inputs
%   kind  char, one of:
%           'raw'      raw I-24 MOTION segments   <data>/i24motion/2022-11-DD
%           'slim'     westbound processed JSON   <results>/slim/2022-11-DD
%           'full'     full processed JSON        <results>/full/2022-11-DD
%           'figures'  .mat + figures             <results>/figures/2022-11-DD
%           'gps'      assembled GPS JSON         <results>/gps        (day-less)
%           'cache'    derived caches             <results>/.mvt/cache/2022-11-DD
%   day   numeric day of Nov. 2022 (16, 17, or 18)
%
% Outputs
%   d  char, absolute path (not created)
%
% Dependencies
%   mvt.paths
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

p = mvt.paths();
mvt.assertDay(day);
dayName = ['2022-11-', num2str(day)];

switch lower(kind)
    case 'raw'
        d = fullfile(p.dataDir, 'i24motion', dayName);
    case 'slim'
        d = fullfile(p.resultsDir, 'slim', dayName);
    case 'full'
        d = fullfile(p.resultsDir, 'full', dayName);
    case 'figures'
        d = fullfile(p.resultsDir, 'figures', dayName);
    case 'gps'
        d = fullfile(p.resultsDir, 'gps');
    case 'cache'
        d = fullfile(p.cacheDir, dayName);
    otherwise
        error('mvt:dayDir:unknownKind', ...
            ['Unknown kind ''%s''. Expected one of: raw, slim, full, ', ...
            'figures, gps, cache.'], kind);
end
end
