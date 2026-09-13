function report = status(varargin)
% MVT.STATUS  Report which stages would rebuild, and why.
%
% Purpose
%   Answers "what will happen if I run the pipeline now?" without running
%   anything or starting a MATLAB per stage. Useful after editing a script:
%   the reason column names the file that made the target stale.
%
% Inputs
%   ...  options struct and/or name/value pairs (see mvt.options); Days selects
%        which days to report on (default [16 17 18])
%
% Outputs
%   report  struct array with fields stage, day, stale, reason. Printed as a
%           table when called without an output argument.
%
% Example
%   mvt.status('Days', 17)
%
% Dependencies
%   mvt.options, mvt.expectedOutputs, mvt.isStale, mvt.sources, mvt.dayDir
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

opts = mvt.options(varargin{:});
p = mvt.paths();

stages = { ...
    'gps',     'assemble_data_GPS'; ...
    'slim',    'generate_data_mvt_slim'; ...
    'full',    'generate_data_mvt_full'; ...
    'samples', 'generate_data_samples'; ...
    'fields',  'generate_macroscopic_fields'; ...
    'macro',   'plot_macroscopic_fields'; ...
    'micro',   'plot_microscopic_trajectories'};

report = struct('stage', {}, 'day', {}, 'stale', {}, 'reason', {});
for day = opts.Days
    for iStage = 1:size(stages, 1)
        stage = stages{iStage, 1};
        fcn = stages{iStage, 2};
        try
            outputs = mvt.expectedOutputs(stage, day, opts);
            [stale, reason] = mvt.isStale(outputs, stageInputs(stage, day, p), ...
                mvt.sources(fcn, opts), opts);
        catch err
            stale = true;
            reason = sprintf('cannot evaluate (%s)', err.message);
        end
        report(end+1) = struct('stage', stage, 'day', day, ...
            'stale', stale, 'reason', reason); %#ok<AGROW>
    end
end

% Cross-day AV figures
sampleFiles = arrayfun(@(d) fullfile(mvt.dayDir('figures', d), ...
    ['samples_for_distance_analysis_' num2str(d) '.mat']), opts.Days, ...
    'UniformOutput', false);
try
    [stale, reason] = mvt.isStale(mvt.expectedOutputs('av', opts.Days(1), opts), ...
        sampleFiles, mvt.sources('plot_AV_analysis', opts), opts);
catch err
    stale = true;
    reason = sprintf('cannot evaluate (%s)', err.message);
end
report(end+1) = struct('stage', 'av', 'day', NaN, 'stale', stale, 'reason', reason);

if nargout == 0
    fprintf('%-8s %-6s %-7s %s\n', 'STAGE', 'DAY', 'STATE', 'REASON');
    for iRow = 1:numel(report)
        if isnan(report(iRow).day)
            dayStr = 'all';
        else
            dayStr = num2str(report(iRow).day);
        end
        if report(iRow).stale
            state = 'BUILD';
        else
            state = 'fresh';
        end
        fprintf('%-8s %-6s %-7s %s\n', report(iRow).stage, dayStr, state, ...
            report(iRow).reason);
    end
    clear report
end
end

% ---------------------------------------------------------------------------
function inputs = stageInputs(stage, day, p)
rawGlob = fullfile(mvt.dayDir('raw', day), '*_0_*.json');
gpsFile = fullfile(mvt.dayDir('gps', day), ...
    sprintf('CIRCLES_GPS_10Hz_2022-11-%d.json', day));
slimGlob = fullfile(mvt.dayDir('slim', day), 'I-24MOTION_*.json');

switch stage
    case 'gps'
        inputs = {rawGlob, ...
            fullfile(p.dataDir, 'cars', 'cars_gps', 'circles_v2_1_car*.csv'), ...
            fullfile(p.dataDir, 'cars', 'cars_vins.csv'), ...
            fullfile(p.dataDir, 'cars', sprintf('veh_ping_202211%d.csv', day))};
    case {'slim', 'full'}
        inputs = {rawGlob, gpsFile};
    case {'samples', 'fields', 'micro'}
        inputs = {slimGlob};
    case 'macro'
        inputs = {fullfile(mvt.dayDir('figures', day), ...
            sprintf('fields_motion_2022-11-%d.mat', day)), gpsFile};
    otherwise
        inputs = {};
end
end
