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
%   report  struct array with fields stage, day, stale, optional, reason.
%           Printed as a table when called without an output argument.
%           `optional` marks a stage that no default target builds - only
%           `full`, and only while it has never been built - so a caller can
%           tell "nothing will build this" from "this is pending".
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
    'lanes',   'generate_orig_dist_lanes'; ...
    'lc',      'extract_lane_changes_v_dist_to_av'; ...
    'relspeed','relative_speed_histogram'; ...
    'lcplot',  'plotting_LC_analysis'; ...
    'relspeedplot', 'plot_relative_speed'; ...
    'samples', 'generate_data_samples'; ...
    'fields',  'generate_macroscopic_fields'; ...
    'macro',   'plot_macroscopic_fields'; ...
    'micro',   'plot_microscopic_trajectories'};

report = struct('stage', {}, 'day', {}, 'stale', {}, 'optional', {}, ...
    'blocked', {}, 'reason', {});
for day = opts.Days
    for iStage = 1:size(stages, 1)
        stage = stages{iStage, 1};
        fcn = stages{iStage, 2};
        optional = false;
        blocked = false;
        try
            outputs = mvt.expectedOutputs(stage, day, opts);
            inputs = stageInputs(stage, day, p);
            [stale, reason] = mvt.isStale(outputs, inputs, ...
                mvt.sources(fcn, opts), opts);
            % A stage whose inputs are simply absent is not "stale": there is
            % nothing it could do. Saying so, and naming the file, is what tells
            % someone their download is incomplete rather than broken.
            absent = missingInputs(inputs);
            if stale && ~isempty(absent)
                blocked = true;
                reason = sprintf('needs %s', strjoin(absent, ', '));
            end
            % `full` is not part of `all` in either make implementation, so a
            % missing `full` tree is the normal state, not work that is pending.
            % Reporting it as BUILD alongside stages that `all` really does
            % build invents 24 missing files a day that nothing will ever
            % create. Once someone has built it, staleness matters again and is
            % reported as usual.
            if strcmp(stage, 'full') && ~any(cellfun(@isfile, outputs))
                optional = true;
                reason = 'not built by ''all''; opt in with `make full`';
            end
        catch err
            stale = true;
            switch err.identifier
                case {'mvt:expectedOutputs:noManifestNoProduct', ...
                        'mvt:manifest:noRawFolder'}
                    % No raw recordings and no product to describe instead:
                    % this stage is simply not available on this download.
                    blocked = true;
                    reason = 'needs the raw I-24 MOTION recordings';
                otherwise
                    reason = sprintf('cannot evaluate (%s)', err.message);
            end
        end
        report(end+1) = struct('stage', stage, 'day', day, ...
            'stale', stale, 'optional', optional, 'blocked', blocked, ...
            'reason', reason); %#ok<AGROW>
    end
end

% Cross-day AV figures
sampleFiles = arrayfun(@(d) fullfile(mvt.dayDir('analysis', d), ...
    ['samples_for_distance_analysis_' num2str(d) '.mat']), opts.Days, ...
    'UniformOutput', false);
try
    [stale, reason] = mvt.isStale(mvt.expectedOutputs('av', opts.Days(1), opts), ...
        sampleFiles, mvt.sources('plot_AV_analysis', opts), opts);
catch err
    stale = true;
    reason = sprintf('cannot evaluate (%s)', err.message);
end
report(end+1) = struct('stage', 'av', 'day', NaN, 'stale', stale, ...
    'optional', false, 'blocked', false, 'reason', reason);

if nargout == 0
    fprintf('%-8s %-6s %-7s %s\n', 'STAGE', 'DAY', 'STATE', 'REASON');
    for iRow = 1:numel(report)
        if isnan(report(iRow).day)
            dayStr = 'all';
        else
            dayStr = num2str(report(iRow).day);
        end
        if report(iRow).optional
            state = 'opt-in';
        elseif report(iRow).blocked
            state = 'needs';
        elseif report(iRow).stale
            state = 'BUILD';
        else
            state = 'fresh';
        end
        fprintf('%-8s %-6s %-7s %s\n', report(iRow).stage, dayStr, state, ...
            report(iRow).reason);
    end
    if any([report.blocked])
        fprintf(['\n''needs'' means an input is not in this download, not that ', ...
            'anything is broken.\nThose stages will be skipped; the rest build ', ...
            'normally (use `make -k`).\nSee the download routes in README.md.\n']);
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
    case {'samples', 'fields', 'micro', 'relspeed'}
        inputs = {slimGlob};
    case 'relspeedplot'
        inputs = {fullfile(mvt.dayDir('analysis', day), ...
            sprintf('relspeed_data_%d.mat', day))};
    case 'lanes'
        inputs = {rawGlob};
    case 'lc'
        inputs = {slimGlob, fullfile(mvt.dayDir('analysis', day), ...
            'I-24MOTION_*_orig_dist_lane.mat')};
    case 'lcplot'
        inputs = {fullfile(mvt.dayDir('analysis', day), ...
            sprintf('LC_data_%d.mat', day)), gpsFile};
    case 'macro'
        inputs = {fullfile(mvt.dayDir('analysis', day), ...
            sprintf('fields_motion_2022-11-%d.mat', day)), gpsFile};
    otherwise
        inputs = {};
end
end

% ---------------------------------------------------------------------------
function absent = missingInputs(inputs)
% Declared inputs that are not on disk. Globs count as present when they match
% anything, because a stage that reads "the day's segments" needs some, not a
% particular one.
absent = {};
p = mvt.paths();
for iSpec = 1:numel(inputs)
    spec = inputs{iSpec};
    if contains(spec, '*')
        listing = dir(spec);
        found = ~isempty(listing(~[listing.isdir]));
    else
        found = isfile(spec);
    end
    if ~found
        short = spec;
        prefix = [p.dataRoot, filesep];
        if startsWith(short, prefix)
            short = short(numel(prefix)+1:end);
        end
        absent{end+1} = short; %#ok<AGROW>
    end
end
end
