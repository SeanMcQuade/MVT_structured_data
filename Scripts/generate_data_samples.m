function [] = generate_data_samples(processingDay, varargin)
% GENERATE_DATA_SAMPLES  Collect per-sample data used by the fuel-consumption figures.
%
% Purpose
%   Walks one day of processed MOTION trajectories and collects, for every
%   sample within MAX_DIST of an engaged CIRCLES control vehicle, the distance
%   to that vehicle together with speed, fuel rate, fuel consumption, vehicle
%   class, lane, position, and time. plot_AV_analysis bins these samples to
%   produce the article's fuel-consumption-versus-distance results.
%
% Inputs
%   processingDay  16, 17, or 18 (November 2022)
%   varargin       options struct and/or name/value pairs (see mvt.options);
%                  Force, Clean, DryRun, Verbose
%   Files:
%     <results>/slim/2022-11-DD/I-24MOTION_*.json   (falls back to full/)
%
% Outputs
%   <results>/figures/2022-11-DD/samples_for_distance_analysis_DD.mat
%     variables samples_dist, samples_speed, samples_fr, samples_fcons,
%     samples_class, samples_xpos, samples_lane, samples_t   (saved as -v7.3)
%
% Algorithm
%   1. Prefer the slim data set; fall back to full when fewer than 24 slim
%      segments are present.
%   2. For each segment, skip files with no engaged AV present.
%   3. For every trajectory sample, compute signed distance to the nearest
%      engaged AV (upstream negative, downstream positive) and keep samples
%      within Max_Dist (1000 m).
%   4. Append per-file samples into day-level column vectors and save.
%
% Parallel safety
%   Single output file per day; not shardable. Days are independent.
%
% Dependencies
%   mvt.options, mvt.paths, mvt.dayDir, mvt.isStale, mvt.sources,
%   mvt.ensureDir, mvt.log
%
% (C) 2025-2026 CIRCLES Consortium. Authors: Sulaiman Almatrudi (energy team),
% adapted by Sean McQuade (scenario team). BSD-3-Clause.
if nargin < 1
    error(['Specify the day of Nov. 2022 MVT to collect samples from '...
        'MVT data files (from 16 to 18)']);
end
mvt.assertDay(processingDay)
opts = mvt.options(varargin{:});
%========================================================================
% Generate samples_for_distance_analysis
%========================================================================
% Get file path of slim data. mvt.paths resolves the layout from the location
% of the code, so this no longer depends on the current folder being Scripts/.
p = mvt.paths();
parentDirectory = p.repoRoot; %#ok<NASGU> % retained for local edits/debugging
dataRootDirectory = p.dataRoot;
dataFolderPath = mvt.dayDir('slim', processingDay);

I24FilesInDir = dir(fullfile(dataFolderPath, ...
    ['I-24*' char(num2str(processingDay)) '*.json']));
nrFiles = length(I24FilesInDir);
if nrFiles < 24
    % try the full data if slim is not available
    dataFolderPath = mvt.dayDir('full', processingDay);
    I24FilesInDir = dir(fullfile(dataFolderPath, ...
        ['I-24*' char(num2str(processingDay)) '*.json']));
    nrFiles = length(I24FilesInDir);
    if nrFiles < 24
        error(['Processed I-24 data files not found. Please generate files'...
            ' by running generate_data_mvt_slim.m first.'])
    end
end

% create the save output folder if it does not already exist

outputPath = mvt.dayDir('figures', processingDay);

mvt.ensureDir(outputPath)

% Rebuild when the output is missing, older than the processed segments it
% summarizes, or older than the code that produced it; Force overrides.
filenameSave = fullfile(outputPath, ['samples_for_distance_analysis_' char(num2str(processingDay)) ...
    '.mat']);
[stale, staleReason] = mvt.isStale(filenameSave, ...
    fullfile(dataFolderPath, ['I-24*' char(num2str(processingDay)) '*.json']), ...
    mvt.sources('generate_data_samples', opts), opts);
if ~stale
    mvt.log(opts, 'skip %s: %s', ...
        ['samples_for_distance_analysis_' char(num2str(processingDay)) '.mat'], staleReason);
    return % return fr this function
end
mvt.log(opts, 'build samples_for_distance_analysis_%d.mat: %s', ...
    processingDay, staleReason);
if opts.Clean && isfile(filenameSave) && ~opts.DryRun
    delete(filenameSave)
end
if opts.DryRun
    return
end

%% script outputs
samples_dist = []; %distance relative to engaged av (positive = behind an av, negative = ahead of an av)
samples_speed = []; %speed samples with downstream/upstream av
samples_fr = []; %fuel rate (g/s) samples with downstream/upstream av
samples_fcons = []; %fuel consumption (g/m) samples with downstream/upstream av
samples_class = []; %vehicle class for samples with downstream/upstream  av
samples_xpos = []; %x position for samples with downstream/upstream  av
samples_lane = []; %lane number for samples with downstream/upstream  av
samples_t = []; %timestamp for samples with downstream/upstream  av (seconds after 6am)
%% data collection
for file_nr =1:nrFiles % loop over data files and append samples
    fprintf('Loading original data file %d / %d ...', file_nr,nrFiles )
    tic
    % Load MOTION data file
    filenameLoad = fullfile(dataFolderPath , I24FilesInDir(file_nr).name);
    data = jsondecode(fileread(filenameLoad));
    fprintf('Done (%0.0fsec).\n',toc)
    % Skip file if no AVs are active
    if ~isfield(data,'distance_to_upstream_engaged_av_meters')
        fprintf('no AVs on the road\n');
        continue;
    end
    fprintf('Processing data...')
    tic
    % collect samples from trajectories
    [dist, speed, class, fr, xpos, t, lane, fcons] = stats_to_av_dist(data,processingDay);
    fprintf('Done (%0.0fsec).\n',toc)
    % Append samples from file to day samples
    fprintf('Appending data...');tic
    samples_dist = [samples_dist;dist];
    samples_speed = [samples_speed;speed];
    samples_class = [samples_class;class];
    samples_fr = [samples_fr;fr];
    samples_xpos = [samples_xpos;xpos];
    samples_t = [samples_t;t];
    samples_lane = [samples_lane;lane];
    samples_fcons = [samples_fcons;fcons];
    fprintf('Done (%0.0fsec).\n',toc)
end
% filenameSave comes from up above
save(filenameSave,'samples_*','-v7.3')
end
%%%%%%%%%%%%%%%%%%%%%% local function %%%%%%%%%%%%%%%%%%%%%%
function [all_d, all_v, all_vc, all_fr, all_x, all_t, all_lane, all_fcons] = ...
    stats_to_av_dist(dataTemp,dayDate)
sixAM18 = 1668772800; % [s] epoch time corresponding to 6 am 11/18/2022 EDT
Max_Dist = 1000; %[m] maximum distance from an engaged AV for samples to be collected
% initialize with appropriate data type to optimize for memory use
all_v = inf*ones(size(vertcat(dataTemp.timestamp)));
all_fr = inf*ones(size(vertcat(dataTemp.timestamp)));
all_d = inf*ones(size(vertcat(dataTemp.timestamp)));
all_x = int16(inf*ones(size(vertcat(dataTemp.timestamp))));
all_vc = uint8(all_x);
all_lane = uint8(all_x);
all_t = uint16(all_x);
sampleCounter = 1;
for vehInd=1:length(dataTemp)
    veh = dataTemp(vehInd);
    if isfield(veh,'direction')
        % only aggregate for westbound traffic
        if veh.direction>0
            continue
        end
    end
    nemp_p = ~isempty(veh.distance_to_downstream_engaged_av_meters);
    nemp_n = ~isempty(veh.distance_to_upstream_engaged_av_meters);
    if nemp_p || nemp_n
        % if there's an active AV up or downstream of trajectory
        for pointInd =1:length(veh.timestamp) % loop over all points in the trajectory
            %aggregate data for points that have an AV downstream within Max_Dist
            if nemp_p && (~isnan(veh.distance_to_downstream_engaged_av_meters(pointInd))) && ...
                    (veh.distance_to_downstream_engaged_av_meters(pointInd) <= Max_Dist)
                % append sample to appropriate vectors
                temp_d = veh.distance_to_downstream_engaged_av_meters(pointInd);
                all_v(sampleCounter) = veh.speed_meters_per_second(pointInd);
                all_fr(sampleCounter) = veh.fuel_rate_grams_per_second(pointInd);
                all_vc(sampleCounter) = uint8(veh.coarse_vehicle_class);
                all_x(sampleCounter) = int16(veh.x_position_meters(pointInd));
                all_d(sampleCounter) =   temp_d;
                all_lane(sampleCounter) = veh.lane_number;
                %record time in (s) after 6am of each test day
                all_t(sampleCounter) = uint16(veh.timestamp(pointInd)-(sixAM18-(18-dayDate)*24*60*60));
                sampleCounter =sampleCounter+1;
            end
            %aggregate data for points that have an AV upstream within Max_Dist
            if nemp_n && (~isnan(veh.distance_to_upstream_engaged_av_meters(pointInd))) && ...
                    (veh.distance_to_upstream_engaged_av_meters(pointInd) >= -Max_Dist)
                % append sample to appropriate vectors
                temp_d = veh.distance_to_upstream_engaged_av_meters(pointInd);
                all_v(sampleCounter) = veh.speed_meters_per_second(pointInd);
                all_fr(sampleCounter) = veh.fuel_rate_grams_per_second(pointInd);
                all_vc(sampleCounter) = uint8(veh.coarse_vehicle_class);
                all_x(sampleCounter) = int16(veh.x_position_meters(pointInd));
                all_d(sampleCounter) =   temp_d;
                all_lane(sampleCounter) = veh.lane_number;
                %record time in (s) after 6am of each test day
                all_t(sampleCounter) = uint16(veh.timestamp(pointInd)-(sixAM18-(18-dayDate)*24*60*60));
                sampleCounter =sampleCounter+1;
            end
        end
    end
end
% claer extra rserved size in all vectors
all_d(sampleCounter:end) = [];
all_v(sampleCounter:end) = [];
all_vc(sampleCounter:end) = [];
all_fr(sampleCounter:end) = [];
all_x(sampleCounter:end) = [];
all_t(sampleCounter:end) = [];
all_lane(sampleCounter:end) = [];
% calculate instantanious fuel consumption (add 1e-6 to avoid division by 0)
all_fcons = all_fr./ (1e-6+all_v);
end
