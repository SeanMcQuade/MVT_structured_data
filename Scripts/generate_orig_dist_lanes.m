function [] = generate_orig_dist_lanes(processingDay, varargin)
% GENERATE_ORIG_DIST_LANES  Origin and destination lane of every slim trajectory.
% (C) 2026 CIRCLES Energy team
%
% Purpose
%   Processes base I-24 MOTION data from the MVT to record, for every
%   trajectory in the slim data set, the lane it came from and the lane it went
%   to. The slim JSON itself carries only the lane the (clipped) trajectory was
%   driven in, so the lane-change analysis needs this sidecar alongside it.
%
% Inputs
%   processingDay  16, 17, or 18 (November 2022)
%   varargin       options struct and/or name/value pairs (see mvt.options);
%                  Force, Clean, DryRun, Verbose, Shard
%
% Outputs
%   <results>/analysis/2022-11-DD/
%     I-24MOTION_2022-11-DD_HH-MM-SS_orig_dist_lane.mat
%   one per raw segment, holding dataTemp_lane_orig_dist: a struct array with
%   origin_lane and destination_lane, in the same order as the trajectories in
%   the matching slim JSON, which is what lets the consumer pair them by index.
%
% Algorithm
%   1. Resolve paths and the raw-segment manifest (raw file -> output name).
%   2. Take this shard's segments (mvt.shardIndices; worker k of N takes
%      k, k+N, ...), and for each one rebuild only when the sidecar is older
%      than the raw input or the code (mvt.isStale) - checked before the
%      multi-GB decode, because the manifest already knows the output name.
%   3. Keep westbound trajectories, assign lanes, and clip lane changes; the
%      clipping records the origin and destination lane of each clipped piece.
%   4. Write the sidecar through a temporary name (mvt.atomicSave).
%
% Notes
%   FORK: the local assign_lanes and clip_lane_changes below are a fork of the
%   ones in generate_data_mvt_slim.m, extended to track origin_lane and
%   destination_lane through the clipping. The two copies must stay in step:
%   the sidecar is paired with the slim JSON *by index*, so any divergence in
%   the clipping silently mispairs lanes with trajectories rather than failing.
%   mvt.sources cannot see this coupling, so editing one copy will not
%   invalidate the other's outputs.
%   TODO: factor the shared clipping into one file both stages call, and have
%   this stage keep only the origin/destination bookkeeping.
%
% Dependencies
%   mvt.options, mvt.paths, mvt.dayDir, mvt.manifest, mvt.expectedOutputs,
%   mvt.isStale, mvt.sources, mvt.shardIndices, mvt.segmentName,
%   mvt.atomicSave, mvt.ensureDir, mvt.progress, mvt.log
%
if nargin < 1
error(['Specify the day of Nov. 2022 MVT to generate slim'...
        'MVT data files (from 16 to 18)']);
end
mvt.assertDay(processingDay)
opts = mvt.options(varargin{:});

%========================================================================
% Parameters
%========================================================================
% Set processing parameters
processingOpts.laneIdentificationOpts.LaneWidth = 12;  % Lane width in feet
% Number of cells in the x direction used to identify driving line
processingOpts.laneIdentificationOpts.Nr_XCells = 200;      
% Minimum number of samples per cell used to identify driving line
processingOpts.laneIdentificationOpts.MinCellSamples = 20; 
% samples upper bound as a multiplier of lane width to remove outliers
processingOpts.laneIdentificationOpts.Y_Up_Lim = 5;        
% samples lower bound as a multiplier of lane width to remove outliers
processingOpts.laneIdentificationOpts.Y_Low_Lim = 0.5;        
% Shift parameters for y position correction (west{east}bound) s.t. y_corrected =  
% Sw{e} * (y  - drivingLineShift_w{e}) + Cw{e};
processingOpts.laneIdentificationOpts.Se = 0.97;            
processingOpts.laneIdentificationOpts.Ce = 1;
processingOpts.laneIdentificationOpts.Sw = 0.98;          
processingOpts.laneIdentificationOpts.Cw = 1;
% Threshold to identify lane change as a multiplier of lane width
processingOpts.laneChangeClippingOpts.LaneChangeThresh = 0.5;        
% Threshold for maximum rate of change in the identified lane for which 
% trajectory is assumed to not change lane as a multiplier of lane width/s
processingOpts.laneChangeClippingOpts.MaxLaneChangeRate = 0.1;      
% Minimum duration for a clipped trajectory in seconds
processingOpts.laneChangeClippingOpts.MinClipTime = 0.5;            
% Buffer threshold to ensure full lane change as a multiplier of lane width
processingOpts.laneChangeClippingOpts.ChangeBufferThresh = 0.2; 
%========================================================================
% Initilize
%========================================================================
% Resolve the layout: the repository is a sibling of data/ and results/.
% mvt.paths derives this from the location of the code, so the stage no longer
% depends on the current folder being Scripts/.
p = mvt.paths();
parentDirectory = p.repoRoot;

% The sidecars live with the other analysis .mat products, not in slim/, which
% holds the released data set. mvt.expectedOutputs is the single declaration of
% both the folder and the 24 names.
outputPath = mvt.dayDir('analysis', processingDay);
mvt.ensureDir(outputPath)

% Map each raw segment to the file it produces, without decoding it. This is
% what lets the staleness check below run before the expensive jsondecode.
segments = mvt.manifest(processingDay, opts);
outputs = mvt.expectedOutputs('lanes', processingDay, opts);
sourceFiles = mvt.sources('generate_orig_dist_lanes', opts);

if numel(segments) < 24
    error('I24 base files for the day: %d, Nov. 2022 are missing or incomplete.'...
        ,processingDay)
end

%========================================================================
% Process each I24 MOTION file 
%========================================================================
addpath(fullfile(parentDirectory, 'Models'));
% Restores the path even when a segment throws, which the trailing rmpath the
% loop used to end with did not.
restoreModelsPath = onCleanup(@() rmpath(fullfile(parentDirectory, 'Models'))); %#ok<NASGU>
% Segments owned by this shard: worker k of N takes k, k+N, k+2N, ...
shardSegments = mvt.shardIndices(numel(segments), opts.Shard);
reportProgress = mvt.progress(numel(shardSegments), ...
    sprintf('lanes 2022-11-%d', processingDay), 'Opts', opts);
segmentsDone = 0;
for fileNr = shardSegments % loop over base data files
    segment = segments(fileNr);
    filenameLoad = segment.rawPath;
    % The output name is known from the manifest, so freshness is decided
    % before the (multi-minute, multi-GB) decode rather than after it.
    filenameSave = outputs{fileNr};
    [~, sidecarName, sidecarExt] = fileparts(filenameSave);
    sidecarName = [sidecarName sidecarExt];
    [stale, staleReason] = mvt.isStale(filenameSave, filenameLoad, sourceFiles, opts);
    if ~stale
        mvt.log(opts, 'skip %s: %s', sidecarName, staleReason);
        continue  % advance this part of the loop
    end
    mvt.log(opts, 'build %s: %s', sidecarName, staleReason);
    if opts.Clean && isfile(filenameSave) && ~opts.DryRun
        delete(filenameSave)
    end
    if opts.DryRun
        continue
    end

    % Load MOTION data file
    fprintf('Loading and decoding MOTION data file, %d/%d ... ', ...
        fileNr, numel(segments)); tic
    dataTemp = jsondecode(fileread(filenameLoad));
    fprintf('Done (%0.0fsec).\n',toc)

    % Guard against a stale manifest: the name derived from the decoded data
    % must match the one the manifest predicted.
    expectedName = mvt.segmentName(dataTemp(1).first_timestamp);
    if ~strcmp(mvt.laneSidecarName(expectedName), sidecarName)
        error('mvt:generate_orig_dist_lanes:manifestMismatch', ...
            ['Manifest predicted %s for %s but the data says %s. ', ...
            'Delete %s and re-run.'], sidecarName, segment.rawName, ...
            mvt.laneSidecarName(expectedName), fullfile(p.manifestDir, ...
            sprintf('segments_2022-11-%d.json', processingDay)));
    end

    % remove eastbound trajectories
    dataTemp = dataTemp([dataTemp.direction]<0);
    % delete extra fields
    dataTemp = rmfield(dataTemp,{'flags','compute_node_id','fragment_ids','merged_ids',...
        'configuration_id','fine_vehicle_class','x_score','y_score','road_segment_ids'});
    fprintf('Identifying lanes and clipping lane changes ...'),tic
    % Identify and assign lanes
    dataLanes = assign_lanes(dataTemp,processingOpts.laneIdentificationOpts);
    % Clip data to remove lane switches
    dataTemp = clip_lane_changes(dataTemp,dataLanes,processingOpts.laneChangeClippingOpts);
    fprintf('Done (%0.0fsec).\n',toc)
    % Preallocated fresh for every segment: reusing the variable across
    % iterations left the previous segment's trailing entries in place whenever
    % this one yielded fewer trajectories, which also made the result depend on
    % which segments a shard happened to own.
    dataTemp_lane_orig_dist = repmat( ...
        struct('origin_lane', [], 'destination_lane', []), 1, length(dataTemp));
    for i = 1:length(dataTemp)
        dataTemp_lane_orig_dist(i).origin_lane =  dataTemp(i).origin_lane;
        dataTemp_lane_orig_dist(i).destination_lane = dataTemp(i).destination_lane;
    end
    % Written through a temporary name: a sidecar truncated by an interrupted
    % run would still satisfy isfile and be taken for a finished product.
    mvt.atomicSave(filenameSave, ...
        struct('dataTemp_lane_orig_dist', dataTemp_lane_orig_dist));
    clear dataTemp dataTemp_lane_orig_dist dataLanes
    segmentsDone = segmentsDone + 1;
    reportProgress(segmentsDone, sidecarName);
end
reportProgress();
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Local Functions Definitions %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function dataLanes = assign_lanes(baseFileName,laneIdentificationOpts)

% Function that takes in original I24 v2 data and outputs lane identification at
% each time step for all trajectories and adjusted y_position (ft)
%
% Input: baseFileName: json file name [char array] or data file [struct].
%        laneIdentificationOpts: Constants set for processing lane identification (refer to head of the script)
%
% Example: assign_lanes('2022-11-14_I-24MOTION_00.json')
%
% Output: dataLanes: struct with two fields:
%       lane: assigned lane for all trajectories at each timestep
%       y_corr: y adjusted for the wiggle in the data for all trajectories  (ft)
%
%
% Lane = 5 -> vehicle is on a ramp
% Lane = 0 -> vehicle is off the highway
%

% Parse input, load data if input is a file path
if isa(baseFileName,'struct')
    data = baseFileName;
elseif ischar(baseFileName)
    fprintf(['Loading data ' baseFileName '\n'])
    filenameLoad = baseFileName;
    data = jsondecode(fileread(filenameLoad));
else
    fprintf('Wrong input \n')
    return
end
% sample data in both directions
N = length(data);
xWest = []; yWest = [];
for trajInd=1:N  % Loop over all trajectories 
    veh = data(trajInd);
    t = veh.timestamp;
    ts = t(2)-t(1); % [s] timestep
    ss = ceil(1 / ts); % sample at 1 Hz
    direction = veh.direction;
    if (t(end)-t(1))>5     %sample only trajectories > 5s , each second
        if direction == -1   %west bound
            xWest = [xWest, veh.x_position(1 :ss: end)'];
            yWest = [yWest, veh.y_position(1 :ss: end)'];
        end
    end
end
% Copy defined processing parameters 
laneWidth = laneIdentificationOpts.LaneWidth;   
nrXCells = laneIdentificationOpts.Nr_XCells;   
minCellsSamples = laneIdentificationOpts.MinCellSamples;
yUpLim = laneIdentificationOpts.Y_Up_Lim;
yLowLim = laneIdentificationOpts.Y_Low_Lim;
upBoundY = laneWidth*yUpLim;
lowBoundY = laneWidth*yLowLim;
% Remove outliers from the sampled data based on bounds
xWest = xWest(yWest<=upBoundY & yWest>=lowBoundY);
yWest = yWest(yWest<=upBoundY & yWest>=lowBoundY);
% Initialize x position grid in both directions
xCellsWest = linspace(min(xWest),max(xWest),nrXCells+1);
% Initialize driving line in both directions
drivingLineWest = zeros(1,nrXCells);
% Loop over x cells to calculate the shift in driving line
for cellInd = 1:nrXCells
    % Identify points in the cells
    pointsInCellWest = xCellsWest(cellInd)<xWest&xWest<=xCellsWest(cellInd+1);
    % Aggregate samples in the cell. Shift down so leftmost lane edge = 0 
    yCellWest = yWest(pointsInCellWest) - laneWidth/2; 
    % Calculate driving line shift if sufficient samples are found in cell
    if length(yCellWest) > minCellsSamples 
        % Mapping shift to a unit circle
        p_j_w = cos(2*pi*yCellWest/laneWidth); 
        q_j_w = sin(2*pi*yCellWest/laneWidth); 
        % Find mean shift
        dlw_t = atan2(mean(q_j_w), mean(p_j_w))*laneWidth/2/pi; 
        drivingLineWest(cellInd) = dlw_t;
    else % assume no shift if not enough samples are found in cell
        drivingLineWest(cellInd) = 0;
    end
end
% Remove maximum edge in the grid
xCellsWest = xCellsWest(1:end-1);
%%% Start Lane assignment
% Set shifting factors used for y-position correction
Sw = laneIdentificationOpts.Sw; Cw = laneIdentificationOpts.Cw;
% Loop over trajectories to correct y position and estimate lane position
dataLanes = struct([]);
for trajInd = 1: N 
    veh = data(trajInd);
    y = veh.y_position;
    x = veh.x_position;
    laneTemp = zeros(size(y)); % Pre-filtering estimated lane position
    lane = zeros(size(y)); % Output lane estimate
    n = length(laneTemp); % trajectory length
    % correct for wiggle in y
    if veh.direction <0 % westbound
        dataLanes(trajInd).y_corr = ...
            Sw * (y - interp1(xCellsWest,drivingLineWest,x,'linear','extrap'))+Cw;
    end
    % Estimate lane position from corrected y position 
    laneTemp = ((abs(dataLanes(trajInd).y_corr) - laneWidth/2) ./laneWidth) ;
    laneTemp = max(0,min(5,(laneTemp)));
    % Median filter estimated lane
    ww = 10;      %median filter window width (num. time steps)
    buff = ceil(ww/2);
    if n > ww
        for ii = buff : n-buff
            lane(ii,1) = median(laneTemp(ii-buff+1: ii+buff));
        end
        lane(1:buff-1) = lane(buff) * ones(buff-1,1);
        lane(n-buff+1:n) = lane(n-buff) * ones(buff,1);
    else
        lane = (median(laneTemp)*ones(size(laneTemp)));
    end
    dataLanes(trajInd).lane = lane;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newData = clip_lane_changes(data,dataLanes,opts)
% Function to clip v2 MVT data during lane changes and split trajectories into lanes
%
% Input: - data: original not clipped I24MOTION v2 data
%        - dataLanes: struct with two fields:
%                lane: assigned lane for all trajectories at each timestep
%                y_corr: y adjusted for the wiggle in the data for all trajectories  (ft)
%
%        - opts: struct containing constants used in clipping lanes
%
% Output:  - mewData: processed data with clipped lane changes
%

%========================================================================
% Parameters 
%========================================================================
% Threshold to identify lane change as a multiplier of lane width
LaneChangeThresh = opts.LaneChangeThresh;         
% Threshold for maximum rate of change in the identified lane for which trajectory 
% is assumed to not change lane as a multiplier of lane width
MaxLaneChangeRate = opts.MaxLaneChangeRate;     
% Minimum duration for a clipped trajectory in seconds
MinClipTime = opts.MinClipTime;                 
% Buffer threshold to ensure full lane change as a multiplier of lane width
ChangeBufferThresh = opts.ChangeBufferThresh;   
nrTrajs = length(data);
dataFields = fieldnames(data);
% Initialize struct array to for output data
newData = {};
for fieldInd = 1:length(dataFields)  % Loop over all data fields to initialize
    newData(4*nrTrajs).(cell2mat(dataFields(fieldInd))) = [];
end
newData(4*nrTrajs).lane =[];
newData(4*nrTrajs).origin_lane =[];
newData(4*nrTrajs).destination_lane =[];
trajCounter =1;
% Loop over original trajectories and clip lane changes
for trjInd = 1:nrTrajs
    veh = data(trjInd);
    x = veh.x_position;
    t = veh.timestamp;
    y = dataLanes(trjInd).y_corr;
    lane = dataLanes(trjInd).lane;
    trajectoryLength = length(veh.timestamp);
    tempLane = round(lane(1));
    pointer = 1; % pointer to loop over trajctory points 
    clippedPart = -1; %clipped part number (starting from 0)
    origin_lane = round(lane(1));
    while pointer<trajectoryLength % loop over tajectory with a pointer
        %start from when lane is stable
        stInd = find((abs(diff(lane)./diff(t))) < MaxLaneChangeRate,1); 
        if isempty(stInd)
            break %if lane is not stable at all, delete full trajectory
        end
        lane = lane(stInd:end);
        tempLane = round(lane(1));
        trajectoryLength = length(lane);
        x = x(stInd:end);
        y = y(stInd:end);
        t = t(stInd:end);        
        % Find index of lane changing
        indChange = find(abs(lane-tempLane) > LaneChangeThresh,1);
        if isempty(indChange) 
            indChange =  trajectoryLength;
            dest_lane = tempLane;
        else
            dest_lane = tempLane +  sign(lane(indChange) - tempLane);
        end
        % Remove data before lane changing (up to stable lane)
        mirroredLane0 = lane(indChange:-1:1);
        absDiffMirroredLane0 = abs(diff(mirroredLane0)./diff(t(indChange:-1:1)));
        absDiffMirroredLane0 = [absDiffMirroredLane0;absDiffMirroredLane0(end)];
        secPointer = find(abs(mirroredLane0 - tempLane) < ChangeBufferThresh &...
            absDiffMirroredLane0 < MaxLaneChangeRate,1);
        if ~isempty(secPointer) && (t(indChange-secPointer+1)-t(1)) >= MinClipTime % only MinClipTime [s] minimum
            % Log clipped trajectory 
            clippedPart = clippedPart+1;
            secPointer = indChange - secPointer;
            % Clip and append data to output struct array
            newVeh = veh;
            newVeh.timestamp = t(1:secPointer);
            newVeh.x_position = x(1:secPointer);
            newVeh.y_position = y(1:secPointer);
            newVeh.ending_x = x(secPointer);
            newVeh.first_timestamp = t(1);
            newVeh.last_timestamp = t(secPointer);
            newVeh.starting_x = x(1);
            newVeh.ending_x = x(secPointer);
            newVeh.lane = tempLane;
            newVeh.x_id.x_oid = [newVeh.x_id.x_oid '-' num2str(clippedPart)];
            newVeh.origin_lane = origin_lane;
            newVeh.destination_lane = dest_lane;

            newData(trajCounter) = newVeh;
            trajCounter = trajCounter+1;
            origin_lane = tempLane;
        end
        % Discard processed part
        lane = lane(indChange:end);
        x = x(indChange:end);
        y = y(indChange:end);
        t = t(indChange:end);
        tempLane = round(lane(1));
        % If clipping happen, remove starting part of the remaining trajectory
        % until within (ChangeBufferThresh * lane width) of new lane
        if indChange<trajectoryLength
            ii_n2 = find(abs(lane-tempLane)<ChangeBufferThresh,1);
            if isempty(ii_n2)
                ii_n2 =  length(lane);
            end
            ii_notfull = find(abs(lane-tempLane)>1-ChangeBufferThresh,1);
            if ~isempty(ii_notfull)&& ii_notfull< ii_n2  
                %false lane change: clip but don't start new trajectory
                ii_n2 =  ii_notfull;
                tempLane = round(lane(ii_n2));
            end
            pointer = ii_n2+1;
            lane = lane(pointer:end);
            x = x(pointer:end);
            y = y(pointer:end);
            t = t(pointer:end);
            pointer = 1;
            trajectoryLength = length(x);
        end
    end
end
newData(trajCounter:end) = [];
end


