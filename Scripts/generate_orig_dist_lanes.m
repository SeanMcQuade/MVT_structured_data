function [] = generate_orig_dist_lanes(processingDay)
% (C) 2026 CIRCLES Energy team
%
% Function that process base I-24 MOTION data from the MVT to generate
% a .mat file containing origin and distention of all trajectories in 
% the slim version of CIRCLES' v2.1 of the data, used in the team nature
% paper submission, and saves it to json files.
%
% Generated .mat files will be saved in ..\results\slim\{DATE} folder.
if nargin < 1
error(['Specify the day of Nov. 2022 MVT to generate slim'... 
        'MVT data files (from 16 to 18)']);
end
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
% Get file path of base GPS data
[parentDirectory, ~, ~] = fileparts(pwd);
% directory above contains only the git repository
[dataRootDirectory, ~, ~] = fileparts(parentDirectory);
% directory above that contains the data/ folder
dataFolderPath = fullfile(dataRootDirectory, 'data', 'i24motion', ...
    ['2022-11-', num2str(processingDay)]);

% Build the output path and filename
outputPath = fullfile(dataRootDirectory, 'results', 'slim', ...
    ['2022-11-', num2str(processingDay)]);

% Create output directory if needed
if ~isfolder(outputPath)
    mkdir(outputPath)
end

dayAbbrvs = ["mon","tue","wed","thu","fri"];
dayAbbrv = dayAbbrvs(processingDay-13);
dataFiles = dir(fullfile(dataFolderPath ,['*_' num2str(dayAbbrv) '_0_*.json']));

% avoid processing files that start with .
is_dotfile = startsWith({dataFiles.name},'.');
dataFiles = dataFiles(~is_dotfile);

if length(dataFiles) < 24
    error('I24 base files for the day: %d, Nov. 2022 are missing or incomplete.'...
        ,processingDay)
end

%========================================================================
% Process each I24 MOTION file 
%========================================================================
addpath(fullfile(parentDirectory, 'Models'));
for fileNr = 1:24 % loop over base data files
    % Load MOTION data file
    filenameLoad = fullfile(dataFolderPath,dataFiles(fileNr).name);
    fprintf('Loading and decoding MOTION data file, %d/24 ... ', ...
        fileNr); tic
    dataTemp = jsondecode(fileread(filenameLoad));
    fprintf('Done (%0.0fsec).\n',toc)

        % determine if the file already exists or not...and skip if it does
    % using dataTemp here, since it is the most recent file, to determine
    % what the output file name should be
    fileStartT = (datetime(dataTemp(1).first_timestamp, 'convertfrom', 'posixtime', ...
    'Format', 'HH:mm:ss.SSS','TimeZone' ,'America/Chicago'));
    fileStartT = datestr(fileStartT,'YYYY-mm-dd_HH-MM-SS');
    % outputFolder comes from the top of the file
    filenameSave = fullfile(outputPath,...
        ['I-24MOTION_',fileStartT,'.json']);
    % Save the processed data to a file
    
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
    for i = 1:length(dataTemp)
        dataTemp_lane_orig_dist(i).origin_lane =  dataTemp(i).origin_lane;
        dataTemp_lane_orig_dist(i).destination_lane = dataTemp(i).destination_lane;
    end
    save([filenameSave(1:end-5) '_orig_dist_lane.mat'],'dataTemp_lane_orig_dist')
   
end
rmpath(fullfile(parentDirectory, 'Models'));
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


