function [] = generate_data_mvt_slim(processingDay)
% (C) 2025 CIRCLES Energy team
%
% Function that process base I-24 MOTION data from the MVT to generate
% a slim version of CIRCLES' v2.1 of the data, used in the team nature
% paper submission, and saves it to json files.
%
% Generated json files will be saved in .\Data\{DATE}__MVT_Data_Slim folder.
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
originXPosition = 309804.0625; % [ft] location of the origin for the x-coordinates 
ft2meterFactor = 0.3048; % [m/ft] conversion factor from feet to meter
meter2mileFactor = 6.213712e-04; % [mile/m] conversion factor from meter to mile
g2gallonsFactor = 3.522294e-04; % [gallon/g] conversion factor from fuel gram to gallon
mcDist = 0.225 ; %[mile] the distance between mill creek origin (MM58.675) and MM58.9
% the estimated origin of the road grade map
%========================================================================
% Initilize
%========================================================================
% Get file path of base GPS data
[parentDirectory, ~, ~] = fileparts(pwd);
% directory above contains only the git repository
[dataRootDirectory, ~, ~] = fileparts(parentDirectory);
% directory above that contains the data/ folder
dataFolderPath = fullfile(dataRootDirectory, 'data', '0_base', ...
    ['2022-11-', num2str(processingDay)]);

%%
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
if length(dataFiles) < 24
    error('I24 base files for the day: %d, Nov. 2022 are missing or incomplete.'...
        ,processingDay)
end

%========================================================================
% Load GPS and road grade data
%========================================================================
% Load road grade map
gradeData = readmatrix(fullfile(parentDirectory,'Models','Eastbound_grade_fit.csv'));
gradeDataStart = gradeData(:,2);
gradeDataEnd = gradeData(:,3);
gradeDataPoints = [gradeDataStart; gradeDataEnd(end)];
gradeDataSlope = gradeData(:,4);
gradeDataIntercept = gradeData(:,5);
% GPS Data is assumed to be processed. Run create_data_GPS.m to produce processed GPS files
fprintf('\nLoading and decoding AVs GPS data file ...'); tic
dataGPS = jsondecode(fileread(fullfile(dataRootDirectory,...
    'results','gps',['CIRCLES_GPS_10Hz_2022-11-' num2str(processingDay) '.json'])));
fprintf('Done (%0.0fsec).\n',toc)
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
    
    % Check if file already exists
    if isfile(filenameSave)
        fprintf('Output file already exists: %s\nSkipping processing.\n', filenameSave);
        continue  % advance this part of the loop
    else
        fprintf('File %s does not exist...starting processing.', filenameSave);
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
    fprintf('Calculating Distance to AVs ...'),tic
    distToAvsData = calculate_distance_to_avs(dataTemp,dataGPS);
    fprintf('Done (%0.0fsec).\n',toc)
    fprintf('Generating v2.1 version of data ... '),tic
    n = length(dataTemp);
    data = init_data_struct(n);
    
    for trjInd=1:floor(n) %loop over all data(vehicle) segments
        veh = dataTemp(trjInd);               % get each vehicle
        data(trjInd).trajectory_id = veh.x_id; % id for i-th vehicle trajectory
        data(trjInd).timestamp = veh.timestamp; % time array for i-th vehicle trajectory
        % numbers corresponding to vehicle types 0:sedan, 1:midSUV, 2:van, 3:pickup, 4:semi, 5:truck
        data(trjInd).coarse_vehicle_class = veh.coarse_vehicle_class; 
        data(trjInd).first_timestamp = veh.first_timestamp; % first timestamp in the trajectory fragment
        data(trjInd).last_timestamp = veh.last_timestamp; % last timestamp in the trajectory fragment
        data(trjInd).length = veh.length; % length of vehicle in feet
        data(trjInd).width = veh.width; % width of vehicle in feet
        data(trjInd).height = veh.height; % height of vehicle in feet
        %lane number of the vehicle: 1 leftmost lane(HOV) , 4 right most lane. 0 off the highway, 5 on/off ramp.
        data(trjInd).lane_number = veh.lane; 
        % shift x position to v2 local coordinate system (origin: Mill Creek Bridge) and convert to meters.
        data(trjInd).x_position_meters = ft2meterFactor*(veh.x_position-originXPosition); 
        data(trjInd).y_position_corrected_meters = veh.y_position*ft2meterFactor; % convert to meters after correction and clipping
        data(trjInd).starting_x = data(trjInd).x_position_meters(1);
        data(trjInd).ending_x = data(trjInd).x_position_meters(end);
        x = veh.x_position;          % original x-position for i-th vehicle
        x = ft2meterFactor*abs(x-x(1));      % change position from feet to meters
        t = veh.timestamp;           % original timestamp data for i-th vehicle
        t = t-t(1);                  % start from time zero for each vehicle
        % Calculate velocity (m/s)
        v = (x([2:end,end])-x([1,1:end-1]))./(t([2:end,end])-t([1,1:end-1])); 
        % Calculate acceleration central (m/s/s)
        a = (x(1:end-2)-2*x(2:end-1)+x(3:end))./((t(3:end)-t(1:end-2))/2).^2; 
        a = a([1,1:end,end]);
        data(trjInd).total_distance_traversed_meters = abs(x(end)-x(1));
        data(trjInd).speed_meters_per_second = v;
        data(trjInd).acceleration_meters_per_second_per_second = a;
        % approximate road grade theta
        % 0.225miles is the distance between mill creek origin (MM58.675) and 
        % MM58.9 the estimated origin of the road grade map
        xNew = data(trjInd).x_position_meters*meter2mileFactor - mcDist;
        cellInd = xNew*0;
        for j=1:length(gradeDataPoints)
            cellInd(xNew-gradeDataPoints(j)>=0)=j;
        end
        cellInd = min(max(cellInd,1),length(gradeDataPoints)-1);
        xNewLocal = min(max(xNew,gradeDataPoints(1)),gradeDataPoints(end));
        if veh.direction>0
            theta = asin((gradeDataSlope(cellInd).*xNewLocal(:)/100+gradeDataIntercept(cellInd)/100));
        else
            theta = -asin((gradeDataSlope(cellInd).*xNewLocal(:)/100+gradeDataIntercept(cellInd)/100));
        end
        data(trjInd).road_grade_radians = theta;
        % Identify the class of vehicles, in the original data file
        % 0=sedan 1=midsize 2=van 3=pickup 4=semi 5=truck 6=motorcycle(not present in data)
        C = veh.coarse_vehicle_class;
        switch C
            case 0
                vehType = 'midBase';
            case 1
                vehType = 'midSUV';
            case 2
                vehType = 'Pickup';
            case 3
                vehType = 'Pickup';
            case 4
                vehType = 'Class8Tractor';
            case 5
                vehType = 'Pickup';
        end
        data(trjInd).energy_model = vehType;
        % Calculate total fuel
        integrate = @(t,v) dot( t(2:end)-t(1:end-1) , (v(1:end-1)+v(2:end))/2); % quadrature
        % Fuel consumption rate of drive on road
        eval(sprintf(['[fS,~,infeasible] = fuel_model_'...
            vehType '_simplified(v,a,theta,true);']))
        % Log fuel consumption data
        data(trjInd).fuel_rate_grams_per_second = fS;
        data(trjInd).percent_infeasibility = length(infeasible(infeasible>0))/...
            length(infeasible)*100; % percentage of infeasible drive
        data(trjInd).total_fuel_consumed_grams = integrate(t,fS); % total fuel
        data(trjInd).total_fuel_consumed_gallons = ...
            g2gallonsFactor*data(trjInd).total_fuel_consumed_grams;
        data(trjInd).total_fuel_economy_mpg = ...
            ((data(trjInd).total_distance_traversed_meters)*meter2mileFactor)./...
            (data(trjInd).total_fuel_consumed_gallons);
        % Distance to av engaged US
        if ~isempty(distToAvsData(trjInd).AvIdUSEng)
            data(trjInd).upstream_engaged_av_id = distToAvsData(trjInd).AvIdUSEng;
            data(trjInd).distance_to_upstream_engaged_av_meters = -distToAvsData(trjInd).distanceToAvUSEng;
        else
            data(trjInd).upstream_engaged_av_id = [];
            data(trjInd).distance_to_upstream_engaged_av_meters = [];
        end
        % Distance to av US
        if ~isempty(distToAvsData(trjInd).AvIdUS)
            data(trjInd).upstream_av_id = distToAvsData(trjInd).AvIdUS;
            data(trjInd).distance_to_upstream_av_meters = -distToAvsData(trjInd).distanceToAvUS;
        else
            data(trjInd).upstream_av_id = [];
            data(trjInd).distance_to_upstream_av_meters = [];
        end
        % Distance to av engaged DS
        if ~isempty(distToAvsData(trjInd).AvIdDSEng)
            data(trjInd).downstream_engaged_av_id = distToAvsData(trjInd).AvIdDSEng;
            data(trjInd).distance_to_downstream_engaged_av_meters = distToAvsData(trjInd).distanceToAvDSEng;
        else
            data(trjInd).downstream_engaged_av_id = [];
            data(trjInd).distance_to_downstream_engaged_av_meters = [];
        end
        % Distance to av DS
        if ~isempty(distToAvsData(trjInd).AvIdDS)
            data(trjInd).downstream_av_id = distToAvsData(trjInd).AvIdDS;
            data(trjInd).distance_to_downstream_av_meters = distToAvsData(trjInd).distanceToAvDS;
        else
            data(trjInd).downstream_av_id = [];
            data(trjInd).distance_to_downstream_av_meters = [];
        end
        %Set number of decimal points to be saved for each processed variable
        dataFields = fields(data);        
        for iField =2:length(dataFields)
            if isnumeric(data(trjInd).(string(dataFields(iField))))
                data(trjInd).(string(dataFields(iField))) = round(data(trjInd).(...
                    string(dataFields(iField))),4,'decimals');
            end
        end
        
    end
    clear dataTemp
    pause(5)
    fprintf('Done (%0.0fsec).\n',toc)
    fprintf('Encoding and Writing json file ... '),tic
    % fileStartT = (datetime(data(1).first_timestamp, 'convertfrom', 'posixtime',...
    %     'Format', 'HH:mm:ss.SSS','TimeZone' ,'America/Chicago'));
    % fileStartT = datestr(fileStartT,'YYYY-mm-dd_HH-MM-SS');
    % filenameSave = fullfile(parentDirectory, 'Data',...
    %     ['Data_2022-11-' num2str(processingDay) '__MVT_Slim'],...
    %     ['I-24MOTION_slim_',fileStartT]);    
    jsonStr = jsonencode(data);
    clear data
    fid = fopen(filenameSave, 'w');
    fwrite(fid, jsonStr, 'char');
    fclose(fid);
    fprintf('Done (%0.0fsec).\n',toc)
    clear jsonStr
    toc
end
rmpath(fullfile(parentDirectory, 'Models'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Local Functions Definitions %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = init_data_struct(n)
data(1,n) = struct('trajectory_id', [],'timestamp',[],'x_position_meters',[],'y_position_corrected_meters',[],...
    'coarse_vehicle_class',[],'first_timestamp',[],'last_timestamp',[],...
    'starting_x',[],'ending_x',[],'length',[],'width',[],'height',[],...
    'total_distance_traversed_meters',[],'speed_meters_per_second',[],'acceleration_meters_per_second_per_second',[],...
    'road_grade_radians',[],'lane_number',[],'energy_model',[],'fuel_rate_grams_per_second',[],...
    'percent_infeasibility',[],'total_fuel_consumed_grams',[],'total_fuel_consumed_gallons',[],...
    'total_fuel_economy_mpg',[],'downstream_av_id',[],...
    'distance_to_downstream_av_meters',[],'downstream_engaged_av_id',[],...
    'distance_to_downstream_engaged_av_meters',[]);
end

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
        if isempty(indChange); indChange =  trajectoryLength; end
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
            newData(trajCounter) = newVeh;
            trajCounter = trajCounter+1;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DistToAvs = calculate_distance_to_avs(dataTemp,dataGPS)
% Function to identify closest AVs and the distance to AVs in meters for all trajectories
%
% Input: - dataTemp: I24MOTION data with clipped lane changes
%        - dataGPS: GPS data from MVT
%
% Output: - dataGPS: struct array with the following fields for each data trajectory
%             AvIdUS: nearest AV upstream at each timestamp
%             distanceToAvUS: distance to nearest AV upstream at each timestamp
%             AvIdUSEng: nearest engaged AV upstream at each timestamp
%             distanceToAvUSEng: distance to nearest engaged AV upstream at each timestamp
%             AvIdDS: nearest AV downstream at each timestamp
%             distanceToAvDS: distance to nearest AV downstream at each timestamp
%             AvIdDSEng: nearest engaged AV downstream at each timestamp
%             distanceToAvDSEng: distance to nearest engaged AV downstream at each timestamp
%

originXPosition = 309804.0625; % [ft] location of the origin for the x-coordinates 
ft2meterFactor = 0.3048; % [m/ft] conversion factor from feet to meter
% Intialize output struct array
n = length(dataTemp);
DistToAvs(1,n) = struct('distanceToAvUS',[],'AvIdUS',[],'distanceToAvUSEng',[],'AvIdUSEng',[],...
    'distanceToAvDS',[],'AvIdDS',[],'distanceToAvDSEng',[],'AvIdDSEng',[]);
% Setup field names to loop over
distKeys = ["distanceToAvUS","distanceToAvDS"]; % dist to avs fields
distEngKeys = ["distanceToAvUSEng","distanceToAvDSEng"]; % dist to engaged avs fields
idsKeys = ["AvIdUS","AvIdDS"]; % ids of avs fields
idsEngKeys = ["AvIdUSEng","AvIdDSEng"]; % ids of engaged avs fields
distDir = ["upstream","downstream"]; % position of av to ego vehicle fields
distDirSign = [-1 ,1]; % relevant dist direction fields
% Loop over trajectories in the data and calculate distance to avs
num_trajs = length(dataTemp);
for trjInd = 1:num_trajs
    veh = dataTemp(trjInd);
    tempT = veh.timestamp;
    tempX = veh.x_position;
    tempX = ft2meterFactor*(tempX-originXPosition);
    tempLane = veh.lane;
    tempDir = veh.direction;
    % Find AVs on the road concurrent to the trajectory
    avsOnRoad = find([dataGPS.assigned_lane]==tempLane & ...
        [dataGPS.direction] == tempDir & tempT(1)< [dataGPS.last_timestamp] & ...
        tempT(end)>[dataGPS.first_timestamp]);
    if ~isempty(avsOnRoad)
        for keyDir = 1:2 % evaluate for down/upstream
            sn = distDirSign(keyDir);
            key = distDir(keyDir);
            % Initialize distance to AVs to inf
            vehDistToAvs.(distDir(keyDir)) = distDirSign(keyDir)*inf*...
                (ones(length(avsOnRoad),length(tempT)));
            vehDistToAvsEngaged.(distDir(keyDir)) = distDirSign(keyDir)*...
                inf*(ones(length(avsOnRoad),length(tempT)));
            tempDistX = [];
            for onRoadInd = 1:length(avsOnRoad) % loop over trajectories and log distance
                avInd = avsOnRoad(onRoadInd);
                % Calculate distance to AV
                projAvsX = interp1(dataGPS(avInd).timestamp,dataGPS(avInd).x_position,tempT) ;
                tempDistX = (projAvsX - tempX )* tempDir;
                if ~isempty(dataGPS(avInd).control_car)
                    temp_controller_active = ...
                        min(1,max(0,round(interp1(dataGPS(avInd).timestamp,...
                        double(dataGPS(avInd).control_car),tempT)))) ;
                else
                    temp_controller_active = ...
                        min(1,max(0,round(interp1(dataGPS(avInd).timestamp,...
                        double(dataGPS(avInd).controller_engaged),tempT)))) ;
                end
                % Log calculated distance
                vehDistToAvs.(key)(onRoadInd,:) =  tempDistX.';
                vehDistToAvs.(key)(onRoadInd,isnan(tempDistX)) = sn*inf;
                vehDistToAvsEngaged.(key)(onRoadInd,:) =  tempDistX.';
                vehDistToAvsEngaged.(key)(onRoadInd,isnan(tempDistX)) = sn*inf;
                % Set distance to inf when controller is not engaged for
                % engaged distance
                vehDistToAvsEngaged.(key)(onRoadInd,temp_controller_active==0) = sn*inf;
            end
        end
        for keyDir = 1:2 % find minimum distance for down/upstream
            sn = distDirSign(keyDir);
            key = distDir(keyDir);
            % set distance to inf if av is not in the appropriate direction
            vehDistToAvs.(key)(sn*vehDistToAvs.(key)<0) = sn*inf;
            vehDistToAvsEngaged.(key)(sn*vehDistToAvsEngaged.(key)<0) = sn*inf;
            % Find AV with minimum distance to trajectory
            [nearestDist.(key),Ind.(key)] = min(sn*vehDistToAvs.(key),[],1);
            [nearestDistEngaged.(key),indEngaged.(key)] = min(sn*vehDistToAvsEngaged.(key),[],1);
            % Log output distance and AV id
            DistToAvs(trjInd).(distKeys(keyDir)) = nearestDist.(key);
            DistToAvs(trjInd).(idsKeys(keyDir)) = [dataGPS(avsOnRoad(Ind.(key))).av_id];
            infFlag = isinf(nearestDist.(key));
            if sum(infFlag) == length(infFlag)
                DistToAvs(trjInd).(distKeys(keyDir)) = [];
                DistToAvs(trjInd).(idsKeys(keyDir)) = [];
            else
                DistToAvs(trjInd).(distKeys(keyDir))(infFlag) = nan;
                DistToAvs(trjInd).(idsKeys(keyDir))(infFlag) = nan;
            end
            % Log output engaged distance and AV id
            DistToAvs(trjInd).(distEngKeys(keyDir)) = nearestDistEngaged.(key);
            DistToAvs(trjInd).(idsEngKeys(keyDir)) = [dataGPS(avsOnRoad(indEngaged.(key))).av_id];
            infFlag = isinf(nearestDistEngaged.(key));
            if sum(infFlag) == length(infFlag)
                DistToAvs(trjInd).(distEngKeys(keyDir)) = [];
                DistToAvs(trjInd).(idsEngKeys(keyDir)) = [];
            else
                DistToAvs(trjInd).(distEngKeys(keyDir))(infFlag) = nan;
                DistToAvs(trjInd).(idsEngKeys(keyDir))(infFlag) = nan;
            end
        end
    end
end
end


