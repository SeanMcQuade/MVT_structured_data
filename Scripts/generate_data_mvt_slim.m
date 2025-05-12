% save_i24motiondata_v2_point_1_slim - (C) CIRCLES Energy team, 4/15/2025
%
% This script process v2.0 of the I-24 MOTION data from the MVT to generate
% a slim version of CIRCLES' v2.1 of the data, used in the team nature
% paper submission, and saves it to json files.
%
% Generated json files will be saved in \{DATE}__MVT_Data_Slim folder.

clearvars -except DAY_TO_PROCESS
if ~exist('DAY_TO_PROCESS', 'var')
    global DAY_TO_PROCESS
    DAY_TO_PROCESS = input('Enter the day of Nov. 2022 MVT to generate slim MVT data files (from 16 to 18): ');
end
%========================================================================
% Parameters
%========================================================================
% Set processing parameters
processingOpts.laneIdentificationOpts.LANEWIDTH = 12;       % Lane width in feet
% Number of cells in the x direction used to identify driving line
processingOpts.laneIdentificationOpts.Nr_XCELLS = 200;      
% Minimum number of samples per cell used to identify driving line
processingOpts.laneIdentificationOpts.MINCELLSAMPLES = 20; 
% samples upper bound as a multiplier of lane width to remove outliers
processingOpts.laneIdentificationOpts.Y_UP_LIM = 5;        
% samples lower bound as a multiplier of lane width to remove outliers
processingOpts.laneIdentificationOpts.Y_LOW_LIM = 0.5;        
% Shift parameters for y position correction (west{east}bound) s.t. y_corrected =  
% Sw{e} * (y  - drivingLineShift_w{e}) + Cw{e};
processingOpts.laneIdentificationOpts.Se = 0.97;            
processingOpts.laneIdentificationOpts.Ce = 1;
processingOpts.laneIdentificationOpts.Sw = 0.98;          
processingOpts.laneIdentificationOpts.Cw = 1;
% Threshold to identify lane change as a multiplier of lane width
processingOpts.laneChangeClippingOpts.LANECHANGETHRSH = 0.5;        
% Threshold for maximum rate of change in the identified lane for which 
% trajectory is assumed to not change lane as a multiplier of lane width/s
processingOpts.laneChangeClippingOpts.MAXLANECHANGERATE = 0.1;      
% Minimum duration for a clipped trajectory in seconds
processingOpts.laneChangeClippingOpts.MINCLIPTIME = 0.5;            
% Buffer threshold to ensure full lane change as a multiplier of lane width
processingOpts.laneChangeClippingOpts.CHANGEBUFFERTHRESH = 0.2;     
%========================================================================
% Initilize
%========================================================================
[parentDirectory, ~, ~] = fileparts(pwd);
dataFolderPath = fullfile(parentDirectory,['Data_2022-11-' num2str(DAY_TO_PROCESS) '__I24_Base']);
dayAbbrvs = ["mon","tue","wed","thu","fri"];
dayAbbrv = dayAbbrvs(DAY_TO_PROCESS-13);
dataFiles = dir([dataFolderPath '\*_' char(dayAbbrv) '_0_*.json']);
if length(dataFiles) < 24
    error('I24 base files for the day: %d, Nov. 2022 are missing or incomplete.',DAY_TO_PROCESS)
end
%========================================================================
% Load GPS and road grade data
%========================================================================
% Load road grade map
gradeData = readmatrix([parentDirectory '\Models_Energy\Eastbound_grade_fit.csv']);
% 0.225miles is the distance between mill creek origin (MM58.675) and MM58.9
% the estimated origin of the road grade map
gradeDataStart = gradeData(:,2);
gradeDataEnd = gradeData(:,3);
gradeDataPoints = [gradeDataStart; gradeDataEnd(end)];
gradeDataSlope = gradeData(:,4);
gradeDataIntercept = gradeData(:,5);
% GPS Data is assumed to be processed. Run create_data_GPS.m to produce processed GPS files
fprintf('\nLoading and decoding AVs GPS data file ...'); tic
dataGPS = jsondecode(fileread(fullfile('GPS_Data',['CIRCLES_GPS_10Hz_2022-11-' num2str(DAY_TO_PROCESS) '.json'])));
fprintf('Done (%0.0fsec).\n',toc)
%========================================================================
% Process each I24 MOTION file 
%========================================================================
addpath([parentDirectory '\Models_Energy']);
for fileNr = 1:24
    % Load MOTION data file
    filenameLoad = fullfile(dataFolderPath,dataFiles(fileNr).name);
    fprintf('Loading and decoding MOTION data file, %d/24 ... ', ...
        fileNr); tic
    dataTemp = jsondecode(fileread(filenameLoad));
    fprintf('Done (%0.0fsec).\n',toc)
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
    for i=1:floor(n) %loop over all data(vehicle) segments
        veh = dataTemp(i);               % get each vehicle
        data(i).trajectory_id = veh.x_id; % id for i-th vehicle trajectory
        data(i).timestamp = veh.timestamp; % time array for i-th vehicle trajectory
        % numbers corresponding to vehicle types 0:sedan, 1:midSUV, 2:van, 3:pickup, 4:semi, 5:truck
        data(i).coarse_vehicle_class = veh.coarse_vehicle_class; 
        data(i).first_timestamp = veh.first_timestamp; % first timestamp in the trajectory fragment
        data(i).last_timestamp = veh.last_timestamp; % last timestamp in the trajectory fragment
        data(i).length = veh.length; % length of vehicle in feet
        data(i).width = veh.width; % width of vehicle in feet
        data(i).height = veh.height; % height of vehicle in feet
        %lane number of the vehicle: 1 leftmost lane(HOV) , 4 right most lane. 0 off the highway, 5 on/off ramp.
        data(i).lane_number = veh.lane; 
        % shift x position to v2 local coordinate system (origin: Mill Creek Bridge) and convert to meters.
        data(i).x_position_meters = 0.3048*(veh.x_position-309804.125); 
        data(i).y_position_corrected_meters = veh.y_position*0.3048; % convert to meters after correction and clipping
        data(i).starting_x = data(i).x_position_meters(1);
        data(i).ending_x = data(i).x_position_meters(end);
        x = veh.x_position;          % original x-position for i-th vehicle
        x = 0.3048*abs(x-x(1));      % change position from feet to meters
        t = veh.timestamp;           % original timestamp data for i-th vehicle
        t = t-t(1);                  % start from time zero for each vehicle
        % Calculate velocity (m/s)
        v = (x([2:end,end])-x([1,1:end-1]))./(t([2:end,end])-t([1,1:end-1])); 
        % Calculate acceleration central (m/s/s)
        a = (x(1:end-2)-2*x(2:end-1)+x(3:end))./((t(3:end)-t(1:end-2))/2).^2; 
        a = a([1,1:end,end]);
        data(i).total_distance_traversed_meters = abs(x(end)-x(1));
        data(i).speed_meters_per_second = v;
        data(i).acceleration_meters_per_second_per_second = a;
        % approximate road grade theta
        % shift = 6.01*1609.344 + 309804.5*0.3048; % approx distance between symernia and origin in nashville city 
        % (6.01mil is distance between Sym and mm60 and 309804.5 is distance
        % from mm to origin of I24motion system near Nashville city)
        % xNew = shift-0.3048*veh.x_position;     
        % position of vehicles wrt to symernia origin (end of Sam Ridely Pkwy on ramp)
        % xNew = x_position_meters_symernia;
        xNew = data(i).x_position_meters/1609.344 - 0.225;
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
        data(i).road_grade_radians = theta;
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
        
        data(i).energy_model = vehType;
        % Calculate total fuel
        integrate = @(t,v) dot( t(2:end)-t(1:end-1) , (v(1:end-1)+v(2:end))/2); % quadrature
        % Fuel consumption rate of drive on road
        eval(sprintf(['[fS,~,infeasible] = fuel_model_'...
            vehType '_simplified(v,a,theta,true);']))
        
        data(i).fuel_rate_grams_per_second = fS;
        data(i).percent_infeasibility = length(infeasible(infeasible>0))/length(infeasible)*100; % percentage of infeasible drive
        data(i).total_fuel_consumed_grams = integrate(t,fS); % total fuel
        data(i).total_fuel_consumed_gallons = (1/2839.0588)*data(i).total_fuel_consumed_grams;
        data(i).total_fuel_economy_mpg = ((data(i).total_distance_traversed_meters)/1609.34)./...
            (data(i).total_fuel_consumed_gallons);
        
        % Distance to av engaged US
        if ~isempty(distToAvsData(i).AvIdUSEng)
            data(i).upstream_engaged_av_id = distToAvsData(i).AvIdUSEng;
            data(i).distance_to_upstream_engaged_av_meters = -distToAvsData(i).distanceToAvUSEng;
        else
            data(i).upstream_engaged_av_id = [];
            data(i).distance_to_upstream_engaged_av_meters = [];
        end
        
        % Distance to av US
        if ~isempty(distToAvsData(i).AvIdUS)
            data(i).upstream_av_id = distToAvsData(i).AvIdUS;
            data(i).distance_to_upstream_av_meters = -distToAvsData(i).distanceToAvUS;
        else
            data(i).upstream_av_id = [];
            data(i).distance_to_upstream_av_meters = [];
        end
        
        % Distance to av engaged DS
        if ~isempty(distToAvsData(i).AvIdDSEng)
            data(i).downstream_engaged_av_id = distToAvsData(i).AvIdDSEng;
            data(i).distance_to_downstream_engaged_av_meters = distToAvsData(i).distanceToAvDSEng;
        else
            data(i).downstream_engaged_av_id = [];
            data(i).distance_to_downstream_engaged_av_meters = [];
        end
        
        % Distance to av DS
        if ~isempty(distToAvsData(i).AvIdDS)
            data(i).downstream_av_id = distToAvsData(i).AvIdDS;
            data(i).distance_to_downstream_av_meters = distToAvsData(i).distanceToAvDS;
        else
            data(i).downstream_av_id = [];
            data(i).distance_to_downstream_av_meters = [];
        end
        %Set number of decimal points to be saved for each processed variable
        dataFields = fields(data);
        
        for iField =2:length(dataFields)
            if isnumeric(data(i).(string(dataFields(iField))))
                data(i).(string(dataFields(iField))) = round(data(i).(string(dataFields(iField))),4,'decimals');
            end
        end
        
    end
    clear dataTemp
    pause(5)
    fprintf('Done (%0.0fsec).\n',toc)
    fprintf('Encoding and Writing json file ... '),tic
    fileStartT = (datetime(data(1).first_timestamp, 'convertfrom', 'posixtime', 'Format', 'HH:mm:ss.SSS','TimeZone' ,'America/Chicago'));
    fileStartT = datestr(fileStartT,'YYYY-mm-dd_HH-MM-SS');
    filenameSave = [parentDirectory '\Data_2022-11-' num2str(DAY_TO_PROCESS)...
        '__MVT_Slim\I-24MOTION_slim_',fileStartT];    
    jsonStr = jsonencode(data);
    clear data
    fid = fopen([filenameSave '.json'], 'w');
    fwrite(fid, jsonStr, 'char');
    fclose(fid);
    fprintf('Done (%0.0fsec).\n',toc)
    clear jsonStr
    toc
end
rmpath([parentDirectory '\Models_Energy']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Local Functions Definitions %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = init_data_struct(n)
data(1,n) = struct('trajectory_id', [],'timestamp',[],'x_position_meters',[],'y_position_corrected_meters',[],...
    'coarse_vehicle_class',[],'direction',[],'first_timestamp',[],'last_timestamp',[],...
    'starting_x',[],'ending_x',[],'length',[],'width',[],'height',[],...
    'total_distance_traversed_meters',[],'speed_meters_per_second',[],'acceleration_meters_per_second_per_second',[],...
    'road_grade_radians',[],'lane_number',[],'energy_model',[],'fuel_rate_grams_per_second',[],...
    'percent_infeasibility',[],'total_fuel_consumed_grams',[],'total_fuel_consumed_gallons',[],...
    'total_fuel_economy_mpg',[],'fuel_rate_flat_road_grams_per_second',[],'downstream_av_id',[],...
    'distance_to_downstream_av_meters',[],'downstream_engaged_av_id',[],...
    'distance_to_downstream_engaged_av_meters',[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dataLanes] = assign_lanes(baseFileName,laneIdentificationOpts)

% v1.2
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

if isa(baseFileName,'struct')
    data = baseFileName;
elseif ischar(baseFileName)
    fprintf('Loading data \n')
    filenameLoad = baseFileName;
    data = jsondecode(fileread(filenameLoad));
else
    fprintf('Wrong input \n')
    return
end


%%% define Driving Line
% sample data
N = length(data);
xWest = []; yWest = [];
for i=1:N
    veh = data(i);
    t = veh.timestamp;
    ts = t(2)-t(1);
    ss = ceil(1 / ts); %sample at 1 Hz
    if (t(end)-t(1))>5     %sample only trajectories > 5s , each second
        xWest = [xWest, veh.x_position(1 :ss: end)'];
        yWest = [yWest, veh.y_position(1 :ss: end)'];
    end
end

%%% Process data
LANEWIDTH = laneIdentificationOpts.LANEWIDTH;    %ft
Nr_XCELLS = laneIdentificationOpts.Nr_XCELLS;
MINCELLSAMPLES = laneIdentificationOpts.MINCELLSAMPLES;
% samples upper/lower bound as a multiplier of lane width
Y_UP_LIM = laneIdentificationOpts.Y_UP_LIM;
Y_LOW_LIM = laneIdentificationOpts.Y_LOW_LIM;
upBoundY = LANEWIDTH*Y_UP_LIM;
lowBoundY = LANEWIDTH*Y_LOW_LIM;

%remove outliers from the sampled data based on bounds
xWest = xWest(yWest<=upBoundY & yWest>=lowBoundY);
yWest = yWest(yWest<=upBoundY & yWest>=lowBoundY);

xCellsWest = linspace(min(xWest),max(xWest),Nr_XCELLS+1);

drivingLineWest = zeros(1,Nr_XCELLS);

for i = 1:Nr_XCELLS    
    yCellWest = yWest(xCellsWest(i)<xWest&xWest<=xCellsWest(i+1)) - LANEWIDTH/2; % points in this cell
    
    if length(yCellWest) > MINCELLSAMPLES
        p_j_w = cos(2*pi*yCellWest/LANEWIDTH); q_j_w = sin(2*pi*yCellWest/LANEWIDTH); %mapping shift to unit circle
        dlw_t = atan2(mean(q_j_w), mean(p_j_w))*LANEWIDTH/2/pi;   %find mean shift
        drivingLineWest(i) = dlw_t;
    else
        drivingLineWest(i) = 0;
    end
    
end

xCellsWest = xCellsWest(1:end-1);


%%% Lane assignment
%set shifting factors used for y-position correction
Se = laneIdentificationOpts.Se; Ce = laneIdentificationOpts.Ce; Sw = laneIdentificationOpts.Sw; Cw = laneIdentificationOpts.Cw;
dataLanes = struct([]);
for i = 1: N
    veh = data(i);
    y = veh.y_position;
    x = veh.x_position;
    laneTemp = zeros(size(y));
    lane = zeros(size(y));
    n = length(laneTemp);
    
    drivingLine = drivingLineWest;
    xCells = xCellsWest;
    
    % correct for wiggle in y
    
    dataLanes(i).y_corr = Sw * (y  - interp1(xCellsWest,drivingLineWest,x,'linear','extrap')) + Cw;
    
    laneTemp = ((abs(dataLanes(i).y_corr) - LANEWIDTH/2) ./LANEWIDTH) ;
    laneTemp = max(0,min(5,(laneTemp)));
    
    % median filter the lane
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
    dataLanes(i).lane = lane;
    
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

LANECHANGETHRSH = opts.LANECHANGETHRSH;         % Threshold to identify lane change as a multiplier of lane width
MAXLANECHANGERATE = opts.MAXLANECHANGERATE;     % Threshold for maximum rate of change in the identified lane for which trajectory is assumed to not change lane as a multiplier of lane width
MINCLIPTIME = opts.MINCLIPTIME;                 % Minimum duration for a clipped trajectory in seconds
CHANGEBUFFERTHRESH = opts.CHANGEBUFFERTHRESH;   % Buffer threshold to ensure full lane change as a multiplier of lane width


nrTrajs = length(data);
dataFields = fieldnames(data);
newData = {};
for i = 1:length(dataFields) %init
    newData(4*nrTrajs).(cell2mat(dataFields(i))) = [];
end
newData(4*nrTrajs).lane =[];

trajCounter =1;
for i = 1:nrTrajs
       
    veh = data(i);
    x = veh.x_position;
    t = veh.timestamp;
    y = dataLanes(i).y_corr;
    lane = dataLanes(i).lane;
    trajectoryLength = length(veh.timestamp);
    tempLane = round(lane(1));
    ii=1;
    clippedPart = -1; %clipped part number (starting from 0)
    
    while ii<trajectoryLength
        ii0 = find((abs(diff(lane)./diff(t))) < MAXLANECHANGERATE,1); %start from when lane is stable
        if isempty(ii0)
            break %if lane is not stable at all, delete full trajectory
        end
        
        lane = lane(ii0:end);
        tempLane = round(lane(1));
        trajectoryLength = length(lane);
        x = x(ii0:end);
        y = y(ii0:end);
        t = t(ii0:end);
        
        % find index of lane changing
        indChange = find(abs(lane-tempLane) > LANECHANGETHRSH,1);
        if isempty(indChange); indChange =  trajectoryLength; end
        % remove data before lane changing (up to stable lane)
        mirroredLane0 = lane(indChange:-1:1);
        absDiffMirroredLane0 = abs(diff(mirroredLane0)./diff(t(indChange:-1:1)));
        absDiffMirroredLane0 = [absDiffMirroredLane0;absDiffMirroredLane0(end)];
        ii_p2 = find(abs(mirroredLane0 - tempLane) < CHANGEBUFFERTHRESH & absDiffMirroredLane0 < MAXLANECHANGERATE,1);
        % log clipped trajectory
        if ~isempty(ii_p2) && (t(indChange-ii_p2+1)-t(1)) >= MINCLIPTIME % only MINCLIPTIME [s] minimum
            
            
            clippedPart = clippedPart+1;
            ii_p2 = indChange - ii_p2;
            
            newVeh = veh;
            newVeh.timestamp = t(1:ii_p2);
            newVeh.x_position = x(1:ii_p2);
            newVeh.y_position = y(1:ii_p2);
            newVeh.ending_x = x(ii_p2);
            newVeh.first_timestamp = t(1);
            newVeh.last_timestamp = t(ii_p2);
            newVeh.starting_x = x(1);
            newVeh.ending_x = x(ii_p2);
            newVeh.lane = tempLane;
            newVeh.x_id.x_oid = [newVeh.x_id.x_oid '-' num2str(clippedPart)];
            newData(trajCounter) = newVeh;
            trajCounter = trajCounter+1;
        end
        
        lane = lane(indChange:end);
        x = x(indChange:end);
        y = y(indChange:end);
        t = t(indChange:end);
        tempLane = round(lane(1));
        
        % if clipping happen, remove starting part of the remaining trajectory
        % until within (CHANGEBUFFERTHRESH * lane width) of new lane
        if indChange<trajectoryLength
            ii_n2 = find(abs(lane-tempLane)<CHANGEBUFFERTHRESH,1);
            if isempty(ii_n2)
                ii_n2 =  length(lane);
            end
            ii_notfull = find(abs(lane-tempLane)>1-CHANGEBUFFERTHRESH,1);
            if ~isempty(ii_notfull)&& ii_notfull< ii_n2  %false lane change: clip but don't start new trajectory
                ii_n2 =  ii_notfull;
                tempLane = round(lane(ii_n2));
            end
            ii = ii_n2+1;
            lane = lane(ii:end);
            x = x(ii:end);
            y = y(ii:end);
            t = t(ii:end);
            ii = 1;
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

n = length(dataTemp);
DistToAvs(1,n) = struct('distanceToAvUS',[],'AvIdUS',[],'distanceToAvUSEng',[],'AvIdUSEng',[],...
    'distanceToAvDS',[],'AvIdDS',[],'distanceToAvDSEng',[],'AvIdDSEng',[]);
distKeys = ["distanceToAvUS","distanceToAvDS"];
distEngKeys = ["distanceToAvUSEng","distanceToAvDSEng"];
idsKeys = ["AvIdUS","AvIdDS"];
idsEngKeys = ["AvIdUSEng","AvIdDSEng"];
distDir = ["upstream","downstream"];
distDirSign = [-1 ,1];

num_trajs = length(dataTemp);
for i = 1:num_trajs
    
    veh = dataTemp(i);
    tempT = veh.timestamp;
    tempX = veh.x_position;
    tempX = 0.3048*(tempX-309804.125);
    tempLane = veh.lane;
    tempDir = veh.direction;
    avsOnRoad = find([dataGPS.assigned_lane]==tempLane& [dataGPS.direction] == tempDir & tempT(1)< [dataGPS.last_timestamp] & tempT(end)>[dataGPS.first_timestamp]);
    
    if ~isempty(avsOnRoad)
        
        for keyDir = 1:2
            sn = distDirSign(keyDir);
            key = distDir(keyDir);
            vehDistToAvs.(distDir(keyDir)) = distDirSign(keyDir)*inf*(ones(length(avsOnRoad),length(tempT)));
            vehDistToAvsEngaged.(distDir(keyDir)) = distDirSign(keyDir)*inf*(ones(length(avsOnRoad),length(tempT)));
            tempDistX = [];
            for ii = 1:length(avsOnRoad)
                
                avInd = avsOnRoad(ii);
                projAvsX = interp1(dataGPS(avInd).timestamp,dataGPS(avInd).x_position,tempT) ;
                tempDistX = (projAvsX - tempX )* tempDir;
                
                if ~isempty(dataGPS(avInd).control_car)
                    temp_controller_active = min(1,max(0,round(interp1(dataGPS(avInd).timestamp,double(dataGPS(avInd).control_car),tempT)))) ;
                else
                    temp_controller_active = min(1,max(0,round(interp1(dataGPS(avInd).timestamp,double(dataGPS(avInd).controller_engaged),tempT)))) ;
                end
                
                vehDistToAvs.(key)(ii,:) =  tempDistX.';
                vehDistToAvs.(key)(ii,isnan(tempDistX)) = sn*inf;
                
                vehDistToAvsEngaged.(key)(ii,:) =  tempDistX.';
                vehDistToAvsEngaged.(key)(ii,isnan(tempDistX)) = sn*inf;
                vehDistToAvsEngaged.(key)(ii,temp_controller_active==0) = sn*inf;
            end
        end
        
        for keyDir = 1:2
            sn = distDirSign(keyDir);
            key = distDir(keyDir);
            vehDistToAvs.(key)(sn*vehDistToAvs.(key)<0) = sn*inf;
            vehDistToAvsEngaged.(key)(sn*vehDistToAvsEngaged.(key)<0) = sn*inf;
            
            [nearestDist.(key),Ind.(key)] = min(sn*vehDistToAvs.(key),[],1);
            [nearestDistEngaged.(key),indEngaged.(key)] = min(sn*vehDistToAvsEngaged.(key),[],1);
            
            DistToAvs(i).(distKeys(keyDir)) = nearestDist.(key);
            DistToAvs(i).(idsKeys(keyDir)) = [dataGPS(avsOnRoad(Ind.(key))).av_id];
            infFlag = isinf(nearestDist.(key));
            if sum(infFlag) == length(infFlag)
                DistToAvs(i).(distKeys(keyDir)) = [];
                DistToAvs(i).(idsKeys(keyDir)) = [];
            else
                DistToAvs(i).(distKeys(keyDir))(infFlag) = nan;
                DistToAvs(i).(idsKeys(keyDir))(infFlag) = nan;
            end
            
            DistToAvs(i).(distEngKeys(keyDir)) = nearestDistEngaged.(key);
            DistToAvs(i).(idsEngKeys(keyDir)) = [dataGPS(avsOnRoad(indEngaged.(key))).av_id];
            infFlag = isinf(nearestDistEngaged.(key));
            if sum(infFlag) == length(infFlag)
                DistToAvs(i).(distEngKeys(keyDir)) = [];
                DistToAvs(i).(idsEngKeys(keyDir)) = [];
            else
                DistToAvs(i).(distEngKeys(keyDir))(infFlag) = nan;
                DistToAvs(i).(idsEngKeys(keyDir))(infFlag) = nan;
            end
            
        end
        
        
    end
    
end
end


