% Script to process AV GPS data from the MVT and produce JSON file for a
% test day in the same directory.
% (C) 2025/04/15 by Sulaiman Almatrudi

clearvars -except DAY_TO_PROCESS
if ~exist('DAY_TO_PROCESS', 'var')
    global DAY_TO_PROCESS
    DAY_TO_PROCESS = input('Enter the day of Nov. 2022 MVT to generate AVs GPS data files (from 16 to 18): ');
end
fprintf('Starting the assembly of AV GPS data for %dth Nov. 2022\n', DAY_TO_PROCESS);

%========================================================================
% Parameters
%========================================================================
%%%% Threshold choices for matching AVs GPS to I24 MOTION trajectories 
maxMatchDist = 6; %[m] maximum distance between AVs GPS trajectories and potential matching I24 trajectories
maxMatchSpdDiff = 2; % [m/s] maximum speed difference 
maxMatchLaneDiff = 0.5; %[lane] % maximum estimated lane difference 
minMatchTime = 3; % [s] minimum matching time 
%========================================================================
% Load, parse and preprocess AVs GPS data
%========================================================================
[parentDirectory, ~, ~] = fileparts(pwd);
gpsFolderPath = fullfile(parentDirectory,'\Data_GPS') ;
if ~isfolder(gpsFolderPath)
    error('Folder %s does not exist.\n',gpsFolderPath)
end
% Parse av files 
fprintf('Parsing AV data... ') ; tic
dataGPS0 = parse_gps_data(gpsFolderPath);
nAVruns = length(dataGPS0);
fprintf('Done (%0.0fsec).\n',toc)

% Pre-process av files
fprintf('Pre-processing AV data... ');tic
dataGPS0 = preproc_gps(dataGPS0,gpsFolderPath);
fprintf('Done (%0.0fsec).\n',toc)

% Find AV activity limits
dataTLimits = datetime(2022,11,DAY_TO_PROCESS,[6 10],0,0,'TimeZone','America/Chicago');
dataTLimits = convertTo(dataTLimits,'epochtime');
minAVStart = min([dataGPS0([dataGPS0.starting_time]>dataTLimits(1)).starting_time]);
minFileNr = max(1,floor((minAVStart-dataTLimits(1))/600));
maxAVStart = max([dataGPS0([dataGPS0.ending_time]<dataTLimits(2)).ending_time]);
maxFileNr = min(24,floor((maxAVStart-dataTLimits(1))/600)+1);

% Find I24 MOTION data files in folder
dataFolderPath = fullfile(parentDirectory,['Data_2022-11-' num2str(DAY_TO_PROCESS) ...
    '__I24_Base']) ;
if ~isfolder(dataFolderPath)
    error('Folder %s does not exist.\n',dataFolderPath)
end
dayAbbrvs = ["wed","thu","fri"];
dayAbbrv = dayAbbrvs(DAY_TO_PROCESS-15);
dataFiles = dir([dataFolderPath '\*_' char(dayAbbrv) '_0_*.json']);
if length(dataFiles) < 24
    error('I24 base files for the day: %d, Nov. 2022 are missing or incomplete.',DAY_TO_PROCESS)
end
dateData = ['2022-11-' num2str(DAY_TO_PROCESS)] ;

%========================================================================
% Find I24 trajectories matching AV gps trajectories
%========================================================================
d = 1;
matchedSegments(10*length(dataGPS0),1).av_trj_i = [];
matchedSegments(10*length(dataGPS0),1).x_diff =   [];
matchedSegments(10*length(dataGPS0),1).timestamp = [];
for fileNr = minFileNr:maxFileNr
    % Load MOTION data file
    filenameLoad = fullfile(dataFolderPath,dataFiles(fileNr).name);
    fprintf('Loading and decoding MOTION data file, %d/%d ... ', ...
        fileNr-minFileNr+1,maxFileNr-minFileNr+1); tic
    dataTemp = jsondecode(fileread(filenameLoad));
    fprintf('Done (%0.0fsec).\n',toc)
    % Find AV runs concurrent to MOTION data trajectories in file
    dataAVsInFile = dataGPS0( [dataGPS0.starting_time] <= max([dataTemp.first_timestamp]) ...
        & [dataGPS0.ending_time] >= min([dataTemp.first_timestamp]) );
    if isempty(dataAVsInFile)
        fprintf('Possible AV trajctories not found in file \n') ; tic
        continue;
    else
        fprintf('Processing file for possible AV matches... ') ; tic
    end
    % Delete unnecessary fields in MOTION data
    dataTemp = rmfield(dataTemp,{'x_id','flags','length','width','height',...
        'merged_ids','road_segment_ids','coarse_vehicle_class','fine_vehicle_class'...
        ,'fragment_ids','compute_node_id','local_fragment_id','configuration_id',...
        'x_score','y_score'});
    % Identify and assign lanes to MOTION data
    dataLanes = assign_lanes(dataTemp);
    % Find matched segments from I24 MOTION to each AV
    for ii =1:length(dataAVsInFile)
        avVeh =   dataAVsInFile(ii);
        avX = avVeh.x_position;
        avLane = avVeh.assigned_lane;
        avT = avVeh.timestamp;
        % Find MOTION trajectories that could possibly match AV
        possTrajsMatch = find([dataTemp.first_timestamp] <= avT(end) ...
            & [dataTemp.last_timestamp] >= avT(1) & [dataTemp.direction] == avVeh.direction);
        for ii2 = 1:length(possTrajsMatch)
            traj = dataTemp(possTrajsMatch(ii2));
            trajLane =  dataLanes(possTrajsMatch(ii2)).lane;
            trajX = (traj.x_position-309804.125);
            trajT = traj.timestamp;
            trajV = diff(trajX)./diff(trajT);
            trajV = [trajV(1);trajV]*0.3048;
            % Calculate distance to AV
            distToAV = (interp1(avT,avX,traj.timestamp) - trajX) * 0.3048;
            avgDistToAVX = distToAV;
            avgDistToAVX = avgDistToAVX(~isnan(avgDistToAVX));
            medianDistToAVX = median(abs(avgDistToAVX));
            avgDistToAVX = mean(abs(avgDistToAVX));
            % Calculate difference in lane to AV
            diffToAVLane = avLane-trajLane;
            avgdiffToAVLane = diffToAVLane;
            avgdiffToAVLane = avgdiffToAVLane(~isnan(avgdiffToAVLane));
            mediandiffToAVLane = median(abs(avgdiffToAVLane));
            avgdiffToAVLane = mean(abs(avgdiffToAVLane));
            % Ignore trajectory if it has an average absolute x position difference
            % of > 200 m or average abs. lane difference of >= 2 lanes, for
            % faster processing
            if avgDistToAVX <= 200 && avgdiffToAVLane < 2
                vDiff = diff(distToAV)./diff(trajT);
                vDiff = [vDiff(1);vDiff];
                vDiffS = abs(smoothdata(vDiff, 'gaussian', 3, 'SamplePoints', trajT));
                timeMatching = 0;
                for k = 2:length(trajT)
                    if abs(distToAV(k)) <= maxMatchDist ...
                            && abs(diffToAVLane(k)) <= maxMatchLaneDiff ...
                            && vDiffS(k) <= maxMatchSpdDiff
                        if timeMatching == 0
                            matchStartInd = k;
                        end
                        timeMatching = timeMatching + trajT(k) - trajT(k-1);
                    else
                        if timeMatching >= minMatchTime
                            matchedSegments(d).av_trj_i = avVeh.index;
                            matchedSegments(d).x_diff =   distToAV(matchStartInd:k);
                            matchedSegments(d).timestamp = trajT(matchStartInd:k);
                            d = d+1;
                        end
                        matchStartInd = [];
                        timeMatching = 0;
                    end
                end
                if timeMatching >= minMatchTime
                    matchedSegments(d).av_trj_i = avVeh.index;
                    matchedSegments(d).x_diff =   distToAV(matchStartInd:k);
                    matchedSegments(d).timestamp = trajT(matchStartInd:k);
                    d = d+1;
                end
            end

        end
    end
fprintf('Done (%0.0fsec).\n',toc)
end
matchedSegments(d:end) = [];
% Calculate median offset between GPS and I24 matched signals
xd_match = [];
xd_match(length(dataGPS0),1).xd = [];
xd_match(length(dataGPS0),1).median_xd = [];
for av_ri =1 :length(dataGPS0)
    A = find([matchedSegments.av_trj_i]==av_ri);
    xd_match(av_ri).xd = vertcat( matchedSegments(A).x_diff);
    if ~isempty(xd_match(av_ri).xd)
        xd_match(av_ri).median_xd = median(xd_match(av_ri).xd(~isnan(xd_match(av_ri).xd)));
    else
        xd_match(av_ri).median_xd = 0;
    end
end
%%
%========================================================================
% Add controller status and clean data structure
%========================================================================
fprintf('Adding controller status and preparing output data... ');tic
avPingsData = readtable(fullfile(gpsFolderPath, ['veh_ping_202211' num2str(DAY_TO_PROCESS) '.csv']));
avVINs = readtable([gpsFolderPath '\veh_vins.csv']);
avConnectionStatus  = get_connection_status(avPingsData,avVINs);
nAVs = length(dataGPS0);
clear dataGPS
dataGPS(nAVs,1).av_id = [];
dataGPS(nAVs,1).assigned_lane = [];
dataGPS(nAVs,1).direction = [];
dataGPS(nAVs,1).timestamp = [];
dataGPS(nAVs,1).latitude = [];
dataGPS(nAVs,1).longitude = [];
dataGPS(nAVs,1).x_position = [];
dataGPS(nAVs,1).y_position = [];
dataGPS(nAVs,1).controller_engaged = [];
dataGPS(nAVs,1).speed = [];
dataGPS(nAVs,1).is_server_connected = [];
dataGPS(nAVs,1).first_timestamp = [];
dataGPS(nAVs,1).last_timestamp = [];
dataGPS(nAVs,1).control_car = [];
dataGPS(nAVs,1).control_last30 = [];
for avInd =1:nAVs
    avVeh = dataGPS0(avInd);
    k1 = find(avConnectionStatus(avVeh.vin).timestamp >= dataGPS0(avInd).timestamp(1),1);
    isConnected = (0*dataGPS(avInd).timestamp);
    if ~isempty(k1)
        runTimeConnected = avConnectionStatus(avVeh.vin).timestamp(k1:end);
        for k2 = 1:length(dataGPS(avInd).timestamp)
            tempInd = find(runTimeConnected>avVeh.timestamp(k2),1);
            if isempty(tempInd); tempInd = length(runTimeConnected)+1;end
            isConnected(k2) = single((runTimeConnected(tempInd-1) ...
                - avVeh.timestamp(k2)) >= -2); %service is connected if a ping is received in the past 2 s
        end

    end

    dataGPS(avInd).av_id = avVeh.vin;
    dataGPS(avInd).assigned_lane = avVeh.assigned_lane;
    dataGPS(avInd).direction = avVeh.direction;

    dataGPS(avInd).timestamp = round(avVeh.timestamp,6,'decimals');
    dataGPS(avInd).latitude = round(avVeh.Lat,6,'decimals');
    dataGPS(avInd).longitude = round(avVeh.Long,6,'decimals');
    dataGPS(avInd).x_position = round(0.3048 * avVeh.x_position ...
        - xd_match(avInd).median_xd,6,'decimals'); %remove bias detected by  matching to I24MOTION data
    dataGPS(avInd).y_position = round(0.3048 * avVeh.y_position,6,'decimals');

    dataGPS(avInd).controller_engaged = avVeh.control_active;
    dataGPS(avInd).speed = avVeh.can_speed;

    dataGPS(avInd).is_server_connected = isConnected;
    dataGPS(avInd).first_timestamp = avVeh.timestamp(1);
    dataGPS(avInd).last_timestamp = avVeh.timestamp(end);

    dataGPS(avInd) = get_control_car_status(dataGPS(avInd));
end

% Clip trajectories at testbed limits and prepare data for writitng
XLIMS = [-400 2.54e4] *0.3048; % testbed limits 
gpsDataFields = string(fieldnames(dataGPS));
for avInd = 1:length(dataGPS)
    inTestBed = dataGPS(avInd).x_position > XLIMS(1) & dataGPS(avInd).x_position < XLIMS(2);
    if sum(inTestBed)>0
        for fieldNum = 1:length(gpsDataFields)
            if length(dataGPS(avInd).(gpsDataFields(fieldNum))) > 1
                dataGPS(avInd).(gpsDataFields(fieldNum)) = ...
                    dataGPS(avInd).(gpsDataFields(fieldNum))(inTestBed);
            end
        end
        dataGPS(avInd).first_timestamp = dataGPS(avInd).timestamp(1);
        dataGPS(avInd).last_timestamp = dataGPS(avInd).timestamp(end);
    end

end
fprintf('Done (%0.0fsec).\n',toc)

%========================================================================
% Encode and write to file
%========================================================================
filenameWrite = ['CIRCLES_GPS_10Hz_2022-11-' num2str(DAY_TO_PROCESS) '.json'];
fprintf('Encoding and writing file %s to file... ',filenameWrite);tic
jsonStr = jsonencode(dataGPS);
fid = fopen(fullfile(gpsFolderPath, filenameWrite), 'w');
fwrite(fid, jsonStr, 'char');
fclose(fid);
clear jsonStr
fprintf('Done (%0.0fsec).\n',toc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Local Functions Definitions %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dataGPS0 = parse_gps_data(gpsDataFolder)
global DAY_TO_PROCESS
%========================================================================
% Parameters
%========================================================================
maxRcsx = 2.91e4; % [ft] maximum x position
minRcsx = -3000; % [ft] minimum x position
maxRcsy = 150; % [ft] maximum y position
maxGpsInactive = 0.01; % maximum percentage of run length with inactive GPS signal
minRunLength = 1.2; % [km] minimum run length along x axis
minRunTime = 60; % [s] minimum run time

%========================================================================
% Initialize
%========================================================================
dayLimits = datetime(2022,11,DAY_TO_PROCESS,[3 18],0,0,'TimeZone','America/Chicago');
dayLimits = convertTo(dayLimits,'epochtime');
dataGPS0 = {};
i = 1;
vins = [1:103];
vehsFiles = dir([gpsDataFolder '\circles_v2_1_car*.csv']);
opts1 = detectImportOptions([gpsDataFolder '\' vehsFiles(1).name]);
opts1.VariableTypes(13) = {'char'};
nrFiles = length(vehsFiles);
vinsWithFiles = zeros(1,nrFiles);
for iFile = 1:nrFiles
    tempVin = vehsFiles(iFile).name(end-6:end-4);
    tempVin =  tempVin(tempVin>='0' & tempVin<='9');
    vinsWithFiles(iFile) = str2double(tempVin);
end

%========================================================================
% Load and parse GPS data for each AV for the specified day
%========================================================================
for vin = 1:length(vins)
    runNum = 0;
    % Load data file
    if ismember(vin,vinsWithFiles)
        vehTable = readtable([gpsDataFolder '\circles_v2_1_car' num2str(vin) '.csv'],opts1);
    else
        fprintf(['data file for car ' num2str(vin) ' not found \n'])
        continue
    end
    % Filter data for the day
    vehTable = vehTable(vehTable.Systime > dayLimits(1) & vehTable.Systime < dayLimits(2),:);
    if isempty(vehTable)
        fprintf(['no data found for car ' num2str(vin) ' on day ' num2str(DAY_TO_PROCESS) ...
            ', car assumed inactive for the day \n'])
        continue
    end
    % Parse runs of the AV in both directions
    diffY = diff(vehTable.rcs_y);
    diffY = [diffY(1);diffY];
    runStart = find(sign(diffY) .* sign( vehTable.rcs_y) < 0 ...
        & abs(vehTable.rcs_y) < maxRcsy & vehTable.rcs_x > minRcsx ...
        & vehTable.rcs_x < maxRcsx ,1);
    vehTable = vehTable(runStart:end,:);
    runEnd = find(abs(vehTable.rcs_y) > maxRcsy | vehTable.rcs_x < minRcsx ...
        | vehTable.rcs_x > maxRcsx,1);
    runTime =  vehTable(1:runEnd,:).Systime;
    runX =  vehTable(1:runEnd,:).rcs_x;
    status = cell2mat(vehTable(1:runEnd,:).Status);
    while size(vehTable,1)>1
        if runTime(end)-runTime(1) > minRunTime && ...
                abs(runX(end)-runX(1)) > minRunLength * 3281 && ...
                sum(double(status)>80) < maxGpsInactive*length(runX)

            runNum = runNum+1;
            runDirection = sign(runX(end) - runX(1));
            runY = - vehTable(1:runEnd,:).rcs_y;
            dataGPS0(i).vin = vin;
            dataGPS0(i).run_num = runNum;
            dataGPS0(i).timestamp = runTime;
            dataGPS0(i).direction = runDirection;
            dataGPS0(i).x_position = runX;
            dataGPS0(i).y_position = runY;
            dataGPS0(i).starting_time = runTime(1);
            dataGPS0(i).ending_time = runTime(end);
            dataGPS0(i).Long = vehTable(1:runEnd,:).Long;
            dataGPS0(i).Lat = vehTable(1:runEnd,:).Lat;
            dataGPS0(i).state_x = vehTable(1:runEnd,:).state_x;
            dataGPS0(i).state_y = vehTable(1:runEnd,:).state_y;
            dataGPS0(i).can_speed = vehTable(1:runEnd,:).can_speed;
            cont_st = string(vehTable(1:runEnd,:).control_active);
            dataGPS0(i).control_active = (cont_st == "True");
            i = i+1;
        end
        if runEnd < size(vehTable,1)
            vehTable = vehTable(runEnd+1:end,:);
            diffY = diff(vehTable.rcs_y);
            diffY = [diffY(1);diffY];
            runStart = find(sign(diffY) .* sign( vehTable.rcs_y) < 0 ...
                & abs(vehTable.rcs_y) < maxRcsy & vehTable.rcs_x > minRcsx & vehTable.rcs_x < maxRcsx ,1);
            vehTable = vehTable(runStart:end,:);
            runEnd = find(abs(vehTable.rcs_y) > maxRcsy | vehTable.rcs_x < minRcsx ...
                | vehTable.rcs_x > maxRcsx,1);
            if isempty(runEnd); runEnd =size(vehTable,1);end
            runTime =  vehTable(1:runEnd,:).Systime;
            runX = vehTable(1:runEnd,:).rcs_x;
            status = cell2mat(vehTable(1:runEnd,:).Status);
        else
            vehTable = [];
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dataGPS0 = preproc_gps(dataGPS0,gpsDataFolder)
nAVs = length(dataGPS0);
lane_width = 12;
avVINs = readtable([gpsDataFolder '\veh_vins.csv']);
for runIdx = 1:nAVs
    t10Hz = floor(dataGPS0(runIdx).timestamp(1)*10)/10:.1:ceil(dataGPS0(runIdx).timestamp(end)*10)/10;
    dataGPS0(runIdx).index = runIdx;
    dataGPS0(runIdx).x_position = dataGPS0(runIdx).x_position - dataGPS0(runIdx).direction * 7;
    dataGPS0(runIdx).x_position = sample_10hz(dataGPS0(runIdx).timestamp, dataGPS0(runIdx).x_position, t10Hz);
    dataGPS0(runIdx).y_position = sample_10hz(dataGPS0(runIdx).timestamp, dataGPS0(runIdx).y_position, t10Hz);
    dataGPS0(runIdx).Lat = sample_10hz(dataGPS0(runIdx).timestamp, dataGPS0(runIdx).Lat, t10Hz);
    dataGPS0(runIdx).Long = sample_10hz(dataGPS0(runIdx).timestamp, dataGPS0(runIdx).Long, t10Hz);
    dataGPS0(runIdx).state_x = sample_10hz(dataGPS0(runIdx).timestamp, dataGPS0(runIdx).state_x, t10Hz);
    dataGPS0(runIdx).state_y = sample_10hz(dataGPS0(runIdx).timestamp, dataGPS0(runIdx).state_y, t10Hz);
    dataGPS0(runIdx).can_speed = sample_10hz(dataGPS0(runIdx).timestamp, dataGPS0(runIdx).can_speed, t10Hz);
    dataGPS0(runIdx).control_active = logical(round(...
        sample_10hz(dataGPS0(runIdx).timestamp, single(dataGPS0(runIdx).control_active), t10Hz)));
    dataGPS0(runIdx).timestamp = t10Hz;
    dataGPS0(runIdx).assigned_lane = avVINs(dataGPS0(runIdx).vin,:).lane_num;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newValue = sample_10hz(orgT, orgValue, newT)
B = orgT;
idx = [true;diff(B)>0.09&diff(B)<0.11];
B = B(idx);
A = orgValue(idx);
newValue = interp1(B,A,newT,'linear','extrap');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dataLanes] = assign_lanes(baseFileName)

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

%%%% Processing optinons for lane identification in the v2 data %%%%
laneIdentificationOpts.LANEWIDTH = 12;       % Lane width in feet
laneIdentificationOpts.Nr_XCELLS = 200;      % Number of cells in the x direction used to identify driving line
laneIdentificationOpts.MINCELLSAMPLES = 20;  % Minimum number of samples per cell used to identify driving line
laneIdentificationOpts.Y_UP_LIM = 5;         % samples upper bound as a multiplier of lane width to remove outliers
laneIdentificationOpts.Y_LOW_LIM = 0.5;      % samples lower bound as a multiplier of lane width to remove outliers
laneIdentificationOpts.Se = 0.97;            % Shift parameters for y position correction (westbound) s.t. y_corrected =  Sw * (y  - drivingLineShift_w) + Cw;
laneIdentificationOpts.Ce = 1;
laneIdentificationOpts.Sw = 0.98;            % Shift parameters for y position correction (eastbound) s.t. y_corrected =  Se * (y  - drivingLineShift_e) + Ce;
laneIdentificationOpts.Cw = 1;


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


%%% define Driving Line
% sample data
N = length(data);
xWest = []; yWest = [];
xEast = []; yEast = [];
for i=1:N
    veh = data(i);
    t = veh.timestamp;
    ts = t(2)-t(1);
    ss = ceil(1 / ts); %sample at 1 Hz
    direction = veh.direction;
    if (t(end)-t(1))>5     %sample only trajectories > 5s , each second
        if direction == -1   %west bound
            xWest = [xWest, veh.x_position(1 :ss: end)'];
            yWest = [yWest, veh.y_position(1 :ss: end)'];
        else                %east bound
            xEast = [xEast, veh.x_position(1 :ss: end)'];
            yEast = [yEast, veh.y_position(1 :ss: end)'];
        end

    end
end

%%% Process data

yEast = -yEast;
LANEWIDTH = laneIdentificationOpts.LANEWIDTH;    %ft
Nr_XCELLS = laneIdentificationOpts.Nr_XCELLS;
MINCELLSAMPLES = laneIdentificationOpts.MINCELLSAMPLES;
% samples upper/lower bound as a multiplier of lane width
Y_UP_LIM = laneIdentificationOpts.Y_UP_LIM;
Y_LOW_LIM = laneIdentificationOpts.Y_LOW_LIM;
upBoundY = LANEWIDTH*Y_UP_LIM;
lowBoundY = LANEWIDTH*Y_LOW_LIM;

%remove outliers from the sampled data based on bounds
xEast = xEast(yEast<=upBoundY & yEast>=lowBoundY);
yEast = yEast(yEast<=upBoundY & yEast>=lowBoundY);
xWest = xWest(yWest<=upBoundY & yWest>=lowBoundY);
yWest = yWest(yWest<=upBoundY & yWest>=lowBoundY);

xCellsEast = linspace(min(xEast),max(xEast),Nr_XCELLS+1);
xCellsWest = linspace(min(xWest),max(xWest),Nr_XCELLS+1);

drivingLineWest = zeros(1,Nr_XCELLS);
drivingLineEast = zeros(1,Nr_XCELLS);

for i = 1:Nr_XCELLS

    yCellWest = yWest(xCellsWest(i)<xWest&xWest<=xCellsWest(i+1)) - LANEWIDTH/2; % points in this cell
    yCellEast = yEast(xCellsEast(i)<xEast&xEast<=xCellsEast(i+1)) - LANEWIDTH/2;


    if length(yCellWest) > MINCELLSAMPLES

        p_j_w = cos(2*pi*yCellWest/LANEWIDTH); q_j_w = sin(2*pi*yCellWest/LANEWIDTH); %mapping shift to unit circle
        dlw_t = atan2(mean(q_j_w), mean(p_j_w))*LANEWIDTH/2/pi;   %find mean shift
        drivingLineWest(i) = dlw_t;
    else

        %             warning(['not enough samples in cell ' num2str(i) '. 0 shift assumed in west driving line'])
        drivingLineWest(i) = 0;
    end

    if length(yCellEast) > MINCELLSAMPLES

        p_j_e = cos(2*pi*yCellEast/LANEWIDTH); q_j_e = sin(2*pi*yCellEast/LANEWIDTH); %mapping shift to unit circle
        dle_t = atan2(mean(q_j_e), mean(p_j_e))*LANEWIDTH/2/pi; %find mean shift
        drivingLineEast(i) = dle_t;
    else
        %             warning(['not enough samples in cell ' num2str(i) '. 0 shift assumed in east driving line'])
        drivingLineEast(i) = 0;
    end

end

xCellsWest = xCellsWest(1:end-1);
xCellsEast = xCellsEast(1:end-1);

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

    if veh.direction >0        %east
        y = -y;
        drivingLine = drivingLineEast;
        xCells = xCellsEast;
    else                        %west
        drivingLine = drivingLineWest;
        xCells = xCellsWest;
    end

    % correct for wiggle in y
    if veh.direction >0
        dataLanes(i).y_corr = - (Se * (y  - interp1(xCellsEast,drivingLineEast,x,'linear','extrap')) + Ce);
    else
        dataLanes(i).y_corr = Sw * (y  - interp1(xCellsWest,drivingLineWest,x,'linear','extrap')) + Cw;
    end

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
function avConnectionStatus  = get_connection_status(avPingsData,avVINs)
avPingsData.av_id = 0*avPingsData.gpstime;

for i = 1:size(avVINs,1)
    VIN = avVINs(i,:).vin;
    vinDataPoints = strcmp(string(avPingsData.vin),string(VIN));
    avPingsData(vinDataPoints,:).av_id = avVINs(i,:).veh_id*ones(sum(vinDataPoints),1);
end

avConnectionStatus = {};
avConnectionStatus.vin = [];
avConnectionStatus.timestamp = [];
avConnectionStatus.status = [];

for i =1:103
    avConnectionStatus(i).vin = i;
    rawStatus = avPingsData(avPingsData.av_id==i,:).acc_status;
    avConnectionStatus(i).status = ~isnan( rawStatus )  ;
    avConnectionStatus(i).timestamp = avPingsData(avPingsData.av_id==i,:).gpstime/1e3;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function veh1  = get_control_car_status(veh)
    veh1 = veh;
    if veh.speed == 0 % calculate speed from position if can speed is unavailable
        x = veh.x_position(:);
        vehSpeed0 = (x(2:end) - x(1:end-1))./(x(2:end) - x(1:end-1));
        vehSpeed0 = round(abs([vehSpeed0;vehSpeed0(end)])*2)/2;
        paddedSpeed = [zeros(5,1);vehSpeed0;zeros(5,1)];
        vehSpeed = movmedian(paddedSpeed,10);
        vehSpeed = vehSpeed(6:end-5);
    end
    vehSpeed = veh.speed(:);
    tempControlCarStatus = -1*ones(size(veh.controller_engaged));
    controlActive = veh.controller_engaged(:);
    shiftedIndControlActive = [controlActive(2:end);controlActive(end)];
    indActiveMoving = controlActive & vehSpeed~=0;
    tempControlCarStatus(indActiveMoving) = 1;
    indNotActiveMoving = ~controlActive & vehSpeed~=0;
    tempControlCarStatus(indNotActiveMoving) = 0;
    indActiveStopped = controlActive & vehSpeed==0;
    indActiveStopped = find(indActiveStopped);

    if ~isempty(indActiveStopped)
        for k = 1:length(indActiveStopped)
            tempIndActiveStopped = indActiveStopped(k);
            tempIndEndActiveStopped = find(controlActive);
            tempIndEndActiveStopped = tempIndEndActiveStopped(find(tempIndEndActiveStopped > tempIndActiveStopped,1))-1;
            if isempty(tempIndEndActiveStopped) ; tempIndEndActiveStopped = length(veh.timestamp);end

            if sum(vehSpeed(tempIndActiveStopped:tempIndEndActiveStopped)>0) < min(10,tempIndEndActiveStopped-tempIndActiveStopped+2)
                tempControlCarStatus(tempIndActiveStopped:tempIndEndActiveStopped) = 1;

            else
                tempControlCarStatus(tempIndActiveStopped:tempIndEndActiveStopped) = 0;
            end
        end
    end

    tempControlCarStatus(indActiveStopped) = 1;
    indNotActiveStopped = find(~controlActive(1:end-1) & vehSpeed(1:end-1)~=0 & vehSpeed(2:end)==0);
    if ~isempty(indNotActiveStopped)
        for k = 1:length(indNotActiveStopped)
            tempIndNotActiveStopped = indNotActiveStopped(k);
            tempIndEndNotActiveStopped = find(controlActive | vehSpeed~=0);
            tempIndEndNotActiveStopped = tempIndEndNotActiveStopped(find(tempIndEndNotActiveStopped > tempIndNotActiveStopped,1))-1;
            if isempty(tempIndEndNotActiveStopped) ; tempIndEndNotActiveStopped = length(veh.timestamp);end

            tempControlCarStatus(tempIndNotActiveStopped:tempIndEndNotActiveStopped) = 0;
        end
    end

    if sum(tempControlCarStatus==-1)>0

        accelerationB = vehSpeed(2:end)-vehSpeed(1:end-1);
        accelerationF = [accelerationB;accelerationB(end)];
        accelerationB = [accelerationB(1); accelerationB];

        indActiveStoppedStart = find(controlActive & vehSpeed==0 & accelerationB < 0 );
        indActiveStoppedEnd = find((controlActive|[controlActive(2:end);controlActive(end)]) & vehSpeed==0 & accelerationF > 0);
        for ii = 1:length(indActiveStoppedStart)
            tempIndStart = indActiveStoppedStart(ii);
            if ii>length(indActiveStoppedEnd)
                tempIndEnd = length(vehSpeed) - 1;
            else
                tempIndEnd = indActiveStoppedEnd(ii);
            end
            if controlActive(tempIndEnd+1)
                tempControlCarStatus(tempIndStart:tempIndEnd) = ones(tempIndEnd-tempIndStart+1,1);
            else
                tempControlCarStatus(tempIndStart:tempIndEnd) = zeros(tempIndEnd-tempIndStart+1,1);
            end
        end
    end
    tempControlCarStatus(tempControlCarStatus==-1) = 0 ; 

    tempLast30 = [];
    for k =1:length(veh.timestamp)
        tempLast30(k,1) = logical(sum(tempControlCarStatus(max(1,k-300):k))>0);
    end
    veh1.control_car = tempControlCarStatus;
    veh1.control_last30 = tempLast30;
end

