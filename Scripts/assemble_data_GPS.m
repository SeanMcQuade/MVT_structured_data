function [] = assemble_data_GPS(processingDay)
% This function processes AV GPS data from the MVT and produces JSON file for a
% single MVT test day.
% (C) 2025 by Sulaiman Almatrudi
if nargin < 1
error(['Specify the day of Nov. 2022 MVT to generate AVs '...
        'GPS data files (from 16 to 18) to run assemble_data_GPS ']);
end
fprintf('Starting the assembly of AV GPS data for %dth Nov. 2022\n', processingDay);
%========================================================================
% Parameters
%========================================================================
%%%% Threshold choices for matching AVs GPS to I24 MOTION trajectories
maxMatchDist = 6; %[m] maximum distance between AVs GPS trajectories and potential matching I24 trajectories
maxMatchSpdDiff = 2; % [m/s] maximum speed difference
maxMatchLaneDiff = 0.5; %[lane] maximum estimated lane difference
minMatchTime = 3; % [s] minimum matching time
originXPosition = 309804.0625; % [ft] location of the origin for the x-coordinates 
ft2meterFactor = 0.3048; % [m/ft] conversion factor from feet to meter
%========================================================================
% Load, parse, and preprocess AVs GPS data
%========================================================================
% Get file path of base GPS data
[parentDirectory, ~, ~] = fileparts(pwd);
gpsFolderPath = fullfile(parentDirectory,'Data','Data_GPS') ;
if ~isfolder(gpsFolderPath)
    error('Folder %s does not exist.\n',gpsFolderPath)
end
% Parse AV files into trajectories of individual runs
fprintf('Parsing AV data... ') ; tic
dataGPS0 = parse_gps_data(gpsFolderPath,processingDay);
nAVruns = length(dataGPS0);
fprintf('Done (%0.0fsec).\n',toc)
% Pre-process av files
fprintf('Pre-processing AV data... ');tic
dataGPS0 = preproc_gps(dataGPS0,gpsFolderPath);
fprintf('Done (%0.0fsec).\n',toc)
% Find AV activity limits to load MOTION data
% Get epoch time of 6am and 10am for the day.
dataTLimits = datetime(2022,11,processingDay,[6 10],0,0,'TimeZone','America/Chicago');
dataTLimits = double(convertTo(dataTLimits,'epochtime'));
% Find start of first AV run during the test
minAVStart = min([dataGPS0([dataGPS0.ending_time]>dataTLimits(1)).starting_time]);
minFileNr = max(1,floor((minAVStart-dataTLimits(1))/60/10));
% Find end of last AV run during the test
maxAVStart = max([dataGPS0([dataGPS0.starting_time]<dataTLimits(2)).ending_time]);
maxFileNr = min(24,floor((maxAVStart-dataTLimits(1))/60/10)+1);
% Find I24 MOTION data files in the folder
dataFolderPath = fullfile(parentDirectory,'Data',['Data_2022-11-' num2str(processingDay) ...
    '__I24_Base']) ;
if ~isfolder(dataFolderPath)
    error('Folder %s does not exist.\n',dataFolderPath)
end
% Use day abbreviation to identify base data files
dayAbbrvs = ["wed","thu","fri"];
dayAbbrv = dayAbbrvs(processingDay-15);
dataFiles = dir(fullfile(dataFolderPath ,['*_' char(dayAbbrv) '_0_*.json']));
if length(dataFiles) < 24
    error('I24 base files for the day: %d, Nov. 2022 are missing or incomplete.',processingDay)
end
%========================================================================
% Find I24 trajectories matching AV GPS trajectories
%========================================================================
% Initialize 
matchedSegments(10*length(dataGPS0),1).av_trj_i = [];
matchedSegments(10*length(dataGPS0),1).x_diff =   [];
matchedSegments(10*length(dataGPS0),1).timestamp = [];
segmentsCounter = 1;
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
        fprintf('Possible AV trajectories not found in file \n') ; tic
        continue;
    else
        fprintf('Processing file for possible AV matches... ') ; tic
    end
    % Delete unnecessary fields in MOTION data for matching 
    dataTemp = rmfield(dataTemp,{'x_id','flags','length','width','height',...
        'merged_ids','road_segment_ids','coarse_vehicle_class','fine_vehicle_class'...
        ,'fragment_ids','compute_node_id','local_fragment_id','configuration_id',...
        'x_score','y_score'});
    % Identify and assign lanes to MOTION data
    dataLanes = assign_lanes(dataTemp);
    % Loop over AVs runs concurrent to the MOTION data in file and find 
    % matched segments each AV run
    for indAVsInFile =1:length(dataAVsInFile)
        avVeh =   dataAVsInFile(indAVsInFile);
        avX = avVeh.x_position;
        avLane = avVeh.assigned_lane;
        avT = avVeh.timestamp;
        % Find MOTION trajectories that could possibly match AV
        possTrajsMatch = find([dataTemp.first_timestamp] <= avT(end) ...
            & [dataTemp.last_timestamp] >= avT(1) & [dataTemp.direction] == avVeh.direction);
        % Loop over MOTION trajectories possibly matching the AV and
        % evaluate the match
        for indPossMatch = 1:length(possTrajsMatch)
            traj = dataTemp(possTrajsMatch(indPossMatch));
            trajLane =  dataLanes(possTrajsMatch(indPossMatch)).lane;
            trajX = (traj.x_position-originXPosition);
            trajT = traj.timestamp;
            trajV = diff(trajX)./diff(trajT);
            trajV = [trajV(1);trajV]*ft2meterFactor;
            % Calculate distance to AV in meters
            distToAV = (interp1(avT,avX,traj.timestamp) - trajX) * ft2meterFactor;
            avgDistToAVX = distToAV;
            avgDistToAVX = avgDistToAVX(~isnan(avgDistToAVX));
            avgDistToAVX = mean(abs(avgDistToAVX));
            % Calculate difference in lane to AV
            diffToAVLane = avLane-trajLane;
            avgdiffToAVLane = diffToAVLane;
            avgdiffToAVLane = avgdiffToAVLane(~isnan(avgdiffToAVLane));
            avgdiffToAVLane = mean(abs(avgdiffToAVLane));
            % Ignore trajectory if it has an average absolute x position difference
            % of > 200 m or average abs. lane difference of >= 2 lanes.
            % assume too far to be matching
            if avgDistToAVX <= 200 && avgdiffToAVLane < 2
                % calculate and smooth acceleration 
                vDiff = diff(distToAV)./diff(trajT);
                vDiff = [vDiff(1);vDiff];
                vDiffS = abs(smoothdata(vDiff, 'gaussian', 3, 'SamplePoints', trajT));
                timeMatching = 0;
                % Loop over trajectory and track time of continuous matching
                % to AV. If a segment is matching for over minMatchTime,
                % append it to matched segments.
                for stepInTraj = 2:length(trajT)
                    if abs(distToAV(stepInTraj)) <= maxMatchDist ...
                            && abs(diffToAVLane(stepInTraj)) <= maxMatchLaneDiff ...
                            && vDiffS(stepInTraj) <= maxMatchSpdDiff
                        if timeMatching == 0
                            matchStartInd = stepInTraj;
                        end
                        timeMatching = timeMatching + trajT(stepInTraj) - trajT(stepInTraj-1);
                    else
                        if timeMatching >= minMatchTime
                            matchedSegments(segmentsCounter).av_trj_i = avVeh.index;
                            matchedSegments(segmentsCounter).x_diff =   distToAV(matchStartInd:stepInTraj);
                            matchedSegments(segmentsCounter).timestamp = trajT(matchStartInd:stepInTraj);
                            segmentsCounter = segmentsCounter+1;
                        end
                        matchStartInd = [];
                        timeMatching = 0;
                    end
                end
                if timeMatching >= minMatchTime
                    matchedSegments(segmentsCounter).av_trj_i = avVeh.index;
                    matchedSegments(segmentsCounter).x_diff =   distToAV(matchStartInd:stepInTraj);
                    matchedSegments(segmentsCounter).timestamp = trajT(matchStartInd:stepInTraj);
                    segmentsCounter = segmentsCounter+1;
                end
            end

        end
    end
    fprintf('Done (%0.0fsec).\n',toc)
end
matchedSegments(segmentsCounter:end) = [];
% Calculate median offset between GPS runs and I24 matched signals
xd_match = [];
xd_match(length(dataGPS0),1).xd = [];
xd_match(length(dataGPS0),1).median_xd = [];
for av_ri =1 :length(dataGPS0)
    runMatchedSegments = [matchedSegments.av_trj_i]==av_ri;
    xd_match(av_ri).xd = vertcat( matchedSegments(runMatchedSegments).x_diff);
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
avPingsData = readtable(fullfile(gpsFolderPath, ['veh_ping_202211' ...
    num2str(processingDay) '.csv']));
avVINs = readtable(fullfile(gpsFolderPath ,'veh_vins.csv'));
avConnectionStatus  = get_connection_status(avPingsData,avVINs);
nAVs = length(dataGPS0);
clear dataGPS
% Initialize output data structure
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
% Loop over processed AV runs and append to output data structure 
for avInd =1:nAVs
    avVeh = dataGPS0(avInd);
    % Locate start of the run in the server data
    runStartInServer = find(avConnectionStatus(avVeh.vin).timestamp < avVeh.timestamp(1),1,'last');
    isConnected = (0*avVeh.timestamp);
    if ~isempty(runStartInServer)
        % Get server ping times during the run
        runTimeConnected = avConnectionStatus(avVeh.vin).timestamp(runStartInServer:end);
        % Loop over every timestep in the run and check if there has been a
        % ping to the server in the last (maxPingInterval) seconds
        if runTimeConnected(end) > avVeh.timestamp(1)
            for timeStp = 1:length(avVeh.timestamp)
                tempInd = find(runTimeConnected<=avVeh.timestamp(timeStp),1,'last');
                if isempty(tempInd)
                    tempInd = 1;
                end
                isConnected(timeStp) = single((runTimeConnected(tempInd) ...
                    - avVeh.timestamp(timeStp)) >= -2); 
            end
        end

    end
    % Append parsed data to output structure to write to JSON file
    dataGPS(avInd).av_id = avVeh.vin;
    dataGPS(avInd).assigned_lane = avVeh.assigned_lane;
    dataGPS(avInd).direction = avVeh.direction;
    dataGPS(avInd).timestamp = round(avVeh.timestamp,6,'decimals');
    dataGPS(avInd).latitude = round(avVeh.Lat,6,'decimals');
    dataGPS(avInd).longitude = round(avVeh.Long,6,'decimals');
    %remove bias detected by  matching to I24MOTION data
    dataGPS(avInd).x_position = round(ft2meterFactor * avVeh.x_position ...
        - xd_match(avInd).median_xd,6,'decimals'); 
    dataGPS(avInd).y_position = round(ft2meterFactor * avVeh.y_position,6,'decimals');
    dataGPS(avInd).controller_engaged = avVeh.control_active;
    dataGPS(avInd).speed = avVeh.can_speed;
    dataGPS(avInd).is_server_connected = isConnected;
    dataGPS(avInd).first_timestamp = avVeh.timestamp(1);
    dataGPS(avInd).last_timestamp = avVeh.timestamp(end);
    % Append modified control status signal to correct for controller
    % showing as inactive when the vehicle stops
    dataGPS(avInd) = get_control_car_status(dataGPS(avInd));
end
% Clip trajectories at testbed limits and prepare data for writing
XLIMS = [-400 2.54e4] *ft2meterFactor; % testbed limits
gpsDataFields = string(fieldnames(dataGPS));
% Loop over all runs and clip at testbed limits
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
filenameWrite = ['CIRCLES_GPS_10Hz_2022-11-' num2str(processingDay) '.json'];
fprintf('Encoding and writing file %s to file... ',filenameWrite);tic
jsonStr = jsonencode(dataGPS);
fid = fopen(fullfile(gpsFolderPath, filenameWrite), 'w');
fwrite(fid, jsonStr, 'char');
fclose(fid);
clear jsonStr
fprintf('Done (%0.0fsec).\n',toc)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Local Functions Definitions %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dataGPS0 = parse_gps_data(gpsDataFolder,processingDay)
% Function to parse GPS data from base csv data files
%========================================================================
% Parameters for parsing GPS runs
%========================================================================
maxRcsx = 2.91e4; % [ft] maximum x position for testbed
minRcsx = -3000; % [ft] minimum x position for testbed
maxRcsy = 150; % [ft] maximum y position to clip GPS trajectories
maxGpsInactive = 0.01; % maximum percentage of run length with inactive GPS signal
minRunLength = 1.2; % [km] minimum run length along x axis
minRunTime = 60; % [s] minimum run time
km2ftFactor = 3281; % conversion factor from km to ft
%========================================================================
% Initialize
%========================================================================
% Define test day time limits to extract day data
dayLimits = datetime(2022,11,processingDay,[3 18],0,0,'TimeZone','America/Chicago');
dayLimits = convertTo(dayLimits,'epochtime');
vehsFiles = dir(fullfile(gpsDataFolder,'circles_v2_1_car*.csv'));
% Detect and ensure  correct import options for base AV csv files
opts1 = detectImportOptions(fullfile(gpsDataFolder , vehsFiles(1).name));
opts1.VariableTypes(13) = {'char'};
% Identify AVs with available base data files
nrFiles = length(vehsFiles);
IDsWithFiles = zeros(1,nrFiles);
for iFile = 1:nrFiles
    tempVin = vehsFiles(iFile).name(end-6:end-4); % look at end of .csv file name
    tempVin =  tempVin(tempVin>='0' & tempVin<='9'); % consider only integer char
    IDsWithFiles(iFile) = str2double(tempVin);
end
fprintf(['\n' num2str(nrFiles) ' vehicles GPS files found...\n'])
%========================================================================
% Load and parse GPS data for each AV for the specified day
%========================================================================
% Track the number of extracted runs
runCounter = 1;
dataGPS0 = {};
% Loop over AVs
for avID = 1:103
    runNum = 0;
    % Load data file
    if ismember(avID,IDsWithFiles)
        vehTable = readtable(fullfile(gpsDataFolder,['circles_v2_1_car' num2str(avID) '.csv']),opts1);
    else
        continue
    end
    % Filter data for the day by day time limits
    vehTable = vehTable(vehTable.Systime > dayLimits(1) & vehTable.Systime < dayLimits(2),:);
    if isempty(vehTable)
        fprintf(['no data found for car ' num2str(avID) ' on day ' num2str(processingDay) ...
            ', car assumed inactive for the day \n'])
        continue
    end
    % Parse runs of the AV in east and west bound traffic
    diffY = diff(vehTable.rcs_y);
    diffY = [diffY(1);diffY];
    % Locate beginning  of the run
    runStart = find(sign(diffY) .* sign( vehTable.rcs_y) < 0 ... % moving towards the center of the highway (entering)
        & abs(vehTable.rcs_y) < maxRcsy & vehTable.rcs_x > minRcsx ... % within testbed limits
        & vehTable.rcs_x < maxRcsx ,1);
    vehTable = vehTable(runStart:end,:); % discard data prior to the start of the run
    % Locate end of individual run
    runEnd = find(abs(vehTable.rcs_y) > maxRcsy | vehTable.rcs_x < minRcsx ...
        | vehTable.rcs_x > maxRcsx,1); % vehicle exceeded testbed limits
    runTime =  vehTable(1:runEnd,:).Systime; % timestamps of the run
    runX =  vehTable(1:runEnd,:).rcs_x; % [ft] x position of the run
    status = cell2mat(vehTable(1:runEnd,:).Status); % GPS unit status
    % Loop over the AV data to extract all runs
    while size(vehTable,1)>1
        if runTime(end)-runTime(1) > minRunTime && ...
                abs(runX(end)-runX(1)) > minRunLength * km2ftFactor && ...
                sum(double(status)>80) < maxGpsInactive*length(runX)
            % Append feasible runs
            runNum = runNum+1; % track number of runs per AV
            runDirection = sign(runX(end) - runX(1));  % east(+)/west(-)
            runY = - vehTable(1:runEnd,:).rcs_y; % extract y position and correct sign
            dataGPS0(runCounter).vin = avID;
            dataGPS0(runCounter).run_num = runNum;
            dataGPS0(runCounter).timestamp = runTime;
            dataGPS0(runCounter).direction = runDirection;
            dataGPS0(runCounter).x_position = runX;
            dataGPS0(runCounter).y_position = runY;
            dataGPS0(runCounter).starting_time = runTime(1);
            dataGPS0(runCounter).ending_time = runTime(end);
            dataGPS0(runCounter).Long = vehTable(1:runEnd,:).Long;
            dataGPS0(runCounter).Lat = vehTable(1:runEnd,:).Lat;
            dataGPS0(runCounter).state_x = vehTable(1:runEnd,:).state_x;
            dataGPS0(runCounter).state_y = vehTable(1:runEnd,:).state_y;
            dataGPS0(runCounter).can_speed = vehTable(1:runEnd,:).can_speed;
            % Translate control status from string to boolean
            cont_st = string(vehTable(1:runEnd,:).control_active);
            dataGPS0(runCounter).control_active = (cont_st == "True");
            runCounter = runCounter+1;
        end
        if runEnd < size(vehTable,1)
            % Discard run base data after parsing while there's remaining data
            vehTable = vehTable(runEnd+1:end,:);
            % Identify run start and end (similar to previously)
            diffY = diff(vehTable.rcs_y);
            diffY = [diffY(1);diffY];
            runStart = find(sign(diffY) .* sign( vehTable.rcs_y) < 0 ...
                & abs(vehTable.rcs_y) < maxRcsy & vehTable.rcs_x > minRcsx & vehTable.rcs_x < maxRcsx ,1);
            vehTable = vehTable(runStart:end,:);
            runEnd = find(abs(vehTable.rcs_y) > maxRcsy | vehTable.rcs_x < minRcsx ...
                | vehTable.rcs_x > maxRcsx,1);
            if isempty(runEnd); runEnd =size(vehTable,1);end
            % Get new run timestamps, x position and GPS status
            runTime =  vehTable(1:runEnd,:).Systime;
            runX = vehTable(1:runEnd,:).rcs_x;
            status = cell2mat(vehTable(1:runEnd,:).Status);
        else
            % Discard base data after last run
            vehTable = [];
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dataGPS0 = preproc_gps(dataGPS0,gpsDataFolder)
% Function to preprocess GPS data and sample the base data at 10
% samples/sec
inVehShift = 7; % [ft] x position shift from GPS unit position to AV rear bumper
nAVrunss = length(dataGPS0);
avVINData = readtable(fullfile(gpsDataFolder ,'veh_vins.csv'));
for runIdx = 1:nAVrunss % loop over AV runs and sample data at 10 hz
    % Initialize time grid at 10hz
    t10Hz = floor(dataGPS0(runIdx).timestamp(1)*10)/10:.1:ceil(dataGPS0(runIdx).timestamp(end)*10)/10;
    dataGPS0(runIdx).index = runIdx;
    % Shift x position data from average GPS position to rear bumper of AV
    dataGPS0(runIdx).x_position = dataGPS0(runIdx).x_position - dataGPS0(runIdx).direction * inVehShift;
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
    dataGPS0(runIdx).assigned_lane = avVINData(dataGPS0(runIdx).vin,:).lane_num;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newValue = sample_10hz(orgT, orgValue, newT)
% Function to sample data at 10 hz, Function also identifies and removes
% erroneous data records.
% Input:    orgT: original data timestamps [double array]
%           orgValue: original data values [double array]
%           newT: sampling time grid [double array]
%
% Output:   newValue: data sampled on the new grid [double array]
%
B = orgT;
% Only consider data recorded at < 0.01s from expected timestamps at 10hz 
idx = [true;diff(B)>0.09&diff(B)<0.11];
B = B(idx);
A = orgValue(idx);
newValue = interp1(B,A,newT,'linear','extrap'); % interploate data on 10hz grid
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dataLanes] = assign_lanes(baseFileName)

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
%%%% Processing options for lane identification in the v2 data %%%%
laneIdentificationOpts.LANEWIDTH = 12;       % Lane width in feet
laneIdentificationOpts.Nr_XCELLS = 200;      % Number of cells in the x direction used to identify driving line
laneIdentificationOpts.MINCELLSAMPLES = 20;  % Minimum number of samples per cell used to identify driving line
laneIdentificationOpts.Y_UP_LIM = 5;         % samples upper bound as a multiplier of lane width to remove outliers
laneIdentificationOpts.Y_LOW_LIM = 0.5;      % samples lower bound as a multiplier of lane width to remove outliers
% Shift parameters for y position correction (westbound) s.t. y_corrected =  Sw * (y  - drivingLineShift_w) + Cw;
laneIdentificationOpts.Se = 0.97;           
laneIdentificationOpts.Ce = 1;
% Shift parameters for y position correction (eastbound) s.t. y_corrected =  Se * (y  - drivingLineShift_e) + Ce;
laneIdentificationOpts.Sw = 0.98;            
laneIdentificationOpts.Cw = 1;
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
xEast = []; yEast = [];
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
        else                %east bound
            xEast = [xEast, veh.x_position(1 :ss: end)'];
            yEast = [yEast, veh.y_position(1 :ss: end)'];
        end

    end
end
% Copy defined processing parameters 
laneWidth = laneIdentificationOpts.LANEWIDTH;   
nrXCells = laneIdentificationOpts.Nr_XCELLS;   
minCellsSamples = laneIdentificationOpts.MINCELLSAMPLES;
yUpLim = laneIdentificationOpts.Y_UP_LIM;
yLowLim = laneIdentificationOpts.Y_LOW_LIM;
upBoundY = laneWidth*yUpLim;
lowBoundY = laneWidth*yLowLim;
% Remove outliers from the sampled data based on bounds
yEast = -yEast;  % flip eastbound y position to be positive in sign
xEast = xEast(yEast<=upBoundY & yEast>=lowBoundY);
yEast = yEast(yEast<=upBoundY & yEast>=lowBoundY);
xWest = xWest(yWest<=upBoundY & yWest>=lowBoundY);
yWest = yWest(yWest<=upBoundY & yWest>=lowBoundY);
% Initialize x position grid in both directions
xCellsEast = linspace(min(xEast),max(xEast),nrXCells+1);
xCellsWest = linspace(min(xWest),max(xWest),nrXCells+1);
% Initialize driving line in both directions
drivingLineWest = zeros(1,nrXCells);
drivingLineEast = zeros(1,nrXCells);
% Loop over x cells to calculate the shift in driving line
for cellInd = 1:nrXCells
    % Identify points in the cells
    pointsInCellWest = xCellsWest(cellInd)<xWest&xWest<=xCellsWest(cellInd+1);
    pointsInCellEast = xCellsEast(cellInd)<xEast&xEast<=xCellsEast(cellInd+1);
    % Aggregate samples in the cell. Shift down so leftmost lane edge = 0 
    yCellWest = yWest(pointsInCellWest) - laneWidth/2; 
    yCellEast = yEast(pointsInCellEast) - laneWidth/2;
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
    % Similar process for Eastbound data
    if length(yCellEast) > minCellsSamples
        p_j_e = cos(2*pi*yCellEast/laneWidth); 
        q_j_e = sin(2*pi*yCellEast/laneWidth); 
        dle_t = atan2(mean(q_j_e), mean(p_j_e))*laneWidth/2/pi; 
        drivingLineEast(cellInd) = dle_t;
    else
        drivingLineEast(cellInd) = 0;
    end
end
% Remove maximum edge in the grid
xCellsWest = xCellsWest(1:end-1);
xCellsEast = xCellsEast(1:end-1);
%%% Start Lane assignment
% Set shifting factors used for y-position correction
Se = laneIdentificationOpts.Se; Ce = laneIdentificationOpts.Ce; 
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
    if veh.direction >0        % take positive y value for eastbound
        y = -y;
    end
    % correct for wiggle in y
    if veh.direction >0 % eastbound
        dataLanes(trajInd).y_corr = - ...
            (Se * (y - interp1(xCellsEast,drivingLineEast,x,'linear','extrap'))+Ce);
    else
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
function avConnectionStatus  = get_connection_status(avPingsData,avVINs)
% Function that takes in base vehicle server data and outputs vehicle connection
% status to the server.
%Inputs:    avPingsData: Server ping data [Table]
%           avVINs: AVs VIN data [Table]
%
%Outputs:   avConnectionStatus: Connection status to the server for each Av
%                               [Boolean array]
%
avPingsData.av_id = 0*avPingsData.gpstime;
% Loop over AVs and append av ID to each AV ping data
for i = 1:size(avVINs,1)
    VIN = avVINs(i,:).vin;
    % Locate ping data for the AV
    vinDataPoints = strcmp(string(avPingsData.vin),string(VIN));
    % Append AV id to the ping data logs
    avPingsData(vinDataPoints,:).av_id = avVINs(i,:).veh_id*ones(sum(vinDataPoints),1);
end
% Initialize output table
avConnectionStatus = {};
avConnectionStatus.vin = [];
avConnectionStatus.timestamp = [];
avConnectionStatus.status = [];
for avID =1:103
    avConnectionStatus(avID).vin = avID;
    % Identify where logged pings are not empty
    rawStatus = avPingsData(avPingsData.av_id==avID,:).acc_status;
    avConnectionStatus(avID).status = ~isnan(rawStatus);
    % append time stamp from gpstime given in ms
    avConnectionStatus(avID).timestamp = avPingsData(avPingsData.av_id==avID,:).gpstime/1e3;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function veh1  = get_control_car_status(veh)
% Function that appends updated control status to account for when
% controller is assumed inactive whenever AVs are stopped. 
veh1 = veh;
% calculate speed from position if CAN speed is unavailable
if sum(veh.speed) == 0 
    x = veh.x_position(:);
    t = veh.timestamp(:);
    vehSpeed0 = (x(2:end) - x(1:end-1))./(t(2:end) - t(1:end-1));
    % Round to nearest 0.5 
    vehSpeed0 = round(abs([vehSpeed0;vehSpeed0(end)])*2)/2;
    % Pad data for filtering 
    windowSize = 10; % Window size for filtering
    paddedSpeed = [zeros(windowSize/2,1);vehSpeed0;zeros(windowSize/2,1)];
    % Apply moving median filter 
    vehSpeed = movmedian(paddedSpeed,windowSize);
    % Remove padding
    vehSpeed = vehSpeed(windowSize/2+1:end-windowSize/2);
    veh.speed = vehSpeed;
end
vehSpeed = veh.speed(:);
tempControlCarStatus = -1*ones(size(veh.controller_engaged));
controlActive = veh.controller_engaged(:); % original control status signal
indActive = find(controlActive);
% Set status to active when vehicle is moving and orig. signal is active
tempControlCarStatus(indActive) = 1; 
% Set status to inactive when vehicle is moving and orig. signal is inactive
indNotActiveMoving = ~controlActive & vehSpeed~=0;
tempControlCarStatus(indNotActiveMoving) = 0;
% Process cases where status is turned inactive while vehicle is stopped
indActiveStopped = controlActive & vehSpeed==0;
indActiveStopped = find(indActiveStopped);
% Identify where original control status signal is inactive due to the vehicle 
% being stopped but the controller is turned on 
if ~isempty(indActiveStopped)
    % Find where the original control status signal is turned off after
    % stopping and loop over each instance
    indActiveStoppedEdge = indActiveStopped(diff([indActiveStopped;length(controlActive)])>1);
    for indTemp = indActiveStoppedEdge(:).'
        % Find nearest point where original signal turns active 
        tempIndEndActiveStopped = indActive(find(indActive > indTemp,1))-1;
        if isempty(tempIndEndActiveStopped) 
            tempIndEndActiveStopped = length(veh.timestamp);
        end        
        if tempIndEndActiveStopped == indTemp
            continue
        end
        % If original control status signal is turned active within (maxSteps) 
        % of vehicle starting to move, assume controller was on during stoppage.
        maxSteps = 15;
        if sum(vehSpeed(indTemp:tempIndEndActiveStopped)>0) < ...
                min(maxSteps,tempIndEndActiveStopped-indTemp+2)
            tempControlCarStatus(indTemp:tempIndEndActiveStopped) = 1;
        else
            tempControlCarStatus(indTemp:tempIndEndActiveStopped) = 0;
        end
    end
end
% Find instances where original signal is inactive and the vehicle is
% stopping
indNotActiveStopped = find(~controlActive(1:end-1) & vehSpeed(1:end-1)~=0 & vehSpeed(2:end)==0);
if ~isempty(indNotActiveStopped)
    % Loop over every instance of vehicle coming to a stop already inactive, ensure
    % new signal is inactive 
    for tempIndNotActiveStopped = indNotActiveStopped(:).'
        % Find where vehicle starts moving or signal turns active
        tempIndEndNotActiveStopped = find(controlActive | vehSpeed~=0);
        tempIndEndNotActiveStopped = tempIndEndNotActiveStopped(find(tempIndEndNotActiveStopped > tempIndNotActiveStopped,1))-1;
        if isempty(tempIndEndNotActiveStopped) ; tempIndEndNotActiveStopped = length(veh.timestamp);end
        % Ensure new signal is inactive during vehicle manual stoppage
        tempControlCarStatus(tempIndNotActiveStopped:tempIndEndNotActiveStopped) = 0;
    end
end
tempControlCarStatus(tempControlCarStatus==-1) = 0 ; % assign inactive for unresolved states
% Identify if the new signal is active in the last 30 seconds of each timestep
tempLast30 = [];
for timeS =1:length(veh.timestamp)
    tempLast30(timeS,1) = logical(sum(tempControlCarStatus(max(1,timeS-300):timeS))>0);
end
veh1.control_car = tempControlCarStatus;
veh1.control_last30 = tempLast30;
end