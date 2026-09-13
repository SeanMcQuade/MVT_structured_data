%========================================================================
% Initialize
%========================================================================
clear

processingDay = 16; % TODO: change to an input of the script
[parentDirectory, ~, ~] = fileparts(pwd);
[dataRootDirectory, ~, ~] = fileparts(parentDirectory);

dataFolderPath = fullfile(dataRootDirectory, 'data', 'i24motion', ...
    ['2022-11-', num2str(processingDay)]);

outputPath = fullfile(dataRootDirectory, 'results', 'slim', ...
    ['2022-11-', num2str(processingDay)]);

dataFiles = dir(fullfile(outputPath ,'*.json'));
is_dotfile = startsWith({dataFiles.name},'.');
dataFiles = dataFiles(~is_dotfile);

if length(dataFiles) < 24
    error('I24 base files for the day: %d, Nov. 2022 are missing or incomplete.', processingDay)
end

% Use cell arrays to store results from each file
all_start_cell = cell(24, 1);
all_end_cell = cell(24, 1);

for fileNr = 1:24 
    filenameLoad = fullfile(outputPath, dataFiles(fileNr).name);
    fprintf('Loading and decoding MOTION data file, %d/24 ... ', fileNr); tic
    dataTemp = jsondecode(fileread(filenameLoad));
    fprintf('Done (%0.0fsec).\n',toc)

    lanes_od = load([filenameLoad(1:end-5)  '_orig_dist_lane.mat']);
    lanes_od = lanes_od.dataTemp_lane_orig_dist;

    for i = 1:length(dataTemp)
        dataTemp(i).origin_lane = lanes_od(i).origin_lane;
        dataTemp(i).destination_lane = lanes_od(i).destination_lane;
    end
    
    % Filter for any trajectories where ANY lane change occurred (start or end)
    lane_changes_data = dataTemp([dataTemp.origin_lane] ~= [dataTemp.lane_number] | ...
        [dataTemp.destination_lane] ~= [dataTemp.lane_number]);
    
    num_changes = length(lane_changes_data);
    % PREALLOCATE separate struct arrays for START and END events
    empty_start = struct('lane_number', [], 'lane_change_at_start', [], ...
                          'dist_to_eng_av_at_start', [], 'dist_to_av_at_start', [], ...
                          'x_position_at_start', [], 't_at_start', [], ...
                          'av_at_start',[],'eng_av_at_start',[]);
                          
    empty_end = struct('lane_number', [], 'lane_change_at_end', [], ...
                          'dist_to_eng_av_at_end', [], 'dist_to_av_at_end', [], ...
                          'x_position_at_end', [], 't_at_end', [], ...
                          'av_at_end',[],'eng_av_at_end',[]);
    start_capt = repmat(empty_start, num_changes * 2, 1);
    end_capt = repmat(empty_end, num_changes * 2, 1);
    
    tt_start = 1;
    tt_end = 1;
    
    for trajInd = 1:num_changes
        traj = lane_changes_data(trajInd);
        
        % ==========================================
        % 1. MERGE IN (START OF TRAJECTORY)
        % ==========================================
        if traj.origin_lane ~= traj.lane_number
            
            % Downstream AV Check (Negative Distance)
            if ~isempty(traj.downstream_av_id) 
                start_capt(tt_start).lane_number = traj.lane_number;
                start_capt(tt_start).lane_change_at_start = traj.lane_number - traj.origin_lane;
                start_capt(tt_start).dist_to_av_at_start = traj.distance_to_downstream_av_meters(1);
                if ~isempty(traj.downstream_engaged_av_id) 
                    start_capt(tt_start).dist_to_eng_av_at_start = traj.distance_to_downstream_engaged_av_meters(1);
                    start_capt(tt_start).eng_av_at_start = traj.downstream_engaged_av_id(1);
                end
                start_capt(tt_start).av_at_start = traj.downstream_av_id(1);
                start_capt(tt_start).x_position_at_start = traj.x_position_meters(1);
                start_capt(tt_start).t_at_start = traj.timestamp(1);
                tt_start = tt_start + 1;
            end
            
            % Upstream AV Check (Positive Distance)
            if ~isempty(traj.upstream_av_id) 
                start_capt(tt_start).lane_number = traj.lane_number;
                start_capt(tt_start).lane_change_at_start = traj.lane_number - traj.origin_lane;
                start_capt(tt_start).dist_to_av_at_start = traj.distance_to_upstream_av_meters(1);
                start_capt(tt_start).av_at_start = traj.upstream_av_id(1);
                if ~isempty(traj.upstream_engaged_av_id) 
                    start_capt(tt_start).dist_to_eng_av_at_start = traj.distance_to_upstream_engaged_av_meters(1);
                    start_capt(tt_start).eng_av_at_start = traj.upstream_engaged_av_id(1);
                end
                start_capt(tt_start).x_position_at_start = traj.x_position_meters(1);
                start_capt(tt_start).t_at_start = traj.timestamp(1);
                tt_start = tt_start + 1;
            end
        end
        
        % ==========================================
        % 2. MERGE OUT (END OF TRAJECTORY)
        % ==========================================
        if traj.destination_lane ~= traj.lane_number
            
            % Downstream AV Check (Negative Distance)
            if ~isempty(traj.downstream_av_id) 
                end_capt(tt_end).lane_number = traj.lane_number;
                end_capt(tt_end).lane_change_at_end = traj.destination_lane - traj.lane_number;
                end_capt(tt_end).dist_to_av_at_end = traj.distance_to_downstream_av_meters(end);
                if ~isempty(traj.downstream_engaged_av_id) 
                    end_capt(tt_end).dist_to_eng_av_at_end = traj.distance_to_downstream_engaged_av_meters(end);
                    end_capt(tt_end).eng_av_at_end = traj.downstream_engaged_av_id(end);
                end
                end_capt(tt_end).av_at_end = traj.downstream_av_id(end);
                end_capt(tt_end).x_position_at_end = traj.x_position_meters(end);
                end_capt(tt_end).t_at_end = traj.timestamp(end);
                tt_end = tt_end + 1;
            end
            
            % Upstream AV Check (Positive Distance)
            if ~isempty(traj.upstream_av_id) 
                end_capt(tt_end).lane_number = traj.lane_number;
                end_capt(tt_end).lane_change_at_end = traj.destination_lane - traj.lane_number;
                if ~isempty(traj.upstream_engaged_av_id) 
                    end_capt(tt_end).dist_to_eng_av_at_end = traj.distance_to_upstream_engaged_av_meters(end);
                    end_capt(tt_end).eng_av_at_end = traj.upstream_engaged_av_id(end);
                end
                end_capt(tt_end).av_at_end = traj.upstream_av_id(end);
                end_capt(tt_end).dist_to_av_at_end = traj.distance_to_upstream_av_meters(end);
                end_capt(tt_end).x_position_at_end = traj.x_position_meters(end);
                end_capt(tt_end).t_at_end = traj.timestamp(end);
                tt_end = tt_end + 1;
            end
        end
    end
    
    % Trim the unused preallocated rows for both structs
    start_capt(tt_start:end) = [];
    end_capt(tt_end:end) = [];
    
    if exist('start_capt', 'var') && ~isempty(start_capt)
        all_start_cell{fileNr} = start_capt;
    end
    
    if exist('end_capt', 'var') && ~isempty(end_capt)
        all_end_cell{fileNr} = end_capt;
    end
end

% Combine the cell arrays into two distinct struct arrays
all_lane_changes_start = vertcat(all_start_cell{:});
all_lane_changes_end = vertcat(all_end_cell{:});
save(fullfile(outputPath, ['LC_data_' char(num2str(processingDay)) '.mat'])...
,'all_lane_changes_end','all_lane_changes_start')