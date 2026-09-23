function [] = extract_lane_changes_v_dist_to_av(processingDay, varargin)
% EXTRACT_LANE_CHANGES_V_DIST_TO_AV  Lane-change events and their distance to AVs.
% (C) 2026 CIRCLES Energy team
%
% Purpose
%   Collects every lane change in a day's slim data set, recording where each
%   one happened relative to the nearest (engaged) control vehicle. Merge-in
%   events are taken at the start of a clipped trajectory and merge-out events
%   at its end, which is where the clipping placed the lane change.
%
% Inputs
%   processingDay  16, 17, or 18 (November 2022)
%   varargin       options struct and/or name/value pairs (see mvt.options);
%                  Force, Clean, DryRun, Verbose
%
% Outputs
%   <results>/analysis/2022-11-DD/LC_data_DD.mat, holding
%     all_lane_changes_start  merge-in events (one row per AV the event was
%                             measured against, upstream and downstream)
%     all_lane_changes_end    merge-out events, likewise
%
% Algorithm
%   1. Resolve the day's 24 slim segments and their lane sidecars from the
%      manifest, and rebuild only when the single output is older than any of
%      those 48 inputs or than the code (mvt.isStale).
%   2. Per segment, pair the slim trajectories with the sidecar *by index* and
%      keep the ones whose origin or destination lane differs from the lane
%      they were driven in.
%   3. Record each event's distance to the nearest upstream and downstream AV,
%      and to the nearest engaged ones, at the moment of the change.
%   4. Concatenate the day and write it through a temporary name.
%
% Notes
%   This stage decodes all 24 slim segments (several GB) in one process; it is
%   deliberately not sharded, so the whole day is one unit of work.
%
% Dependencies
%   mvt.options, mvt.dayDir, mvt.expectedOutputs, mvt.isStale,
%   mvt.sources, mvt.laneSidecarName, mvt.atomicSave, mvt.ensureDir,
%   mvt.progress, mvt.log
%
if nargin < 1
    error(['Specify the day of Nov. 2022 MVT to extract lane changes ' ...
        '(from 16 to 18)']);
end
mvt.assertDay(processingDay)
opts = mvt.options(varargin{:});

%========================================================================
% Initialize
%========================================================================
slimPath = mvt.dayDir('slim', processingDay);
outputPath = mvt.dayDir('analysis', processingDay);
mvt.ensureDir(outputPath)

% Both lists come from mvt.expectedOutputs rather than mvt.manifest. This stage
% reads the slim trajectories and their sidecars, never the raw data, so it must
% work for a reader who downloaded those two and has no raw segments to build a
% manifest from.
slimFiles = mvt.expectedOutputs('slim', processingDay, opts);
laneFiles = mvt.expectedOutputs('lanes', processingDay, opts);
if numel(slimFiles) < 24
    error('I24 base files for the day: %d, Nov. 2022 are missing or incomplete.', processingDay)
end
if numel(laneFiles) ~= numel(slimFiles)
    error('mvt:extract_lane_changes:sidecarCount', ...
        ['2022-11-%d has %d slim segments but %d lane sidecars. Run ', ...
        '`make lanes-%d` to rebuild them.'], processingDay, numel(slimFiles), ...
        numel(laneFiles), processingDay);
end
segmentNames = cell(1, numel(slimFiles));
for iSeg = 1:numel(slimFiles)
    [~, segName, segExt] = fileparts(slimFiles{iSeg});
    segmentNames{iSeg} = [segName segExt];
end

outputFile = mvt.expectedOutputs('lc', processingDay, opts);
outputFile = outputFile{1};
[stale, staleReason] = mvt.isStale(outputFile, [slimFiles, laneFiles], ...
    mvt.sources('extract_lane_changes_v_dist_to_av', opts), opts);
if ~stale
    mvt.log(opts, 'skip LC_data_%d.mat: %s', processingDay, staleReason);
    return
end
mvt.log(opts, 'build LC_data_%d.mat: %s', processingDay, staleReason);
if opts.Clean && isfile(outputFile) && ~opts.DryRun
    delete(outputFile)
end
if opts.DryRun
    return
end

% Use cell arrays to store results from each file
all_start_cell = cell(numel(slimFiles), 1);
all_end_cell = cell(numel(slimFiles), 1);

reportProgress = mvt.progress(numel(slimFiles), ...
    sprintf('lc 2022-11-%d', processingDay), 'Opts', opts);
for fileNr = 1:numel(slimFiles)
    filenameLoad = slimFiles{fileNr};
    if ~isfile(filenameLoad)
        error('mvt:extract_lane_changes:missingSlim', ...
            'Missing slim segment %s. Run `make slim-%d` first.', ...
            filenameLoad, processingDay);
    end
    if ~isfile(laneFiles{fileNr})
        error('mvt:extract_lane_changes:missingSidecar', ...
            'Missing lane sidecar %s. Run `make lanes-%d` first.', ...
            laneFiles{fileNr}, processingDay);
    end
    fprintf('Loading and decoding MOTION data file, %d/%d ... ', ...
        fileNr, numel(slimFiles)); tic
    dataTemp = jsondecode(fileread(filenameLoad));
    fprintf('Done (%0.0fsec).\n',toc)

    lanes_od = load(laneFiles{fileNr});
    lanes_od = lanes_od.dataTemp_lane_orig_dist;

    % The sidecar is paired with the slim JSON by index, which only holds while
    % both were built from the same raw segment with the same clipping. A
    % length mismatch is the one symptom that surfaces cheaply, and silently
    % mispaired lanes are worse than a stopped build.
    if numel(lanes_od) ~= numel(dataTemp)
        error('mvt:extract_lane_changes:sidecarMismatch', ...
            ['%s has %d trajectories but %s has %d. They must come from the ', ...
            'same build; re-run `make lanes-%d` (or FORCE=1) and try again.'], ...
            segmentNames{fileNr}, numel(dataTemp), ...
            mvt.laneSidecarName(segmentNames{fileNr}), numel(lanes_od), ...
            processingDay);
    end

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
    clear dataTemp lanes_od lane_changes_data
    reportProgress(fileNr, segmentNames{fileNr});
end
reportProgress();

% Combine the cell arrays into two distinct struct arrays
all_lane_changes_start = vertcat(all_start_cell{:});
all_lane_changes_end = vertcat(all_end_cell{:});
fprintf('Writing %s (%d merge-in, %d merge-out events) ...', ...
    outputFile, numel(all_lane_changes_start), numel(all_lane_changes_end)); tic
mvt.atomicSave(outputFile, struct( ...
    'all_lane_changes_end', all_lane_changes_end, ...
    'all_lane_changes_start', all_lane_changes_start));
fprintf(' Done (%0.0fsec).\n',toc)

end