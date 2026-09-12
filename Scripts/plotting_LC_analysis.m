clear 
processingDay = 16;

[parentDirectory, ~, ~] = fileparts(pwd);
[dataRootDirectory, ~, ~] = fileparts(parentDirectory);
outputPath = fullfile(dataRootDirectory, 'results', 'slim', ...
    ['2022-11-', num2str(processingDay)]);
filenameLoad = fullfile(outputPath, ['LC_data_' char(num2str(processingDay)) '.mat']);
load(filenameLoad,'all_lane_changes_end','all_lane_changes_start')
dataGPS = jsondecode(fileread(fullfile(dataRootDirectory,...
    'results','gps',['CIRCLES_GPS_10Hz_2022-11-' num2str(processingDay) '.json'])));
dataGPS = dataGPS([dataGPS.direction]<0);


max_T = posixtime(datetime(['2022-11-' char(num2str(processingDay)) ' 11:00:00'],...
    'InputFormat', 'yyyy-MM-dd HH:mm:ss', 'TimeZone', 'America/Chicago'));
dataGPS = dataGPS([dataGPS.first_timestamp] < max_T);

lanesToAnalyze = [2, 3, 4];
binEdges = -500:10:500;
binCenters = binEdges(1:end-1) + diff(binEdges)/2;

distFieldsStart = {'dist_to_eng_av_at_start', 'dist_to_av_at_start'};
distFieldsEnd   = {'dist_to_eng_av_at_end', 'dist_to_av_at_end'};
titleSuffix     = {'Engaged AVs', 'All AVs'};
xLabels         = {'Distance to Nearest Engaged AV (m)', 'Distance to Nearest AV (m)'};

% =======================================================================
% SECTION A: RAW EVENT COUNTS (PARTS 1 - 3)
% =======================================================================

% Define data for raw counts (to toggle between isolated or all merges here)
dataStart = all_lane_changes_start; 
dataEnd   = all_lane_changes_end;   

% -----------------------------------------------------------------------
% PART 1: MERGE OUT (END OF TRAJECTORY) - COUNTS
% -----------------------------------------------------------------------
for caseIdx = 1:2
    figure( 'Color', 'w', 'Name', ['Merge Out: ', titleSuffix{caseIdx}]);
    t = tiledlayout(length(lanesToAnalyze), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(t, ['Merge Out (Leaving Lane) Relative to ', titleSuffix{caseIdx}], 'FontSize', 14, 'Color', 'k');
    
    for i = 1:length(lanesToAnalyze)
        currentLane = lanesToAnalyze(i);
        laneMask = [dataEnd.lane_number] == currentLane;
        
        rawDist = {dataEnd(laneMask).(distFieldsEnd{caseIdx})};
        rawDist(cellfun(@isempty, rawDist)) = {NaN};
        laneDistData = cell2mat(rawDist);
        
        laneChangeData = [dataEnd(laneMask).lane_change_at_end];
        
        distPlus1  = laneDistData(laneChangeData > 0);
        distMinus1 = laneDistData(laneChangeData < 0);
        
        distPlus1  = distPlus1(~isnan(distPlus1));
        distMinus1 = distMinus1(~isnan(distMinus1));
        
        nexttile;
        histogram(distPlus1, binEdges, 'FaceColor', [0 0.4470 0.7410], 'FaceAlpha', 0.6, 'EdgeColor', 'none'); 
        hold on;
        histogram(distMinus1, binEdges, 'FaceColor', [0.8500 0.3250 0.0980], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        hold off;
        
        ylabel('Count'); title(sprintf('Lane %d', currentLane), 'Color', 'k'); xlim([-500 500]); grid on;
        set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
        if i == 1, legend('Change to the right', 'Change to the left', 'Location', 'best'); end
    end
    xlabel(t, xLabels{caseIdx}, 'FontSize', 12, 'Color', 'k');
end

% -----------------------------------------------------------------------
% PART 2: MERGE IN (START OF TRAJECTORY) - COUNTS
% -----------------------------------------------------------------------
for caseIdx = 1:2
    figure('Color', 'w', 'Name', ['Merge In: ', titleSuffix{caseIdx}]);
    t = tiledlayout(length(lanesToAnalyze), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(t, ['Merge In (Entering Lane) Relative to ', titleSuffix{caseIdx}], 'FontSize', 14, 'Color', 'k');
    
    for i = 1:length(lanesToAnalyze)
        currentLane = lanesToAnalyze(i);
        laneMask = [dataStart.lane_number] == currentLane;
        
        rawDist = {dataStart(laneMask).(distFieldsStart{caseIdx})};
        rawDist(cellfun(@isempty, rawDist)) = {NaN};
        laneDistData = cell2mat(rawDist);
        
        laneChangeData = [dataStart(laneMask).lane_change_at_start];
        
        distPlus1  = laneDistData(laneChangeData > 0);
        distMinus1 = laneDistData(laneChangeData < 0);
        
        distPlus1  = distPlus1(~isnan(distPlus1));
        distMinus1 = distMinus1(~isnan(distMinus1));
        
        nexttile;
        histogram(distPlus1, binEdges, 'FaceColor', [0 0.4470 0.7410], 'FaceAlpha', 0.6, 'EdgeColor', 'none'); 
        hold on;
        histogram(distMinus1, binEdges, 'FaceColor', [0.8500 0.3250 0.0980], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        hold off;
        
        ylabel('Count'); title(sprintf('Lane %d', currentLane), 'Color', 'k'); xlim([-500 500]); grid on;
        set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
        if i == 1, legend('Change from the left', 'Change from the right', 'Location', 'best'); end
    end
    xlabel(t, xLabels{caseIdx}, 'FontSize', 12, 'Color', 'k');
end

% -----------------------------------------------------------------------
% PART 3: MERGE IN vs MERGE OUT - COUNTS
% -----------------------------------------------------------------------
for caseIdx = 1:2
    figure( 'Color', 'w', 'Name', ['In vs Out: ', titleSuffix{caseIdx}]);
    t = tiledlayout(length(lanesToAnalyze), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(t, ['Merge In vs Merge Out Relative to ', titleSuffix{caseIdx}], 'FontSize', 14, 'Color', 'k');
    
    for i = 1:length(lanesToAnalyze)
        currentLane = lanesToAnalyze(i);
        
        maskStart = [dataStart.lane_number] == currentLane;
        rawDistIn = {dataStart(maskStart).(distFieldsStart{caseIdx})};
        rawDistIn(cellfun(@isempty, rawDistIn)) = {NaN};
        distIn = cell2mat(rawDistIn);
        distIn = distIn(~isnan(distIn));
        
        maskEnd = [dataEnd.lane_number] == currentLane;
        rawDistOut = {dataEnd(maskEnd).(distFieldsEnd{caseIdx})};
        rawDistOut(cellfun(@isempty, rawDistOut)) = {NaN};
        distOut = cell2mat(rawDistOut);
        distOut = distOut(~isnan(distOut));
        
        nexttile;
        histogram(distIn, binEdges, 'FaceColor', [0.4660 0.6740 0.1880], 'FaceAlpha', 0.6, 'EdgeColor', 'none'); 
        hold on;
        histogram(distOut, binEdges, 'FaceColor', [0.4940 0.1840 0.5560], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        hold off;
        
        ylabel('Count'); title(sprintf('Lane %d', currentLane), 'Color', 'k'); xlim([-500 500]); grid on;
        set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
        if i == 1, legend('Merge In (Start)', 'Merge Out (End)', 'Location', 'best'); end
    end
    xlabel(t, xLabels{caseIdx}, 'FontSize', 12, 'Color', 'k');
end


% Define data for rates (toggle between isolated or all merges here)
dataStart = all_lane_changes_start; 
dataEnd   = all_lane_changes_end;   

dt_seconds = mode(diff(dataGPS(1).timestamp));
dt_minutes = dt_seconds / 60;
% -----------------------------------------------------------------------
% 1. INITIALIZE ALL FIGURES UPFRONT
% -----------------------------------------------------------------------
figExp = gobjects(2,1); tExp = gobjects(2,1);
figOut = gobjects(2,1); tOut = gobjects(2,1);
figIn  = gobjects(2,1); tIn  = gobjects(2,1);
figComp= gobjects(2,1); tComp= gobjects(2,1);

for caseIdx = 1:2
    figExp(caseIdx) = figure('Color', 'w', 'Name', ['Exposure: ', titleSuffix{caseIdx}]);
    tExp(caseIdx) = tiledlayout(length(lanesToAnalyze), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tExp(caseIdx), ['Total Exposure Time: ', titleSuffix{caseIdx}], 'FontSize', 14, 'Color', 'k');

    figOut(caseIdx) = figure('Color', 'w', 'Name', ['Rate Merge Out: ', titleSuffix{caseIdx}]);
    tOut(caseIdx) = tiledlayout(length(lanesToAnalyze), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tOut(caseIdx), ['Rate: Merge Out Relative to ', titleSuffix{caseIdx}], 'FontSize', 14, 'Color', 'k');
    
    figIn(caseIdx) = figure('Color', 'w', 'Name', ['Rate Merge In: ', titleSuffix{caseIdx}]);
    tIn(caseIdx) = tiledlayout(length(lanesToAnalyze), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tIn(caseIdx), ['Rate: Merge In Relative to ', titleSuffix{caseIdx}], 'FontSize', 14, 'Color', 'k');
    
    figComp(caseIdx) = figure('Color', 'w', 'Name', ['Rate In vs Out: ', titleSuffix{caseIdx}]);
    tComp(caseIdx) = tiledlayout(length(lanesToAnalyze), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tComp(caseIdx), ['Rate: Merge In vs Merge Out Relative to ', titleSuffix{caseIdx}], 'FontSize', 14, 'Color', 'k');
end

% -----------------------------------------------------------------------
% 2. LANE LOOP (Calculate Exposure ONCE per lane, Plot BOTH cases)
% -----------------------------------------------------------------------
for i = 1:length(lanesToAnalyze)
    currentLane = lanesToAnalyze(i);
    
    % --- MASTER EXPOSURE CALCULATION ---
    maskStart = [dataStart.lane_number] == currentLane;
    maskEnd   = [dataEnd.lane_number] == currentLane;
    
    all_merge_X = [[dataStart(maskStart).x_position_at_start], [dataEnd(maskEnd).x_position_at_end]];
    if isempty(all_merge_X), continue; end
    X_min = min(all_merge_X);
    X_max = max(all_merge_X);
    
    primaryLanes = arrayfun(@(x) mode(x.assigned_lane), dataGPS);
    laneAVs = dataGPS(primaryLanes == currentLane);
    
    grid_res = 1; 
    spatial_grid = floor(X_min):grid_res:ceil(X_max);
    num_grid_points = length(spatial_grid);
    points_per_bin = diff(binEdges(1:2)) / grid_res;
    
    % Initialize parallel exposure arrays
    exposure_min_all = zeros(length(binCenters), 1);
    exposure_min_eng = zeros(length(binCenters), 1);
    
    all_t = []; all_x = []; all_engaged = [];
    for j = 1:length(laneAVs)
        all_t = [all_t; laneAVs(j).timestamp(:)];
        all_x = [all_x; laneAVs(j).x_position(:)];
        all_engaged = [all_engaged; laneAVs(j).control_car(:)];
    end
    
    if ~isempty(all_t)
        [unique_t, ~, idx_t] = unique(all_t);
        total_k = length(unique_t);
        
        wb = waitbar(0, sprintf('Lane %d: Calculating Parallel Exposure (0%%)', currentLane), ...
            'Name', 'Directional Grid Distances');
            
        for k = 1:total_k
            % Extract positions for this frame
            current_x_all = all_x(idx_t == k);
            current_engaged = all_engaged(idx_t == k);
            current_x_eng = current_x_all(current_engaged == 1);
            
            % -- PROCESS ALL AVs --
            if ~isempty(current_x_all)
                dist_matrix = spatial_grid - current_x_all ; 
                if isscalar(current_x_all)
                    all_nearest_dists = dist_matrix;
                else
                    pos_dist = dist_matrix; pos_dist(pos_dist <= 0) = NaN;
                    min_pos = min(pos_dist, [], 1);
                    
                    neg_dist = dist_matrix; neg_dist(neg_dist > 0) = NaN;
                    max_neg = max(neg_dist, [], 1); 
                    
                    all_nearest_dists = [min_pos(~isnan(min_pos)), max_neg(~isnan(max_neg))];
                end
                counts = histcounts(all_nearest_dists, binEdges);
                exposure_min_all = exposure_min_all + (counts(:) ./ points_per_bin) .* dt_minutes;
            end
            
            % -- PROCESS ENGAGED AVs --
            if ~isempty(current_x_eng)
                dist_matrix = spatial_grid - current_x_eng ; 
                if isscalar(current_x_eng)
                    all_nearest_dists = dist_matrix;
                else
                    pos_dist = dist_matrix; pos_dist(pos_dist <= 0) = NaN;
                    min_pos = min(pos_dist, [], 1);
                    
                    neg_dist = dist_matrix; neg_dist(neg_dist > 0) = NaN;
                    max_neg = max(neg_dist, [], 1); 
                    
                    all_nearest_dists = [min_pos(~isnan(min_pos)), max_neg(~isnan(max_neg))];
                end
                counts = histcounts(all_nearest_dists, binEdges);
                exposure_min_eng = exposure_min_eng + (counts(:) ./ points_per_bin) .* dt_minutes;
            end
            
            if mod(k, 100) == 0 || k == total_k
                pct = k / total_k;
                waitbar(pct, wb, sprintf('Lane %d: Calculating Parallel Exposure (%d%%)', currentLane, round(pct * 100)));
            end
        end
        close(wb);
    end
    
    exposure_min_all(exposure_min_all == 0) = NaN;
    exposure_min_eng(exposure_min_eng == 0) = NaN;


    % --- 3. PLOT FOR BOTH CASES ---
    for caseIdx = 1:2
        % Select the correctly matched exposure array for the loop iteration
        if caseIdx == 1
            exposure_minutes = exposure_min_eng;
        else
            exposure_minutes = exposure_min_all;
        end
        % Plot Exposure time 
        set(0, 'CurrentFigure', figExp);
        nexttile(tExp(caseIdx)); hold on;
        plot(binCenters, exposure_minutes, 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
        fill([binCenters, fliplr(binCenters)], [exposure_minutes', zeros(1, length(binCenters))], ...
            [0.8500 0.3250 0.0980], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        hold off;
        ylabel('Exposure (Min)'); title(sprintf('Lane %d', currentLane), 'Color', 'k'); grid on; xlim([-500 500]);
        set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
        if i == length(lanesToAnalyze), xlabel(tExp(caseIdx), xLabels{caseIdx}, 'FontSize', 12, 'Color', 'k'); end

        % -------------------------------------------------------------------
        % APPLY TO PART 4 (MERGE OUT)
        % -------------------------------------------------------------------
        rawDistOut = {dataEnd(maskEnd).(distFieldsEnd{caseIdx})};
        rawDistOut(cellfun(@isempty, rawDistOut)) = {NaN};
        laneDistDataOut = cell2mat(rawDistOut);
        laneChangeDataOut = [dataEnd(maskEnd).lane_change_at_end];
        
        distPlus1Out  = laneDistDataOut(laneChangeDataOut > 0);
        distMinus1Out = laneDistDataOut(laneChangeDataOut < 0);
        
        ratePlus1Out  = histcounts(distPlus1Out, binEdges)' ./ exposure_minutes;
        rateMinus1Out = histcounts(distMinus1Out, binEdges)' ./ exposure_minutes;
        
        set(0, 'CurrentFigure', figOut(caseIdx));
        nexttile(tOut(caseIdx)); hold on;
        bar(binCenters, ratePlus1Out, 1, 'FaceColor', [0 0.4470 0.7410], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        bar(binCenters, rateMinus1Out, 1, 'FaceColor', [0.8500 0.3250 0.0980], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        hold off;
        ylabel('Merges / Min'); title(sprintf('Lane %d', currentLane), 'Color', 'k'); xlim([-500 500]); grid on;
        set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
        if i == 1, legend('Change to the right', 'Change to the left', 'Location', 'best'); end
        if i == length(lanesToAnalyze), xlabel(tOut(caseIdx), xLabels{caseIdx}, 'FontSize', 12, 'Color', 'k'); end
        
        % -------------------------------------------------------------------
        % APPLY TO PART 5 (MERGE IN)
        % -------------------------------------------------------------------
        rawDistIn = {dataStart(maskStart).(distFieldsStart{caseIdx})};
        rawDistIn(cellfun(@isempty, rawDistIn)) = {NaN};
        laneDistDataIn = cell2mat(rawDistIn);
        laneChangeDataIn = [dataStart(maskStart).lane_change_at_start];
        
        distPlus1In  = laneDistDataIn(laneChangeDataIn > 0);
        distMinus1In = laneDistDataIn(laneChangeDataIn < 0);
        
        ratePlus1In  = histcounts(distPlus1In, binEdges)' ./ exposure_minutes;
        rateMinus1In = histcounts(distMinus1In, binEdges)' ./ exposure_minutes;
        
        set(0, 'CurrentFigure', figIn(caseIdx));
        nexttile(tIn(caseIdx)); hold on;
        bar(binCenters, ratePlus1In, 1, 'FaceColor', [0 0.4470 0.7410], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        bar(binCenters, rateMinus1In, 1, 'FaceColor', [0.8500 0.3250 0.0980], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        hold off;
        ylabel('Merges / Min'); title(sprintf('Lane %d', currentLane), 'Color', 'k'); xlim([-500 500]); grid on;
        set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
        if i == 1, legend('Change from the left', 'Change from the right', 'Location', 'best'); end
        if i == length(lanesToAnalyze), xlabel(tIn(caseIdx), xLabels{caseIdx}, 'FontSize', 12, 'Color', 'k'); end
        
        % -------------------------------------------------------------------
        % APPLY TO PART 6 (IN VS OUT AGNOSTIC)
        % -------------------------------------------------------------------
        distIn_agnostic = laneDistDataIn(~isnan(laneDistDataIn));
        distOut_agnostic = laneDistDataOut(~isnan(laneDistDataOut));
        
        rateInAgnostic  = histcounts(distIn_agnostic, binEdges)' ./ exposure_minutes;
        rateOutAgnostic = histcounts(distOut_agnostic, binEdges)' ./ exposure_minutes;
        
        set(0, 'CurrentFigure', figComp(caseIdx));
        nexttile(tComp(caseIdx)); hold on;
        bar(binCenters, rateInAgnostic, 1, 'FaceColor', [0.4660 0.6740 0.1880], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        bar(binCenters, rateOutAgnostic, 1, 'FaceColor', [0.4940 0.1840 0.5560], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        hold off;
        ylabel('Merges / Min'); title(sprintf('Lane %d', currentLane), 'Color', 'k'); xlim([-500 500]); grid on;
        set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
        if i == 1, legend('Merge In (Start)', 'Merge Out (End)', 'Location', 'best'); end
        if i == length(lanesToAnalyze), xlabel(tComp(caseIdx), xLabels{caseIdx}, 'FontSize', 12, 'Color', 'k'); end
        
    end
end