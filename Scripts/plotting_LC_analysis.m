function [] = plotting_LC_analysis(processingDay)
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
titleSuffix     = {'Engaged AVs'};
titleDay        = {'Wednesday 16-Nov-2022' , 'Thursday 17-Nov-2022', 'Friday 18-Nov-2022'};
xLabels         = {'Distance to Nearest Engaged AV (m)'};

twoColWidth  = 7.5*3; 
twoColHeight = 2.8*4;

oneColWidth  = 3.5*4; 
oneColHeight = 2.6*4;

% Define data for rates (toggle between isolated or all merges here)
dataStart = all_lane_changes_start; 
dataEnd   = all_lane_changes_end;   

dt_seconds = mode(diff(dataGPS(1).timestamp));
dt_minutes = dt_seconds / 60;

% INITIALIZE ALL FIGURES UPFRONT


caseIdx = 1;
figExp = figure('Color', 'w', 'Name', ['Exposure: ', titleSuffix{caseIdx}], ...
    'Units', 'inches', 'Position', [1, 1, twoColWidth, twoColHeight]);
tExp = tiledlayout(length(lanesToAnalyze), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tExp, ['Total Exposure Time Relative to ', titleSuffix{caseIdx},' on ',...
    titleDay{processingDay-15}], 'FontSize', 20, 'Color', 'k', 'FontWeight', 'bold');

figComp = figure('Color', 'w', 'Name', ['Rate Out: ', titleSuffix{caseIdx}], ...
    'Units', 'inches', 'Position', [1, 1, twoColWidth, twoColHeight]);
tComp = tiledlayout(length(lanesToAnalyze), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tComp, ['Rate of Merge-Out Relative to ', titleSuffix{caseIdx}, ' on ',...
    titleDay{processingDay-15}], 'FontSize', 20, 'Color', 'k', 'FontWeight', 'bold');

figIntegral = figure('Color', 'w', 'Name', 'Cumulative Excess Integral', ...
    'Units', 'inches', 'Position', [1, 1, twoColWidth, twoColHeight]);
tIntegral = tiledlayout(length(lanesToAnalyze), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tIntegral, ['Cumulative Excess Merge-Out Relative to ', titleSuffix{caseIdx}, ' on ', ...
    titleDay{processingDay-15}], 'FontSize', 20, 'Color', 'k', 'FontWeight', 'bold');

figCombinedG = figure('Color', 'w', 'Name', 'Combined Cumulative Excess g(x)');
figCombinedG.Position = [100, 100, 800, 300]*2; 
axCombined = axes('Parent', figCombinedG);
hold(axCombined, 'on');
colorOrder = [
    0.000, 0.447, 0.741;  
    0.850, 0.325, 0.098;  
    0.494, 0.184, 0.556   
    ];
h_combined_lines = gobjects(length(lanesToAnalyze), 1);

% LANE LOOP (Calculate Exposure ONCE per lane, Plot BOTH cases)
for i = 1:length(lanesToAnalyze)
    currentLane = lanesToAnalyze(i);

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

            % -- PROCESS ALL AVs -- To give the option for further analysis 
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
    

        % Use engaged exposure time
        exposure_minutes = exposure_min_eng;
   
        % Plot Exposure time 
        set(0, 'CurrentFigure', figExp);
        nexttile(tExp); hold on;
        plot(binCenters, exposure_minutes, 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
        fill([binCenters, fliplr(binCenters)], [exposure_minutes', zeros(1, length(binCenters))], ...
            [0.8500 0.3250 0.0980], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        hold off;
        ylabel('Exposure (Min)', 'FontSize', 18, 'FontWeight', 'bold'); title(sprintf('Lane %d', currentLane), 'FontSize', 20, 'Color', 'k', 'FontWeight', 'bold'); grid on; xlim([-500 500]);
        set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontSize', 16, 'FontWeight', 'bold');
        if i == length(lanesToAnalyze), xlabel(tExp, xLabels{caseIdx}, 'FontSize', 20, 'Color', 'k', 'FontWeight', 'bold'); end

        rawDistOut = {dataEnd(maskEnd).(distFieldsEnd{caseIdx})};
        rawDistOut(cellfun(@isempty, rawDistOut)) = {NaN};
        laneDistDataOut = cell2mat(rawDistOut);
        laneChangeDataOut = [dataEnd(maskEnd).lane_change_at_end];
        
        distPlus1Out  = laneDistDataOut(laneChangeDataOut > 0);
        distMinus1Out = laneDistDataOut(laneChangeDataOut < 0);
        
        ratePlus1Out  = histcounts(distPlus1Out, binEdges)' ./ exposure_minutes;
        rateMinus1Out = histcounts(distMinus1Out, binEdges)' ./ exposure_minutes;
        rateTotalOut  = ratePlus1Out + rateMinus1Out;
        
        % Calculate Base Level (average from 350m to 500m)
        base_mask = binCenters >= 350 & binCenters <= 500;
        base_level = mean(rateTotalOut(base_mask), 'omitnan');
        % Split data into left (downstream) and right (upstream) to break the line at 0
        left_mask  = binCenters < 0;
        right_mask = binCenters > 0;

        x_left = binCenters(left_mask);
        y_left = rateTotalOut(left_mask);

        x_right = binCenters(right_mask);
        y_right = rateTotalOut(right_mask);
        % --- Indexing for the Faded Exclusion Zone ---
        idx_faded_left = (length(x_left)-3) : length(x_left); 
        idx_valid_left = 1 : (length(x_left)-3);

        idx_faded_right = 1 : 4; 
        idx_valid_right = 4 : length(x_right);

        % Set Colors
        solidColor = [0.4940 0.1840 0.5560]; % Purple
        fadedColor = solidColor * 0.3 + [1 1 1] * 0.7; % Faded pastel version

        % Set up plot
        set(0, 'CurrentFigure', figComp);
        nexttile(tComp); hold on;

        % 1. Plot the Baseline Average
        yline(base_level, '--', sprintf('Basal rate (350-500m): %.2f Merges / Min', base_level), ...
            'Color', 0.2*[1 1 1], 'LineWidth', 1.5, ...
            'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment', 'bottom','FontSize', 18, 'FontWeight', 'bold');

        % 2. Plot the faded segments near the AV
        plot(x_left(idx_faded_left), y_left(idx_faded_left), '-o', ...
            'Color', fadedColor, 'MarkerFaceColor', fadedColor, 'MarkerEdgeColor', fadedColor, 'LineWidth', 1.5);
        plot(x_right(idx_faded_right), y_right(idx_faded_right), '-o', ...
            'Color', fadedColor, 'MarkerFaceColor', fadedColor, 'MarkerEdgeColor', fadedColor, 'LineWidth', 1.5);

        % 3. Plot the valid solid segments
        plot(x_left(idx_valid_left), y_left(idx_valid_left), '-o', ...
            'Color', solidColor, 'MarkerFaceColor', solidColor, 'MarkerEdgeColor', solidColor, 'LineWidth', 2);
        plot(x_right(idx_valid_right), y_right(idx_valid_right), '-o', ...
            'Color', solidColor, 'MarkerFaceColor', solidColor, 'MarkerEdgeColor', solidColor, 'LineWidth', 2);

        hold off;
        ylabel('Merge-Out Rate (Merges/Min)', 'FontSize', 18, 'FontWeight', 'bold'); title(sprintf('Lane %d', currentLane), 'FontSize', 20, 'Color', 'k'); xlim([-500 500]); grid on;
        set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontSize', 16, 'FontWeight', 'bold');

      
        hold on
        h_solid = plot(nan, nan, '-o', 'Color', solidColor, 'MarkerFaceColor', solidColor, 'LineWidth', 2);
        h_faded = plot(nan, nan, '-o', 'Color', fadedColor, 'MarkerFaceColor', fadedColor, 'LineWidth', 1.5);
        legend([h_solid, h_faded], {'Merge Rate', 'Masked (<30m)'}, 'Location', 'best','FontSize', 16); 

        ylim([0.2,.7])
        if i == length(lanesToAnalyze), xlabel(tComp, xLabels{caseIdx}, 'FontSize', 20, 'Color', 'k'); end
        
       
        
        % SPLIT DOMAIN OUTSIDE 30m EXCLUSION ZONE
        mask_down = binCenters <= -30 & binCenters >= -500;
        x_down = binCenters(mask_down);
        y_down = rateTotalOut(mask_down);

        mask_up = binCenters >= 30 & binCenters <= 500;
        x_up = binCenters(mask_up);
        y_up = rateTotalOut(mask_up);

       
        
        % CALCULATE EXCESS RATE e(x)
        r_b_up   = mean(y_up(x_up >= 350 & x_up <= 500), 'omitnan');
        e_up = y_up - r_b_up;
        e_up(isnan(e_up)) = 0;

        % Fit on the first 10 reliable points (35m to 125m) 
        x_missing_up   = [5; 15; 25];
        x_fit_up = x_up(1:10); x_fit_up = x_fit_up(:);
        e_fit_up = e_up(1:10); e_fit_up = e_fit_up(:);
        p_up = polyfit(x_fit_up, e_fit_up, 1); 
        e_missing_up = polyval(p_up, x_missing_up); % Evaluate at 5, 15, 25
        
        x_up = [x_missing_up; x_up(:)];
        e_up = [e_missing_up; e_up(:)];
        
        % CALCULATE CUMULATIVE FUNCTION g(x)
        g_up = cumsum(e_up) ;
        


        % PLOT CUMULATIVE EXCESS g(x)
        solidColor = [0 0.4470 0.7410]; % Blue
        fadedColor = solidColor * 0.3 + [1 1 1] * 0.7; % Faded pastel blue
        set(0, 'CurrentFigure', figIntegral);
        nexttile(tIntegral); hold on;

        x_up_row   = x_up(:)';
        g_up_row   = g_up(:)';

        g_350_up = interp1(x_up_row, g_up_row, 350, 'linear');
        insert_idx_up = find(x_up_row > 350, 1, 'first');
        x_up_row = [x_up_row(1:insert_idx_up-1), 350, x_up_row(insert_idx_up:end)];
        g_up_row = [g_up_row(1:insert_idx_up-1), g_350_up, g_up_row(insert_idx_up:end)];

        idx_up_boundary = find(x_up_row == 350, 1, 'first'); 
        idx_up_solid = 1 : idx_up_boundary;
        idx_up_faded = idx_up_boundary : length(x_up_row);

        plot(x_up_row(idx_up_faded), g_up_row(idx_up_faded), 'LineWidth', 2, 'Color', fadedColor);
        plot(x_up_row(idx_up_solid), g_up_row(idx_up_solid), 'LineWidth', 2, 'Color', solidColor);
        ylim([0 2.5])
        xline(350, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        hold off;
        ylabel('Excess Merges (Merges / Min)', 'FontSize', 18, 'FontWeight', 'bold'); title(sprintf('Lane %d', currentLane), 'FontSize', 20, 'Color', 'k', 'FontWeight', 'bold'); xlim([0 500]); grid on;
        set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontSize', 16, 'FontWeight', 'bold');

        if i == length(lanesToAnalyze)
            xlabel(tIntegral, xLabels{caseIdx}, 'FontSize', 20, 'Color', 'k', 'FontWeight', 'bold');
        end
        c = colorOrder(i, :);
        faded_c = c * 0.3 + [1 1 1] * 0.7;
        plot(axCombined, x_up_row(idx_up_faded), g_up_row(idx_up_faded), 'LineWidth', 2, 'Color', faded_c, 'HandleVisibility', 'off');

        h_combined_lines(i) = plot(axCombined, x_up_row(idx_up_solid), g_up_row(idx_up_solid), 'LineWidth', 2, 'Color', c, 'DisplayName', sprintf('Lane %d', currentLane));
    
end
xline(axCombined, 350, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
set(axCombined, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontSize', 16, 'FontWeight', 'bold');
ylabel(axCombined, 'Excess Merges (Merges / Min)', 'FontSize', 18, 'FontWeight', 'bold'); 
xlabel(axCombined, xLabels{caseIdx}, 'FontSize', 20, 'Color', 'k', 'FontWeight', 'bold');
xlim(axCombined, [0 500]); 
grid(axCombined, 'on');
title(axCombined, ['Cumulative Excess Merge-Out Relative to ', titleSuffix{caseIdx}, ' on ', titleDay{processingDay-15}], 'FontSize', 20, 'Color', 'k', 'FontWeight', 'bold');

legend(axCombined, h_combined_lines, 'Location', 'best','FontSize', 16);

drawnow;
axOriginal = nexttile(tIntegral, 1); 
originalUnits = get(axOriginal, 'Units');
set(axOriginal, 'Units', 'pixels');
origPos = get(axOriginal, 'Position');
subplotWidth = origPos(3);
subplotHeight = origPos(4);
set(axOriginal, 'Units', originalUnits);


set(axCombined, 'Units', 'pixels');
combAxPos = get(axCombined, 'Position');
combAxPos(3) = subplotWidth;
combAxPos(4) = subplotHeight;
set(axCombined, 'Position', combAxPos);

figPos = get(figCombinedG, 'Position');
figPos(3) = subplotWidth + 300; 
figPos(4) = subplotHeight + 240;
set(figCombinedG, 'Position', figPos);


set(axCombined, 'Units', 'normalized');


 linkaxes(findobj(figExp, 'Type', 'axes'), 'y');
end