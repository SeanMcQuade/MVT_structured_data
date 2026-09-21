function [] = relative_speed_histogram(processingDay, varargin)
% RELATIVE_SPEED_HISTOGRAM  Relative speed of traffic to the nearest engaged AV.
% (C) 2026/05/26 by Benjamin Seibold, added to by Sean McQuade
%
% Purpose
%   Pools, for one day, the relative speed between each vehicle and the engaged
%   AV ahead of it, and saves the result. This is the expensive half of the
%   relative-speed analysis: it decodes all 24 slim segments. plot_relative_speed
%   turns the saved samples into the figures, so those can be remade without the
%   slim tree.
%
% Inputs
%   processingDay  16, 17, or 18 (November 2022)
%   varargin       options struct and/or name/value pairs (see mvt.options);
%                  Force, Clean, DryRun, Verbose
%
% Outputs
%   <results>/analysis/2022-11-DD/relspeed_data_DD.mat, holding
%     filtered_dist_all_files   distance to the nearest engaged AV (m), pooled
%                               over the day's segments
%     filtered_speed_all_files  the matching relative speeds (m/s)
%     plus the bounds and segment range they were pooled under, as provenance
%
% Algorithm
%   1. Rebuild only when the figures are older than the day's slim segments or
%      than the code (mvt.isStale).
%   2. Per segment, difference the distance to the downstream engaged AV to get
%      a relative speed, using one-sided stencils where the tracked AV changes
%      or the distance jumps, then smooth it.
%   3. Pool samples between lower_Bnd and upper_Bnd metres over the whole day
%      and save them for plot_relative_speed.
%
% Notes
%   Needs the Statistics and Machine Learning Toolbox (prctile, adtest).
%
% Dependencies
%   mvt.options, mvt.dayDir, mvt.manifest, mvt.expectedOutputs, mvt.isStale,
%   mvt.sources, mvt.atomicSave, mvt.ensureDir, mvt.progress, mvt.log
%
if nargin < 1
    error(['Specify the day of Nov. 2022 MVT to build the relative speed ' ...
        'figures (from 16 to 18)']);
end
mvt.assertDay(processingDay)
opts = mvt.options(varargin{:});

%========================================================================
% Parameters
%========================================================================
% First and last of the day's 24 segments to pool. These appear in the output
% filename, so mvt.expectedOutputs matches them with a glob.
j_start = 1; j_end = 24;
edges = -30:0.5:30;          % histogram bin edges for the relative speed (m/s)
% Only samples this far from the AV are pooled: 30 m puts the 35-45 m bin first
lower_Bnd = 30;              % (m) nearest relative distance considered
upper_Bnd = 350;             % (m) furthest relative distance considered

%========================================================================
% Initialize
%========================================================================
outputs = mvt.expectedOutputs('relspeed', processingDay, opts);
outputFile = outputs{1};
mvt.ensureDir(fileparts(outputFile))

segments = mvt.manifest(processingDay, opts);
slimPath = mvt.dayDir('slim', processingDay);
slimFiles = cell(1, numel(segments));
for iSeg = 1:numel(segments)
    slimFiles{iSeg} = fullfile(slimPath, segments(iSeg).outputName);
end

[stale, staleReason] = mvt.isStale(outputs, slimFiles, ...
    mvt.sources('relative_speed_histogram', opts), opts);
if ~stale
    mvt.log(opts, 'skip relspeed_data_%d.mat: %s', processingDay, staleReason);
    return
end
mvt.log(opts, 'build relspeed_data_%d.mat: %s', processingDay, staleReason);
if opts.Clean && ~opts.DryRun
    mvt.removeOutputs(outputs, opts);
end
if opts.DryRun
    return
end

tic
day = processingDay;

%preallocate cell array to save rel speed strcat(datafolder,
filtered_dist_all_files = [];
filtered_speed_all_files = [];
all_av_dist = [];

% Segment order comes from the manifest rather than dir(), so j indexes the
% same segment every run regardless of how the folder happens to sort.
reportProgress = mvt.progress(j_end - j_start + 1, ...
    sprintf('relspeed 2022-11-%d', processingDay), 'Opts', opts);
for j=j_start:j_end
        filename = segments(j).outputName;
        fprintf('Loading %s ...',filename), tic
        data = jsondecode(fileread(slimFiles{j}));
        fprintf(' Done (%0.0fsec).\n',toc)

% Process data
av_dist_all_segments = [];
rel_speed_all_segments = [];
for i = 1:length(data) % loop over all segments
    times = data(i).timestamp;
    av_id = data(i).downstream_engaged_av_id;
    av_dist = data(i).distance_to_downstream_engaged_av_meters; 
    
    ind = ~isnan(av_id); % only keep those instances where an AV is present
    if numel(ind)<2, continue, end % if segment too short, move to next segment
    times = times(ind); av_id = av_id(ind); av_dist = av_dist(ind);
    
    ind_av_switch = diff(av_id)~=0; % AV index switching
    ind_large_jump = abs(diff(av_dist)./diff(times))>70; % excessive jumps
    ind_switch = find(ind_av_switch|ind_large_jump); % indices of switching
    
    n_times = length(times);
    % finite difference stencil, usually centered, but one-sided at switchings
    stencil_forward = [2:n_times,n_times];
    stencil_forward(ind_switch) = stencil_forward(ind_switch)-1;
    stencil_backward = [1,1:n_times-1];
    stencil_backward(ind_switch+1) = stencil_backward(ind_switch+1)+1;
    
    rel_speed_raw = (av_dist(stencil_forward)-av_dist(stencil_backward))./...
        (times(stencil_forward)-times(stencil_backward)); % finite differences
    
    rel_speed = rel_speed_raw;
    for k = 1:40 % apply local averaging to get smooth relative speed
        rel_speed = (rel_speed(stencil_backward)+...
            rel_speed+rel_speed(stencil_forward))/3;
    end

    % save rel_speed and av_dist for all segments
    av_dist_all_segments = [av_dist_all_segments; av_dist];
    rel_speed_all_segments = [rel_speed_all_segments; rel_speed];
end

%only keep
clear filtered_ind filtered_dist filtered_speed;
filtered_ind_low = find(lower_Bnd < av_dist_all_segments);
filtered_ind_up = find(av_dist_all_segments < upper_Bnd);
filtered_ind = intersect(filtered_ind_low,filtered_ind_up);
filtered_dist = av_dist_all_segments(filtered_ind);
filtered_speed = rel_speed_all_segments(filtered_ind);

number_of_data = length(filtered_dist);
mean_rel_speed(j) = mean(filtered_speed,'omitnan');
median_rel_speed(j) = prctile(filtered_speed, 50);
stddev_rel_speed(j) = std(filtered_speed, 'omitnan');
first_quartile(j) = prctile(filtered_speed, 25);
third_quartile(j) = prctile(filtered_speed,75);
interquartile = [first_quartile(j), third_quartile(j)];

%save rel distance and speed from each iteration
filtered_dist_all_files = [filtered_dist_all_files; filtered_dist];
filtered_speed_all_files = [filtered_speed_all_files; filtered_speed];
reportProgress(j - j_start + 1, filename);

% %plot relative speed histogram for jth file
% figure;
% hold on
% H = histogram(filtered_speed,edges);
% plot(mean_rel_speed(j)*ones(2,1), [0,max(H.Values)],"LineWidth",2)
% plot(median_rel_speed(j)*ones(2,1),[0,0.5*max(H.Values)],"LineWidth",2)
% std_dev = [mean_rel_speed(j)-stddev_rel_speed(j),...
%     mean_rel_speed(j) + stddev_rel_speed(j)];
% plot(std_dev, [0.5,0.5],"LineWidth",3)
% plot(interquartile, [0.5,0.5],"LineWidth",4)
% fontsz = 24;
% legend("Histogram", "Mean speed", "Median Speed", ...
%                   "Standard Dev", "Interquartile","FontSize",fontsz)
% xlabel("Relative Speed m/s","FontSize",fontsz);
% ylabel_formatSpec = "Frequency of speeds between %d and %d m/s";
% ylabel_string = sprintf(ylabel_formatSpec, lower_Bnd, upper_Bnd);
% ylabel(ylabel_string,"FontSize",fontsz);
% 
% 
% if day == 16 %Write the day in the title
%     formatSpec = "Relative speed histogram Wed," + ...
%         "for file %d, total data points = %d";
%     title_string = sprintf(formatSpec,j,number_of_data);
% elseif day == 17
%     formatSpec = "Relative speed histogram Thurs, " + ...
%         "for file %d, total data points = %d";
%     title_string = sprintf(formatSpec,j,number_of_data);
% elseif day == 18
%     formatSpec = "Relative speed histogram Fri," + ...
%         "for file %d, total data points = %d";
%     title_string = sprintf(formatSpec,j,number_of_data);
% end
% 
% title(title_string,"FontSize",fontsz);
% formatSpecSave= "../../results/figures/Relative speed histogram, day %d for file j = %d.png";
% savename = sprintf(formatSpecSave,day, j);
% saveas(gcf,savename)
end
reportProgress();

%========================================================================
% Save the pooled samples
%========================================================================
% The bounds and segment range travel with the data: the figures are only
% comparable across days if they were pooled the same way.
fprintf('Writing %s (%d samples) ...', outputFile, ...
    numel(filtered_speed_all_files)); tic
mvt.atomicSave(outputFile, struct( ...
    'filtered_dist_all_files', filtered_dist_all_files, ...
    'filtered_speed_all_files', filtered_speed_all_files, ...
    'lower_Bnd', lower_Bnd, 'upper_Bnd', upper_Bnd, ...
    'j_start', j_start, 'j_end', j_end, 'day', processingDay));
fprintf(' Done (%0.0fsec).\n', toc)

% Returns h = 1 if the data is NOT Gaussian, h = 0 if it IS Gaussian
[h,p] = adtest(filtered_speed_all_files);
fprintf('Anderson-Darling on the pooled relative speeds: h = %d, p = %g\n', h, p);
toc
end
