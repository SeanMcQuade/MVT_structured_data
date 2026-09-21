function [] = relative_speed_histogram(processingDay, varargin)
% RELATIVE_SPEED_HISTOGRAM  Relative speed of traffic to the nearest engaged AV.
% (C) 2026/05/26 by Benjamin Seibold, added to by Sean McQuade
%
% Purpose
%   Builds the paper's relative-speed figures for one day: a histogram of the
%   relative speed between each vehicle and the engaged AV ahead of it, and
%   (via binned_relative_speed) the same statistic binned by distance behind
%   that AV.
%
% Inputs
%   processingDay  16, 17, or 18 (November 2022)
%   varargin       options struct and/or name/value pairs (see mvt.options);
%                  Force, Clean, DryRun, Verbose
%
% Outputs
%   <results>/figures/
%     Relative speed histogram, day DD for files j = 1 to 24.pdf
%     Relative speeds behind AV, day DD.pdf
%   The names, spaces included, are the ones the paper already cites, and the
%   figures sit at the figures root rather than in a day folder, which is where
%   this script has always put them. Both are declared in mvt.expectedOutputs.
%
% Algorithm
%   1. Rebuild only when the figures are older than the day's slim segments or
%      than the code (mvt.isStale).
%   2. Per segment, difference the distance to the downstream engaged AV to get
%      a relative speed, using one-sided stencils where the tracked AV changes
%      or the distance jumps, then smooth it.
%   3. Pool samples between lower_Bnd and upper_Bnd metres over the whole day,
%      plot the histogram with mean/median/spread, and hand the pooled samples
%      to binned_relative_speed for the distance-binned figure.
%
% Notes
%   Needs the Statistics and Machine Learning Toolbox (prctile, adtest).
%   Figures are rendered off-screen and closed, so this is safe under
%   `matlab -batch`.
%
% Dependencies
%   mvt.options, mvt.dayDir, mvt.manifest, mvt.expectedOutputs, mvt.isStale,
%   mvt.sources, mvt.ensureDir, mvt.progress, mvt.log, binned_relative_speed
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
% Figure canvas in pixels. Pinned for the same reason as in
% binned_relative_speed: exportgraphics crops to content, and the content
% extent follows a figure size MATLAB otherwise takes from the screen.
figRes = [700 420];

%========================================================================
% Initialize
%========================================================================
outputs = mvt.expectedOutputs('relspeed', processingDay, opts);
histogramFile = outputs{1};
binnedFile = outputs{2};
mvt.ensureDir(fileparts(binnedFile))

segments = mvt.manifest(processingDay, opts);
slimPath = mvt.dayDir('slim', processingDay);
slimFiles = cell(1, numel(segments));
for iSeg = 1:numel(segments)
    slimFiles{iSeg} = fullfile(slimPath, segments(iSeg).outputName);
end

[stale, staleReason] = mvt.isStale(outputs, slimFiles, ...
    mvt.sources('relative_speed_histogram', opts), opts);
if ~stale
    mvt.log(opts, 'skip relative speed figures for 2022-11-%d: %s', ...
        processingDay, staleReason);
    return
end
mvt.log(opts, 'build relative speed figures for 2022-11-%d: %s', ...
    processingDay, staleReason);
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

number_of_data_all_files = length(filtered_speed_all_files);
mean_all_fil_rel_speed = mean(filtered_speed_all_files,'omitnan');
median_fil_all_files = prctile(filtered_speed_all_files,50);
stddev_fil_all_rel_speed = std(filtered_speed_all_files, 'omitnan');
first_quartile_fil_all = prctile(filtered_speed_all_files, 25);
third_quartile_fil_all = prctile(filtered_speed_all_files,75);

%plot relative speed histogram for all files
fig = figure('Visible', 'off', 'Position', [10 50 figRes], ...
    'PaperPositionMode', 'auto');
closeFigure = onCleanup(@() close(fig));
hold on
H = histogram(filtered_speed_all_files,edges);
plot(mean_all_fil_rel_speed*ones(2,1), [0,max(H.Values)],"LineWidth",2)
plot(median_fil_all_files*ones(2,1),[0,0.5*max(H.Values)],"LineWidth",2)
std_dev = [mean_all_fil_rel_speed-stddev_fil_all_rel_speed,...
           mean_all_fil_rel_speed + stddev_fil_all_rel_speed];
interquartile_all = [first_quartile_fil_all, third_quartile_fil_all];
plot(std_dev, [10,10],"LineWidth",3)
plot(interquartile_all, [0.5,0.5],"LineWidth",4)
% fontsz = 24;
fontsz=8;
legendformatSpec_mean = "mean speed = %3.3f (m/s)";
mean_legend = sprintf(legendformatSpec_mean,mean_all_fil_rel_speed);
legendformatSpec_med = "median speed = %3.3f (m/s)";
med_legend = sprintf(legendformatSpec_med,median_fil_all_files);
legend("histogram", mean_legend, med_legend, "Standard Dev", "Interquartile", ...
                                                       "FontSize",fontsz)

if day == 16 %Write the day in the title
    % formatSpec = "Relative speed histogram Wed, " + ...
    %     "for files %d to %d, total data points = %d";
    % title_string = sprintf(formatSpec,j_start,j_end,number_of_data_all_files);
    title_string = "Relative speed histogram Wednesday 16-Nov-2022";
elseif day == 17
    % formatSpec = "Relative speed histogram Thurs, " + ...
    %     "for files %d to %d, total data points = %d";
    % title_string = sprintf(formatSpec,j_start,j_end,number_of_data_all_files);
    title_string = "Relative speed histogram Thursday 17-Nov-2022";
elseif day == 18
    % formatSpec = "Relative speed histogram Fri, " + ...
    %     "for file %d to %d, total data points = %d";
    % title_string = sprintf(formatSpec, j_start, j_end,number_of_data_all_files);
    title_string = "Relative speed histogram Friday 18-Nov-2022";
end

title(title_string,"FontSize",fontsz);
xlim([-15 15]);
xlabel("Relative Speed (m/s)","FontSize",fontsz);
% ylabel_formatSpec = "Frequency of speeds";
% ylabel_string = sprintf(ylabel_formatSpec, lower_Bnd, upper_Bnd);
ylabel_string = "Number of samples";
ylabel(ylabel_string,"FontSize",fontsz);
fprintf('Save figure in %s ...', histogramFile), tic
exportgraphics(fig, histogramFile, 'ContentType', 'vector');
fprintf(' Done (%0.0fsec).\n', toc)

% average speed in standard distance bins from the AV; once a script sharing
% this workspace, now a function taking what it needs
binned_relative_speed(filtered_dist_all_files, filtered_speed_all_files, ...
    day, binnedFile);

% Returns h = 1 if the data is NOT Gaussian, h = 0 if it IS Gaussian
[h,p] = adtest(filtered_speed_all_files);
fprintf('Anderson-Darling on the pooled relative speeds: h = %d, p = %g\n', h, p);
toc
end