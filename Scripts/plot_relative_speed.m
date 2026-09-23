function [] = plot_relative_speed(processingDay, varargin)
% PLOT_RELATIVE_SPEED  Relative-speed figures from the pooled samples.
%
% Purpose
%   The plotting half of the relative-speed analysis. relative_speed_histogram
%   does the expensive part - decoding a day of slim trajectories and pooling
%   the relative speeds - and saves the result; this turns that .mat into the
%   two figures, so the figures can be remade without the slim tree.
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
%   The names, spaces included, are the ones the paper already cites.
%
% Notes
%   Needs the Statistics and Machine Learning Toolbox (prctile). Figures are
%   rendered off-screen at an explicit canvas size and closed, so this is safe
%   under `matlab -batch` and reproduces the published page geometry.
%
% Dependencies
%   mvt.options, mvt.dayDir, mvt.expectedOutputs, mvt.isStale, mvt.sources,
%   mvt.ensureDir, mvt.log, binned_relative_speed
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

if nargin < 1
    error(['Specify the day of Nov. 2022 MVT to plot the relative speed ' ...
        'figures (from 16 to 18)']);
end
mvt.assertDay(processingDay)
opts = mvt.options(varargin{:});

%========================================================================
% Parameters
%========================================================================
edges = -30:0.5:30;          % histogram bin edges for the relative speed (m/s)
% Figure canvas in pixels. Pinned because MATLAB otherwise takes the default
% from the screen, and exportgraphics crops to the drawn content: the same code
% then produced a different page size under `matlab -batch` than interactively.
figRes = [700 420];

%========================================================================
% Initialize
%========================================================================
outputs = mvt.expectedOutputs('relspeedplot', processingDay, opts);
% outputs{1} is declared as a glob, because the segment range is part of the
% name and is a tunable of the data stage. The concrete name is built below,
% once the pooled samples have said which segments they cover.
binnedFile = outputs{2};
figuresRoot = fileparts(binnedFile);
mvt.ensureDir(figuresRoot)

dataFile = mvt.expectedOutputs('relspeed', processingDay, opts);
dataFile = dataFile{1};
if ~isfile(dataFile)
    error('mvt:plot_relative_speed:missingData', ...
        'Missing %s. Run `make relspeed-%d` first.', dataFile, processingDay);
end

[stale, staleReason] = mvt.isStale(outputs, dataFile, ...
    mvt.sources('plot_relative_speed', opts), opts);
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

pooled = load(dataFile, 'filtered_dist_all_files', 'filtered_speed_all_files', ...
    'j_start', 'j_end');
filtered_dist_all_files = pooled.filtered_dist_all_files;
filtered_speed_all_files = pooled.filtered_speed_all_files;
j_start = pooled.j_start;
j_end = pooled.j_end;
day = processingDay;

histogramFile = fullfile(figuresRoot, sprintf( ...
    'Relative speed histogram, day %d for files j = %d to %d.pdf', ...
    processingDay, j_start, j_end));

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

end