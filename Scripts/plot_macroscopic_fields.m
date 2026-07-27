function [] = plot_macroscopic_fields(processingDay, varargin)
% PLOT_MACROSCOPIC_FIELDS  Plot the macroscopic traffic fields (Fig. 3, SM5).
%
% Purpose
%   Renders the time-space fields produced by GENERATE_MACROSCOPIC_FIELDS, one
%   figure per field, optionally overlaid with the CIRCLES control-vehicle GPS
%   traces coloured by whether control was engaged.
%
% Inputs
%   processingDay  16, 17, or 18 (November 2022)
%   varargin       options struct and/or name/value pairs (see mvt.options);
%                  Force, Clean, DryRun, Verbose
%   Files:
%     <results>/figures/2022-11-DD/fields_motion_2022-11-DD.mat
%     <results>/gps/CIRCLES_GPS_10Hz_2022-11-DD.json   (when overlaying AVs)
%
% Outputs
%   <results>/figures/2022-11-DD/
%     fig_field_<yyyyMMdd>_<direction>_<lane>_motion_<field>[_av]_nature_large.png
%   for field in {Rho, Q, F, U, Phi, Psi}
%
% Parameters (constants at the top of this file)
%   direction = -1 (westbound), lane = 0 (all lanes), AV overlay on, and the
%   per-field colorbar limits that make the days comparable to one another.
%
% Parallel safety
%   Writes only into this day's figures folder; days may run concurrently.
%   Requires a working graphics stack (see README for headless notes).
%
% Dependencies
%   mvt.options, mvt.paths, mvt.dayDir, mvt.isStale, mvt.sources,
%   mvt.expectedOutputs, mvt.log
%
% (C) 2025-2026 CIRCLES Consortium. Author: Benjamin Seibold, edited by
% Sulaiman Almatrudi. BSD-3-Clause.
if nargin < 1
    error(['Specify the day of Nov. 2022 MVT to plot'...
        ' macroscopic fileds for (from 16 to 18)']);
end
mvt.assertDay(processingDay)
opts = mvt.options(varargin{:});
%========================================================================
% Parameters
%========================================================================
direction = -1; % -1=Westbound, 1=Eastbound
lane = 0; % 0= all lanes
plotFields = {'Rho','Q','F','U','Phi','Psi'};
flagSaveFigures = 1; % if true, save figure into png file
flagPlotAvTrajectories = 1; % if true, overlay fields with AV traces
skipTPlot = 1; % sub-sample trajectories for plotting
subfieldNameX = 'x_position';
timeZoomWindow = [0610 0950]; % time window of intrest to zoom into (in military time)
% Set maximum value limit for plotting ,  set to 0 to defualt to 99th
% percentile of the field vaules to be plotted
colorbarLimit.Rho = 220; % (veh/km)
colorbarLimit.Q =  8000; % (veh/hr)
colorbarLimit.F =  0.25; % (g/s)
colorbarLimit.U =  35; % (m/s)
colorbarLimit.Phi =  4.5; % (g/s)
colorbarLimit.Psi =  0.22; % (g/m)
%========================================================================
% Load data files
%========================================================================
% mvt.paths resolves the layout from the location of the code, so this no
% longer depends on the current folder being Scripts/.
p = mvt.paths();
parentDirectory = p.repoRoot; %#ok<NASGU> % retained for local edits/debugging
dataRootDirectory = p.dataRoot;
inputPath = mvt.dayDir('figures', processingDay);

filename = fullfile(inputPath ,['fields_motion_2022-11-'...
     num2str(processingDay) '.mat']);

% Skip re-plotting when every expected figure is newer than the fields file,
% the GPS overlay, and this script; Force overrides.
stageInputs = {filename, fullfile(mvt.dayDir('gps', processingDay), ...
    ['CIRCLES_GPS_10Hz_2022-11-' num2str(processingDay) '.json'])};
[stale, staleReason] = mvt.isStale( ...
    mvt.expectedOutputs('macro', processingDay, opts), stageInputs, ...
    mvt.sources('plot_macroscopic_fields', opts), opts);
if ~stale
    mvt.log(opts, 'skip macroscopic field figures for 2022-11-%d: %s', ...
        processingDay, staleReason);
    return
end
mvt.log(opts, 'plot macroscopic fields for 2022-11-%d: %s', processingDay, staleReason);
if opts.Clean && ~opts.DryRun
    mvt.removeOutputs(mvt.expectedOutputs('macro', processingDay, opts), opts);
end
if opts.DryRun
    return
end

fprintf('Loading %s ...',filename), tic
load(filename)
fprintf(' Done (%0.0fsec).\n',toc)
if flagPlotAvTrajectories
    % p.resultsDir, not <dataRoot>/results: the results tree is relocatable
    % via MVT_RESULTS_DIR (make RESULTS=...), and hardcoding it here made this
    % stage fail with "Index exceeds array bounds" on any non-default tree.
    gpsDir = fullfile(p.resultsDir, 'gps');
    dataFiles = dir(fullfile(gpsDir, ['CIRCLES_GPS_10Hz_2022-11-'...
         num2str(processingDay)  '.json']));
    if isempty(dataFiles)
        error('mvt:plot_macroscopic_fields:missingGPS', ...
            'No assembled GPS file for 2022-11-%d in %s. Run the gps stage first.', ...
            processingDay, gpsDir);
    end
    filename = dataFiles(1).name; % if multiple files, use first one
    fprintf('Loading %s ...',filename), tic
    fid = fopen(fullfile(gpsDir, filename));
    data = fread(fid,inf);
    fclose(fid);
    fprintf(' Done (%0.0fsec).\n',toc)
    fprintf('Decoding data structures ...'), tic
    data = jsondecode(char(data'));
    fprintf(' Done (%0.0fsec).\n',toc)
    % Find indices of trajectories in direction and relevant lane(s)
    ind = find([data.direction]*direction>0&...
        ([data.assigned_lane]==lane|lane==0));
end
%========================================================================
% Plot results
%========================================================================
figRes = [2500 800]; figResZoom = [1250 800];
n_xticks = 16; % number of time ticks
close all
reportProgress = mvt.progress(length(plotFields), ...
    sprintf('macro 2022-11-%d', processingDay), 'Opts', opts);
for i = 1:length(plotFields)
    fieldPlot = plotFields{i};
    if direction<0
        directionName = 'Westbound';
        directionFname = 'west';
    else
        directionName = 'Eastound';
        directionFname = 'east';
    end
    if lane>0
        laneName = sprintf('lane %d',lane);
        laneFname = sprintf('lane%d',lane);
    else
        laneName = 'all lanes';
        laneFname = 'laneall';
    end
    fprintf('Produce figure ...'), tic
    t_local = datetime(t,'ConvertFrom','posixtime',...
        'TimeZone','America/Chicago');
    figure
    Val = field.value.(fieldPlot)'*field.factor.(fieldPlot);
    % Interpolate field for plotting
    tq = linspace(min(t),max(t),4*(length(t)-1)+1);
    xq = linspace(min(x),max(x),4*(length(x)-1)+1);
    [T,X] = meshgrid(t,x);
    [Tq,Xq] = meshgrid(tq,xq);
    Valq = interp2(T,X,Val,Tq,Xq);
    % Plot interpolated field
    imagesc(tq,xq/1000,Valq)
    % Set color limits on the heatmap
    cMax = colorbarLimit.(fieldPlot);
    if cMax > 0 
        caxis([0,cMax])     % preset max
    else
        caxis([0,prctile(Val(:),99)]) % max out slightly below maximum
    end
    % Time ticks location and text
    skip_t_data = round(length(t)/n_xticks);
    set(gca,'xtick',t(1:skip_t_data:end),'xticklabel',...
        string(datetime(t_local(1:skip_t_data:end),'Format','HH:mm:ss')))
    % Title of the plot
    [~,weekdayName] = weekday(t_local(1),'long');
    title(sprintf('%s (%s) on %s %s: %s',...
        directionName,laneName,weekdayName,...
        datetime(t_local(1),'Format','d-MMM-y'),...
        field.name.(fieldPlot)))
    hl = colorbar; hl.Label.String = field.unit.(fieldPlot);
    if flagPlotAvTrajectories % add AV trajectories 
        hold on
        for trjInd = ind % add AV trajectories to the plot (inactive/active)
            ind_e = find(data(trjInd).control_car'==0);
            traj_t = data(trjInd).timestamp(ind_e(1:skipTPlot:end));
            traj_x = data(trjInd).(subfieldNameX)(ind_e(1:skipTPlot:end));
            plot(traj_t,traj_x/1000,'w.','MarkerSize',3)
            ind_e = find(data(trjInd).control_car'==1);
            traj_t = data(trjInd).timestamp(ind_e(1:skipTPlot:end));
            traj_x = data(trjInd).(subfieldNameX)(ind_e(1:skipTPlot:end));
            plot(traj_t,traj_x/1000,'r.','MarkerSize',3)
        end
        hl = [plot(0,0,'r.','MarkerSize',26);...
            plot(0,0,'ko','MarkerSize',7)];
        hold off
        legend(hl,'AV engaged','AV disengaged','FontSize',24)
    end
    xlabel(sprintf('time in UTC%s',datetime(t_local(1),'Format','Z')))
    ylabel('position / km')
    fprintf(' Done (%0.0fsec).\n',toc)
    % Zoom in on time window
    lowHour = floor(timeZoomWindow(1)/100); 
    lowMinutes = mod(timeZoomWindow(1),100);
    lowShift = min(240,max(0,(lowHour - 6)*60 + lowMinutes)) ;
    upHour = floor(timeZoomWindow(2)/100); 
    upMinutes = mod(timeZoomWindow(2),100);
    upShift = max(-240,max(0,(upHour - 6)*60 + upMinutes)-240) ;
    axx = xlim;
    xlim(axx+diff(axx)*[lowShift upShift]/240) 
    if flagSaveFigures % save figure
        filename = sprintf('fig_field_%s_%s_%s_motion_%s',...
            datetime(t_local(1),'Format','yyyyMMdd'),...
            directionFname,laneFname,fieldPlot);
        if flagPlotAvTrajectories
            filename = [filename,'_av'];
        end
        filename = [filename,'_nature_large'];
        filename = fullfile(inputPath,filename);
        fprintf('Save figure in %s ...',filename), tic
        set(gcf,'Position',[10 50 figRes],'PaperPositionMode','auto')
        set(gca,'Position',[.023 .093 .92 .866])
        set(gca,'FontSize',20)
        print(filename,'-dpng','-r384');
        fprintf(' Done (%0.0fsec).\n',toc)
    end
    reportProgress(i, plotFields{i});
end
reportProgress();
end