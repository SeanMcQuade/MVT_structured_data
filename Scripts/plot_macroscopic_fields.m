function [] = plot_macroscopic_fields(processingDay)
% Plot macroscopic fields based on I24-MOTION data for Nature.
% (C) 2025 Benjamin Seibold (edited by Sulaiman Almatrudi)
if nargin < 1
    error(['Specify the day of Nov. 2022 MVT to plot'...
        ' macroscopic fileds for (from 16 to 18)']);
end
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
[parentDirectory, ~, ~] = fileparts(pwd);
filename = fullfile(parentDirectory ,'Data','Data_for_Figures',['fields_motion_2022-11-'...
     num2str(processingDay) '.mat']);
fprintf('Loading %s ...',filename), tic
load(filename)
fprintf(' Done (%0.0fsec).\n',toc)
if flagPlotAvTrajectories
    dataFiles = dir(fullfile(parentDirectory, 'Data','Data_GPS',['CIRCLES_GPS_10Hz_2022-11-'...
         num2str(processingDay)  '.json']));
    filename = dataFiles(1).name; % if multiple files, use first one
    fprintf('Loading %s ...',filename), tic
    fid = fopen(fullfile(parentDirectory, 'Data','Data_GPS',filename));
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
        filename = fullfile(parentDirectory, 'Figures',filename);
        fprintf('Save figure in %s ...',filename), tic
        set(gcf,'Position',[10 50 figRes],'PaperPositionMode','auto')
        set(gca,'Position',[.023 .093 .92 .866])
        set(gca,'FontSize',20)
        print(filename,'-dpng','-r384');
        fprintf(' Done (%0.0fsec).\n',toc)
    end
end
end