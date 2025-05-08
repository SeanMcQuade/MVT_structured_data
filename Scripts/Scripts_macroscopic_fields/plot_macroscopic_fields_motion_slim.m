% Construct macroscopic fields based on I24-MOTION data files in the
% same directory.
% Specific for Nature
% (C) Benjamin Seibold (edited by Sulaiman Almatrudi)
clearvars -except DAY_TO_PROCESS
if ~exist('DAY_TO_PROCESS', 'var')
    global DAY_TO_PROCESS
    DAY_TO_PROCESS = input(['Enter the day of Nov. 2022 MVT to generate'...
        ' macroscopic fileds for (from 16 to 18): ']);
end
%========================================================================
% Parameters
%========================================================================
direction = -1; % -1=Westbound, 1=Eastbound
lane = 0; % 0=all lanes
%plot_fields = {'Rho','Q','F','F0','U','Phi','Phi0','Psi','Psi0'};
plot_fields = {'U'};
flag_save_figures = 1; % if true, save figure into png file
flag_plot_av_trajectories = 1; % if true, overlay fields with AV traces
skip_t_plot = 1; % sub-sample trajectories for plotting
subfield_name_x = 'x_position';

%========================================================================
% Load data files
%========================================================================
[parentDirectory, ~, ~] = fileparts(pwd);
filename = [ parentDirectory '\Data_Macroscopic_Fields\fields_motion_2022-11-'...
     num2str(DAY_TO_PROCESS) '.mat'];
fprintf('Loading %s ...',filename), tic
load(filename)
fprintf(' Done (%0.0fsec).\n',toc)
if flag_plot_av_trajectories
    data_files = dir([parentDirectory '\Data_GPS\CIRCLES_GPS_10Hz_2022-11-'...
         num2str(DAY_TO_PROCESS)  '.json']);
    filename = data_files(1).name; % if multiple files, use first one
    fprintf('Loading %s ...',filename), tic
    fid = fopen(filename);
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
fig_res = [2500 800]; fig_res_zoom = [1250 800];
n_xticks = 16; % number of time ticks
close all
for i = 1:length(plot_fields)
    field_plot = plot_fields{i};
    if direction<0
        direction_name = 'Westbound';
        direction_fname = 'west';
    else
        direction_name = 'Eastound';
        direction_fname = 'east';
    end
    if lane>0
        lane_name = sprintf('lane %d',lane);
        lane_fname = sprintf('lane%d',lane);
    else
        lane_name = 'all lanes';
        lane_fname = 'laneall';
    end
    fprintf('Produce figure ...'), tic
    t_local = datetime(t,'ConvertFrom','posixtime',...
        'TimeZone','America/Chicago');
    figure
    Val = field.value.(field_plot)'*field.factor.(field_plot);
    % Interpolate field for plotting
    %tq = (t(floor([1:.5:end]))+t(ceil([1:.5:end])))/2;
    %xq = (x(floor([1:.5:end]))+x(ceil([1:.5:end])))/2;
    tq = linspace(min(t),max(t),4*(length(t)-1)+1);
    xq = linspace(min(x),max(x),4*(length(x)-1)+1);
    [T,X] = meshgrid(t,x);
    [Tq,Xq] = meshgrid(tq,xq);
    Valq = interp2(T,X,Val,Tq,Xq);
    % Plot interpolated field
    imagesc(tq,xq/1000,Valq)
    %clim([0,prctile(Val(:),99)]) % max out slightly below maximum
    clim([0,0.20])
    skip_t_data = round(length(t)/n_xticks);
    set(gca,'xtick',t(1:skip_t_data:end),'xticklabel',...
        string(datetime(t_local(1:skip_t_data:end),'Format','HH:mm:ss')))
    [~,weekday_name] = weekday(t_local(1),'long');
    title(sprintf('%s (%s) on %s %s: %s',...
        direction_name,lane_name,weekday_name,...
        datetime(t_local(1),'Format','d-MMM-y'),...
        field.name.(field_plot)))
    hl = colorbar; hl.Label.String = field.unit.(field_plot);
    if flag_plot_av_trajectories
        hold on
        for j = ind % add trajectories
            ind_e = find(data(j).controller_engaged'==0);
            traj_t = data(j).timestamp(ind_e(1:skip_t_plot:end));
            traj_x = data(j).(subfield_name_x)(ind_e(1:skip_t_plot:end));
            plot(traj_t,traj_x/1000,'w.','MarkerSize',3)
            ind_e = find(data(j).controller_engaged'==1);
            traj_t = data(j).timestamp(ind_e(1:skip_t_plot:end));
            traj_x = data(j).(subfield_name_x)(ind_e(1:skip_t_plot:end));
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

    % Conduct zoom
    axx = xlim;
    xlim(axx+diff(axx)*[10 -10]/240) % 06:10-09:50
    if flag_save_figures % save figure
        filename = sprintf('fig_field_%s_%s_%s_motion_%s',...
            datetime(t_local(1),'Format','yyyyMMdd'),...
            direction_fname,lane_fname,field_plot);
        if flag_plot_av_trajectories
            filename = [filename,'_av'];
        end
        filename = [parentDirectory '\Figures\' filename,'_nature_large'];
        fprintf('Save figure in %s ...',filename), tic
        set(gcf,'Position',[10 50 fig_res],'PaperPositionMode','auto')
        %set(gca,'Position',[.023 .093 .901 .866])
        set(gca,'Position',[.023 .093 .929 .866])
        set(gca,'FontSize',20)
        print(filename,'-dpng','-r384');
        fprintf(' Done (%0.0fsec).\n',toc)
    end

    % Conduct zoom
    %xlim(axx+diff(axx)*[0.42 -0.17])
    %xlim(axx+diff(axx)*[100.8 -40.8]/240)
    xlim(axx+diff(axx)*[100 -40]/240) % 07:40-09:20
    if flag_save_figures % save figure
        filename = sprintf('fig_field_%s_%s_%s_motion_%s',...
            datetime(t_local(1),'Format','yyyyMMdd'),...
            direction_fname,lane_fname,field_plot);
        if flag_plot_av_trajectories
            filename = [filename,'_av'];
        end
        filename = [parentDirectory '\Figures\' filename,'_nature_zoom_large'];
        fprintf('Save figure in %s ...',filename), tic
        set(gcf,'Position',[10 50 fig_res_zoom],'PaperPositionMode','auto')
        %set(gca,'Position',[.045 .093 .852 .866])
        set(gca,'Position',[.045 .093 .858 .866])
        set(gca,'FontSize',20)
        print(filename,'-dpng','-r384');
        fprintf(' Done (%0.0fsec).\n',toc)
    end
end
