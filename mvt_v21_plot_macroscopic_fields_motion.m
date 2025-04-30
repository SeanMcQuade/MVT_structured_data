% Construct macroscopic fields based on I24-MOTION data files in the
% same directory.
% (C) 2023/09/24 by Benjamin Seibold

%========================================================================
% Parameters
%========================================================================
direction = -1; % -1=Westbound, 1=Eastbound
lane = 0; % 0=all lanes
plot_fields = {'Rho','Q','F','F0','U','Phi','Phi0','Psi','Psi0'};
flag_save_figures = 1; % if true, save figure into png file
flag_plot_av_trajectories = 1; % if true, overlay fields with AV traces
skip_t_plot = 1; % sub-sample trajectories for plotting
subfield_name_x = 'x_position';

%========================================================================
% Load data files
%========================================================================
filename = 'fields_motion.mat';
fprintf('Loading %s ...',filename), tic
load(filename)
fprintf(' Done (%0.0fsec).\n',toc)
if flag_plot_av_trajectories
    data_files = dir('CIRCLES_GPS_10Hz_????-??-??.json');
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
fig_res = [2500 800];
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
    imagesc(t,x/1000,Val)
    caxis([0,prctile(Val(:),99)]) % max out slightly below maximum
    skip_t_data = round(length(t)/n_xticks);
    set(gca,'xtick',t(1:skip_t_data:end),'xticklabel',...
        string(datetime(t_local(1:skip_t_data:end),'Format','HH:mm:ss')))
    [~,weekday_name] = weekday(t_local(1),'long');
    title(sprintf('%s (%s) on %s %s in UTC%s: %s',...
        direction_name,lane_name,weekday_name,...
        datetime(t_local(1),'Format','d-MMM-y'),...
        datetime(t_local(1),'Format','Z'),...
        field.name.(field_plot)))
    hl = colorbar; hl.Label.String = field.unit.(field_plot);
    if flag_plot_av_trajectories
        hold on
        for j = ind % add trajectories
            ind_e = find(data(j).controller_engaged'==0);
            traj_t = data(j).timestamp(ind_e(1:skip_t_plot:end));
            traj_x = data(j).(subfield_name_x)(ind_e(1:skip_t_plot:end));
            plot(traj_t,traj_x/1000,'w.','MarkerSize',2)
            ind_e = find(data(j).controller_engaged'==1);
            traj_t = data(j).timestamp(ind_e(1:skip_t_plot:end));
            traj_x = data(j).(subfield_name_x)(ind_e(1:skip_t_plot:end));
            plot(traj_t,traj_x/1000,'r.','MarkerSize',2)
        end
        hl = [plot(0,0,'r.','MarkerSize',26);...
            plot(0,0,'ko','MarkerSize',7)];
        hold off
        legend(hl,'AV engaged','AV disengaged')
    end
    xlabel('time'), ylabel('position / km')
    fprintf(' Done (%0.0fsec).\n',toc)
    if flag_save_figures % save figure
        filename = sprintf('fig_field_%s_%s_%s_motion_%s',...
            datetime(t_local(1),'Format','yyyyMMdd'),...
            direction_fname,lane_fname,field_plot);
        if flag_plot_av_trajectories
            filename = [filename,'_av'];
        end
        fprintf('Save figure in %s ...',filename), tic
        set(gcf,'Position',[10 50 fig_res],'PaperPositionMode','auto')
        set(gca,'Position',[.023 .093 .901 .866]) % fontsize=20
        set(gca,'FontSize',20)
        %print(filename,'-dpng','-r192');
        print(filename,'-dpng','-r384');
        fprintf(' Done (%0.0fsec).\n',toc)
    end
end
