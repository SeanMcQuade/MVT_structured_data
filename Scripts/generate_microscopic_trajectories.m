% function [] = generate_microscopic_trajectories(processingDay)
% Produces time-space plots of I24-MOTION trajectories, visualizing each
% trajectory as a patch whose vertical extension precisely represents the
% vehicle's position and length. Also produces an example zoom window,
% both added to the full trajectory plot, and zoomed in.
% The code can produce low resolution (1600 figure rows) and high
% resolution versions of the figurem, which are of the following sizes,
% which can then be downscaled post-hoc as needed:
% Full: 12800, With zoom window added: 6400, Zoom itself: 6400.
% (C) 2025/06/16 by Benjamin Seibold
% if nargin < 1
%     error(['Specify the day of Nov. 2022 MVT to generate'...
%         ' macroscopic fileds for (from 16 to 18)']);
% end
%========================================================================
% Parameters
%========================================================================
direction = -1; % -1=Westbound, 1=Eastbound
lane = 0; % 0=all lanes
ax_t = {'06:00:00','10:00:00'}; % time interval for fields
ax_x = [0,6500]; % space interval for fields
skip_t_plot = 5; % sub-sample trajectories for plotting
subfield_name_x = 'x_position_meters';
fig_res = [2500 800]; % figure base resolution
n_xticks = 16; % number of time ticks
font_size = 20; % font size in figure
zoom_t_rel = [.44 .49]; zoom_x_rel = [.13 .23]; % zoom window coords %NOTE use clearer paramter
flag_save_lowres = 1; % if true, save low resolution figures
flag_save_highres = 0; % if true, save high resolution figures 
flag_reduce_data_files = 1; % if true, generate .mat reduced data files to
% run the script faster on next iterations 
ft2meterFactor = 0.3048; % [m/ft] conversion factor from feet to meter
%========================================================================
% Initialize
%========================================================================
% Find all data files in folder
[parentDirectory, ~, ~] = fileparts(pwd);
% Use slim version of the processed data if plotting westbound 
if direction < 0 % Westbound
    dataVersion = 'Slim';
else            % Eastbound
    dataVersion = 'Full';
end
dataFolder = fullfile(parentDirectory,'Data',...
    ['Data_2022-11-' char(num2str(processingDay)) '__MVT_' dataVersion]);
if flag_reduce_data_files 
    data_files = dir(fullfile(dataFolder, ...
        '*_reduced.mat'));
    if length(data_files) < 24
        reduce_data(dataFolder);
    end    
else
    data_files = dir(fullfile(dataFolder, ...
        ['I-24*' char(num2str(processingDay)) '*.json']));
end
date_data = ['2022-11-' num2str(processingDay)]; % date of the data files
ax_t = [posixtime(datetime([date_data,' ',ax_t{1}],'TimeZone','America/Chicago')),...
    posixtime(datetime([date_data,' ',ax_t{2}],'TimeZone','America/Chicago'))];
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
if flag_save_highres
    warning('Figures will be saved in high resolution. This may take 30 to 60 minutes.')
end

%========================================================================
% Initialize figure
%========================================================================
close all
fprintf('Initialize figure ...'), tic
figure('Position',[10 50 fig_res],'PaperPositionMode','auto',...
    'InvertHardcopy','off','Color',[1 1 1])%,'Visible','off')
set(gca,'Color',[0 0 0],'FontSize',font_size)
set(gca,'Position',[.023 .093 .957 .866]) % fontsize=20 ; [2500 800]
t = linspace(ax_t(1),ax_t(2),n_xticks+1);
t_local = datetime(t,'ConvertFrom','posixtime','TimeZone','America/Chicago');
set(gca,'xtick',t,'xticklabel',...
    string(datetime(t_local,'Format','HH:mm:ss')))
[~,weekday_name] = weekday(t_local(1),'long');
title(sprintf('%s (%s) on %s %s in UTC%s: trajectories',...
    direction_name,lane_name,weekday_name,...
    datetime(t_local(1),'Format','d-MMM-y'),...
    datetime(t_local(1),'Format','Z')))
xlabel('time'), ylabel('position / km')
axis ij
xlim(ax_t), ylim(ax_x/1000)
caxis([0,35])
cMap = interp1([0;1],[1 0 0;0 1 0],linspace(0,1,256));
colormap(cMap)
hl = colorbar('Fontsize',font_size); hl.Label.String = 'vehicle speed / m/s';
axx = xlim; xlim(axx+diff(axx)*[10 -10]/240) % 06:10-09:50
hold on
fprintf(' Done (%0.0fsec).\n',toc)
%%
%========================================================================
% Go through files/trajectories and plot them
%========================================================================
for fileInd = 1:length(data_files) % loop over relevant files
    % Load data file
    fileName = data_files(fileInd).name;
    fprintf('Loading %s ...',fileName), tic
    filePath = fullfile(data_files(fileInd).folder, fileName);
    if flag_reduce_data_files
        load(filePath) % if using reduced data files, load .mat file
    else
        data = jsondecode(fileread(filePath)); % if not using reduced files, 
                                                %load json files
    end
    fprintf(' Done (%0.0fsec).\n',toc)
    % Find indices of trajectories in direction and relevant lane(s)
    if isfield(data,'direction') %&& ~isempty([data.direction])
        ind = find([data.direction]*direction>0&([data.lane_number]==lane|lane==0));
    else
        ind = find(([data.lane_number]==lane|lane==0));
    end
    veh_lengths = vertcat(data.length); % vector of vehicle lengths
    veh_lengths = veh_lengths(ind); % only those in index set
    [~,si] = sort(veh_lengths,'descend'); % sort by decreasing vehicle length
    ind = ind(si); % reorder index set accordingly
    % Process trajectories
    fprintf('Plotting %d trajectories ...',length(ind)), tic
    for j = ind % loop over used trajectories
        traj_len = data(j).length*ft2meterFactor; % length of vehicle [m]
        traj_t = data(j).timestamp; 
        traj_x = data(j).(subfield_name_x); % vehicle position [m]
        % Calculate velocity
        traj_v = (traj_x([2:end,end])-traj_x([1,1:end-1]))./...
            (traj_t([2:end,end])-traj_t([1,1:end-1]))*direction;
        % Subsample trajectories
        traj_n = length(traj_t); % number of trajectory points
        traj_ns = max(ceil(traj_n/skip_t_plot),2); % number of segments
        traj_ind = round(linspace(1,traj_n,traj_ns)); % subsample vector
        traj_t = traj_t(traj_ind);
        traj_x = traj_x(traj_ind);
        traj_v = traj_v(traj_ind);
        % Plot trajectory as patch
        patch([traj_t;traj_t(end:-1:1)],...
            [traj_x;traj_x(end:-1:1)+traj_len*direction]/1000,... %NOTE  
            [traj_v;traj_v(end:-1:1)],'EdgeColor','None')
    end
    fprintf(' Done (%0.0fsec).\n',toc)
end
clear data  veh_lengths
%========================================================================
% Save figure
%========================================================================
fileName = sprintf('fig_motion_trajectories_%s_%s_%s',...
    datetime(t_local(1),'Format','yyyyMMdd'),...
    direction_fname,lane_fname);
% Save figure
if flag_save_lowres
    print_figure([fileName,'_lowres'],'-r192')
end
if flag_save_highres
    print_figure([fileName,'_highres'],'-r1536')
end

%========================================================================
% Add zoom window
%========================================================================
zoom_t = ax_t(1)+diff(ax_t)*zoom_t_rel;
zoom_x = ax_x(1)+diff(ax_x)*zoom_x_rel;
hb = plot(zoom_t([1 2 2 1 1 2]),zoom_x([1 1 2 2 1 1])/1000,...
    '-','LineWidth',2.5,'Color',[1 1 0]);
% Save figure
if flag_save_lowres
    print_figure([fileName,'_zoomwin_lowres'],'-r192')
end
if flag_save_highres
    print_figure([fileName,'_zoomwin_highres'],'-r768')
end

%========================================================================
% Show zoom itself
%========================================================================
% Remove all trajectories fully outside of zoom window
fprintf('Remove trajectories outside zoom window, and zoom in ...'), tic
h_traj = findobj(gca,'type','patch');
ind_delete = false(1,length(h_traj));
for j = 1:length(h_traj)
    ind_delete(j) = all(h_traj(j).Vertices(:,1)<zoom_t(1))||...
        all(h_traj(j).Vertices(:,1)>zoom_t(2))||...
        all(h_traj(j).Vertices(:,2)<zoom_x(1)/1000)||...
        all(h_traj(j).Vertices(:,2)>zoom_x(2)/1000);
end
delete(h_traj(ind_delete))
% Zoom into zoom window
posgca = get(gca,'Position');
fac = diff(zoom_t_rel)/diff(zoom_x_rel)*posgca(3)/posgca(4);
xlim(zoom_t), ylim(zoom_x/1000)
set(hb,'LineWidth',10)
set(gcf,'Color',[0,0,0],'Position',[10 50 fig_res(1)*fac fig_res(2)])
set(gca,'Visible','off','Position',[0 0 1 1])
fprintf(' Done (%0.0fsec).\n',toc)
% Save figure
if flag_save_lowres
    print_figure([fileName,'_zoom_lowres'],'-r192')
end
if flag_save_highres
    print_figure([fileName,'_zoom_highres'],'-r768')
end
% end

%========================================================================
%========================================================================
% Functions
%========================================================================
%========================================================================
function print_figure(filename,res)
fprintf('Save figure in %s ...',filename), tic
print(filename,'-dpng',res);
fprintf(' Done (%0.0fsec).\n',toc)
end

function [] = reduce_data(dataFolder)

removed_fields = {
    'trajectory_id',...
    'y_position_corrected_meters',...
    'width',...
    'height',...
    'total_distance_traversed_meters',...
    'acceleration_meters_per_second_per_second',...
    'road_grade_radians',...
    'reference_a1_meters_per_second_per_second',...
    'reference_a2_meters_per_second_per_second',...
    'total_fuel_consumed_grams',...
    'total_fuel_consumed_gallons',...
    'total_fuel_economy_mpg',...
    'total_fuel_consumed_flat_road_grams',...
    'total_fuel_consumed_flat_road_gallons',...
    'total_fuel_economy_flat_road_mpg',...
    'reference_fuel_rate_grams_per_second',...
    'percent_reference_infeasibility',...
    'total_reference_fuel_consumed_grams',...
    'reference_fuel_rate_flat_road_grams_per_second',...
    'percent_reference_infeasibility_flat_road',...
    'total_reference_fuel_consumed_flat_road_grams',...
    'downstream_av_id',...
    'distance_to_downstream_av_meters',...
    'downstream_engaged_av_id',...
    'distance_to_downstream_engaged_av_meters'};

data_files = dir(fullfile(dataFolder,'*.json'));
for fileInd = 1:length(data_files)
    filename = data_files(fileInd).name;
    filenameSave = [filename(1:end-5),'_reduced.mat'];
    fprintf('Reducing data and saving as %s ...',filename), tic
    % Load data file
    data = jsondecode(fileread(fullfile(data_files(fileInd).folder,filename)));
    % Reduce data
    fieldsToDelete = intersect(fieldnames(data(1)), removed_fields);
    data = rmfield(data,fieldsToDelete);
    % Save new data file
    save(fullfile(data_files(fileInd).folder,filenameSave),'data')
    fprintf(' Done (%0.0fsec).\n',toc)
end
end