function [] = plot_microscopic_trajectories(processingDay)
% Produces time-space plots of I24-MOTION trajectories, visualizing each
% trajectory as a patch whose vertical extension precisely represents the
% vehicle's position and length. Also produces an example zoom window,
% both added to the full trajectory plot, and zoomed in.
% The code can produce low resolution (1600 figure rows) and high
% resolution versions of the figurem, which are of the following sizes,
% which can then be downscaled post-hoc as needed:
% Full: 12800, With zoom window added: 6400, Zoom itself: 6400.
% (C) 2025 by Benjamin Seibold and Sulaiman Almatrudi
if nargin < 1
    error(['Specify the day of Nov. 2022 MVT to generate'...
        ' macroscopic fileds for (from 16 to 18)']);
end
%========================================================================
% Parameters
%========================================================================
direction = -1; % -1=Westbound, 1=Eastbound
lane = 0; % 0=all lanes
ax_t = {'06:00:00','10:00:00'}; % time interval for trajectories
ax_x = [0,6500]; % space interval for trajectories 
skip_t_plot = 50; % sub-sample trajectories for plotting
% skip_t_plot = 5; % sub-sample trajectories for plotting
subfield_name_x = 'x_position_meters'; % name of x position field in the data
fig_res = [2500 800]; % figure base resolution
timeZoomWin = [0612 0619]; % time window to zoom into (in military time)
xZoomWin = [845 1495]; %845 1495];
n_xticks = 16; % number of time ticks
font_size = 20; % font size in figure
lowResDPI = 192; % [DPI] resolution of low  res. image to be saved
highResDPI = 1536; % [DPI] resolution of high res. image to be saved
zoomHighResDPI = 768; % [DPI] resolution of high res. zoomed image to be saved
flag_zoom_plot = 1; % generate a zoomed version of the plot 
flag_save_lowres = 1; % if true, save low resolution figures
flag_save_highres = 0; % if true, save high resolution figures 
flag_reduce_data_files = 1; % if true, generate .mat reduced data files to
                                % run the script faster on next iterations
flag_use_int32_vars = 1; % if true convert data to int32 to reduce memory usage
toInt32Factor = 1000; % factor to scale date before convering to int32 
ft2meterFactor = 0.3048; % [m/ft] conversion factor from feet to meter
%========================================================================
% Initialize
%========================================================================
% Find all data files in folder
[parentDirectory, ~, ~] = fileparts(pwd);
% directory above contains only the git repository
[dataRootDirectory, ~, ~] = fileparts(parentDirectory);

% Use slim version of the processed data if plotting westbound 
if direction < 0 % Westbound
    dataVersion = 'slim';
else            % Eastbound
    dataVersion = 'full';
end

inputPath = fullfile(dataRootDirectory, 'results', dataVersion, ...
    ['2022-11-', num2str(processingDay)]);

% dataFolder = fullfile(parentDirectory,'Data',...
%     ['Data_2022-11-' char(num2str(processingDay)) '__MVT_' dataVersion]);
if flag_reduce_data_files 
    data_files = dir(fullfile(inputPath, ...
        '*_reduced.mat'));
    
    % avoid processing files that start with .
    is_dotfile = startsWith({data_files.name},'.');
    data_files = data_files(~is_dotfile);

    if length(data_files) < 24
        reduce_data(inputPath);
        % now, the files should be there!
        data_files = dir(fullfile(inputPath, ...
            '*_reduced.mat'));

        % avoid processing files that start with .
        is_dotfile = startsWith({data_files.name},'.');
        data_files = data_files(~is_dotfile);

    end    
else
    data_files = dir(fullfile(inputPath, ...
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
    'InvertHardcopy','off','Color',[1 1 1])
set(gca,'Color',[0 0 0],'FontSize',font_size)
set(gca,'Position',[.023 .093 .957 .866])
t = linspace(ax_t(1),ax_t(2),n_xticks+1);
t_local = datetime(t,'ConvertFrom','posixtime','TimeZone','America/Chicago');
[~,weekday_name] = weekday(t_local(1),'long');
title(sprintf('%s (%s) on %s %s in UTC%s: trajectories',...
    direction_name,lane_name,weekday_name,...
    datetime(t_local(1),'Format','d-MMM-y'),...
    datetime(t_local(1),'Format','Z')))
xlabel('time'), ylabel('position / km')
axis ij
if flag_use_int32_vars
    set(gca,'xtick',(t-ax_t(1))*toInt32Factor,'xticklabel',...
    string(datetime(t_local,'Format','HH:mm:ss')))
    set(gca,'ytick',((0:7)*toInt32Factor),'yticklabel',...
        string((0:7)))
    xlim((ax_t-ax_t(1))*toInt32Factor), ylim(ax_x*toInt32Factor/1000)
    caxis([0,35*toInt32Factor])
else
    set(gca,'xtick',(t-ax_t(1)),'xticklabel',...
    string(datetime(t_local,'Format','HH:mm:ss')))
    set(gca,'ytick',(0:7),'yticklabel',...
        string((0:7)))
    xlim((ax_t-ax_t(1))), ylim(ax_x/1000)
    caxis([0,35])
end
cMap = interp1([0;1],[1 0 0;0 1 0],linspace(0,1,256));
colormap(cMap)
hl = colorbar('Fontsize',font_size); hl.Label.String = 'vehicle speed / m/s';
axx = xlim; xlim(axx+diff(axx)*[10 -10]/240) % Show plot 06:10-09:50
hold on
fprintf(' Done (%0.0fsec).\n',toc)
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
    if isfield(data,'direction')
        ind = find([data.direction]*direction>0&([data.lane_number]==lane|lane==0));
    else
        ind = find(([data.lane_number]==lane|lane==0));
    end
    veh_lengths = vertcat(data.length); % vector of vehicle lengths
    veh_lengths = veh_lengths(ind); % only those in index set
    [~,si] = sort(veh_lengths,'descend'); % sort by decreasing vehicle length
    ind = ind(si); % reorder index set accordingly
    data = data(ind);
    % Process trajectories
    fprintf('Adding %d trajectories to plot ...',length(ind)), tic
    for j = 1:length(data) % loop over used trajectories
        traj_len = data(j).length*ft2meterFactor; % length of vehicle [m]
        traj_t = data(j).timestamp; 
        traj_x = data(j).(subfield_name_x); % vehicle position [m]
        % Calculate velocity
        traj_v = (traj_x([2:end,end])-traj_x([1,1:end-1]))./...
            (traj_t([2:end,end])-traj_t([1,1:end-1]))*direction;
        % Subsample trajectories for plotting
        traj_n = length(traj_t); % number of trajectory points
        traj_ns = max(ceil(traj_n/skip_t_plot),2); % number of segments
        traj_ind = round(linspace(1,traj_n,traj_ns)); % subsample vector
        traj_t = traj_t(traj_ind);
        traj_x = traj_x(traj_ind);
        traj_v = traj_v(traj_ind);
        patch_t = [traj_t;traj_t(end:-1:1)] - ax_t(1) ; 
        patch_x = [traj_x;traj_x(end:-1:1)+traj_len*direction]/1000;
        patch_v = [traj_v;traj_v(end:-1:1)];
        if flag_use_int32_vars
            patch_t = int32(patch_t*toInt32Factor);
            patch_x = int32(patch_x*toInt32Factor);
            patch_v = int32(patch_v*toInt32Factor);
        end
         patch(patch_t, patch_x, patch_v,'EdgeColor','None')
    end
    fprintf(' Done (%0.0fsec).\n',toc)
end
clear data  veh_lengths
% Finish drawing everything 
fprintf('Drawing full plot ...'), tic
drawnow
fprintf(' Done (%0.0fsec).\n',toc)
%========================================================================
% Save figure
%========================================================================
fileName = sprintf('fig_motion_trajectories_%s_%s_%s',...
    datetime(t_local(1),'Format','yyyyMMdd'),...
    direction_fname,lane_fname);
folderName = fullfile(dataRootDirectory, 'results', 'figures', ...
    ['2022-11-', num2str(processingDay)]);
if ~isfolder(folderName)
    print('Creating folder %s', folderName)
    mkdir(folderName);
end
fileName = fullfile(folderName,fileName);
% Save figure
if flag_save_lowres
    print_figure([fileName,'_lowres.png'],lowResDPI)
end
if flag_save_highres
    print_figure([fileName,'_highres.png'],highResDPI)
end

%========================================================================
% Add zoom window
%========================================================================
if flag_zoom_plot
    % Get zoom window position on the plot
    lowHour = floor(timeZoomWin(1)/100); 
    lowMinutes = mod(timeZoomWin(1),100);
    lowShift = min(240,max(0,(lowHour - 6)*60 + lowMinutes)) ;
    upHour = floor(timeZoomWin(2)/100); 
    upMinutes = mod(timeZoomWin(2),100);
    upShift = min(240,max(0,(upHour - 6)*60 + upMinutes)) ;
    zoom_t = [lowShift upShift]*60 ;
    if flag_use_int32_vars 
        zoom_t = zoom_t * toInt32Factor;
    end
    axy = ylim*1000; 
    zoom_x = diff(axy)*xZoomWin/diff(ax_x);
    % Plot window
    hb = plot(zoom_t([1 2 2 1 1 2]),zoom_x([1 1 2 2 1 1])/1000,...
        '-','LineWidth',2.5,'Color',[1 1 0]);
    % Save figure with window plotted 
    print_figure([fileName,'_zoomwin_lowres.png'],lowResDPI)

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
    fac = (diff(zoom_t)/diff(xlim))/(diff(zoom_x/1000)/diff(ylim))*posgca(3)/posgca(4);
    xlim(zoom_t), ylim(zoom_x/1000)
    set(hb,'LineWidth',10)
    set(gcf,'Color',[0,0,0],'Position',[10 50 fig_res(1)*fac fig_res(2)])
    set(gca,'Visible','off','Position',[0 0 1 1])
    fprintf(' Done (%0.0fsec).\n',toc)
    % Save figure
    if flag_save_lowres
        print_figure([fileName,'_zoom_lowres.png'],lowResDPI)
    end
    if flag_save_highres
        print_figure([fileName,'_zoom_highres.png'],zoomHighResDPI)
    end
end
end

%========================================================================
%========================================================================
% Functions
%========================================================================
%========================================================================
function print_figure(filename,res)
fprintf('Save figure in %s ...',filename), tic
fig = gcf; 
exportgraphics(fig, filename, 'Resolution', res);
fprintf(' Done (%0.0fsec).\n',toc)
end

function [] = reduce_data(inputPath)

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

data_files = dir(fullfile(inputPath,'*.json'));

% avoid processing files that start with .
is_dotfile = startsWith({data_files.name},'.');
data_files = data_files(~is_dotfile);

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