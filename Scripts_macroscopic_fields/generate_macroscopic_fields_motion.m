% Construct macroscopic fields based on I24-MOTION data files in the
% same directory.
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
hx = 100; % [m] window size in x
direction = -1; % -1=Westbound, 1=Eastbound
lane = 0; % 0=all lanes
ax_t = {'06:00:00','10:00:00'}; % time interval for fields
ax_x = [0,6500]; % space interval for fields
t_res = 5; % [s] resolution in t
x_res = 50; % [m] resolution in x
v_char = 20; % [m/s] characteristic speed (for window size in t)
skip_t_compute = 5; % sub-sample trajectories for integrals
flag_plot_field = 0; % if true, plot one field in the end
flag_plot_trajectories = 0; % if true, overlay field with trajectory
skip_t_plot = 20; % sub-sample trajectories for plotting
kernel = 'box';
created_fields = {'Rho','Q','F','F0','U','Phi','Phi0','Psi','Psi0'};
subfield_name_x = 'x_position_meters';

%========================================================================
% Initialize
%========================================================================
ht = hx/v_char;
if strcmp(kernel,'box')
    fac = 1/(4*ht*hx); % normalizing factor
elseif strcmp(kernel,'Gaussian')
    wt = ht/sqrt(3); wx = hx/sqrt(3); % kernel widths to have same std as box
    ht = wt*3; hx = wx*3; % make box size 3 times the std
    fac = 1/(2*pi*wt*wx); % normalizing factor
end
% Find all data files in folder
[parentDirectory, ~, ~] = fileparts(pwd);
dataFolderPath = fullfile(parentDirectory,['Data_2022-11-' num2str(DAY_TO_PROCESS) '__MVT_Slim']);
data_files = dir([dataFolderPath '\I-24MOTION_slim_*.json']);
date_data = ['2022-11-' num2str(DAY_TO_PROCESS)]; % date of the data files
ax_t = [posixtime(datetime([date_data,' ',ax_t{1}],'TimeZone','America/Chicago')),...
    posixtime(datetime([date_data,' ',ax_t{2}],'TimeZone','America/Chicago'))];

%========================================================================
% Create grid in (t,x) plane and initialize fields
%========================================================================
clear field
t = ax_t(1):t_res:ax_t(2); x = ax_x(1):x_res:ax_x(2);
[T,X] = ndgrid(t,x);
if any(strcmp(created_fields,'Rho'))||any(strcmp(created_fields,'U'))||...
        any(strcmp(created_fields,'Phi'))||any(strcmp(created_fields,'Phi0'))
    field.name.Rho = 'vehicle density'; field.unit.Rho = 'veh/km';
    field.factor.Rho = 1000; field.value.Rho = T*0;
end
if any(strcmp(created_fields,'Q'))||any(strcmp(created_fields,'U'))||...
        any(strcmp(created_fields,'Psi'))||any(strcmp(created_fields,'Psi0'))
    field.name.Q = 'flow rate'; field.unit.Q = 'veh/h';
    field.factor.Q = 3600; field.value.Q = T*0;
end
if any(strcmp(created_fields,'F'))||any(strcmp(created_fields,'Phi'))||...
        any(strcmp(created_fields,'Psi'))
    field.name.F = 'fuel rate density'; field.unit.F = 'g/s';
    field.factor.F = 1; field.value.F = T*0;
end
if any(strcmp(created_fields,'F0'))||any(strcmp(created_fields,'Phi0'))||...
        any(strcmp(created_fields,'Psi0'))
    field.name.F0 = 'fuel rate density (flat road)'; field.unit.F0 = 'g/s';
    field.factor.F0 = 1; field.value.F0 = T*0;
end

%========================================================================
% Go through files/trajectories and successively increment fields
%========================================================================
for file_j = 1:length(data_files) % loop over relevant files
    % Load data file
    filename = data_files(file_j).name;
    fprintf('Loading and decoding %s ...',filename), tic
    data = jsondecode(fileread(filename));
    fprintf(' Done (%0.0fsec).\n',toc)
    % Find indices of trajectories in direction and relevant lane(s)
    if isfield(data,'direction')
        ind = find([data.direction]*direction>0&([data.lane_number]==lane|lane==0));
    else
        ind = find(([data.lane_number]==lane|lane==0));
    end
    % Process trajectories
    fprintf('Processing %d trajectories ...',length(ind)), tic
    for j = ind % loop over used trajectories
        traj_t = data(j).timestamp;
        traj_x = data(j).(subfield_name_x);
        if isfield(field.value,'Q')
            traj_v = (traj_x([2:end,end])-traj_x([1,1:end-1]))./...
                (traj_t([2:end,end])-traj_t([1,1:end-1]))*direction;
            traj_v = traj_v(1:skip_t_compute:end);
        end
        traj_t = traj_t(1:skip_t_compute:end);
        traj_x = traj_x(1:skip_t_compute:end);
        if isfield(field.value,'F')
            traj_f = data(j).fuel_rate_grams_per_second(1:skip_t_compute:end);
        end
        if isfield(field.value,'F0')
            traj_f0 = data(j).fuel_rate_flat_road_grams_per_second(1:skip_t_compute:end);
        end
        traj_n = length(traj_t);
        if traj_n>1
            traj_dt = diff(traj_t);
            traj_dt = (traj_dt([1,1:end])+traj_dt([1:end,end]))/2;
        else
            traj_dt = (data(j).timestamp(2)-data(j).timestamp(1))*skip_t_compute;
            warning('Trajectory of length 1 encountered.')
        end
        for i = 1:traj_n % loop over entries in given trajectory
            ind_box_t = abs(t-traj_t(i))<=ht; % window size in t
            ind_box_x = abs(x-traj_x(i))<=hx; % window size in x
            if strcmp(kernel,'Gaussian')
                G = exp(-(((T(ind_box_t,ind_box_x)-traj_t(i))/wt).^2+...
                    ((X(ind_box_t,ind_box_x)-traj_x(i))/wx).^2)/2);
            elseif strcmp(kernel,'box')
                G = T(ind_box_t,ind_box_x)*0+1;
            end
            if isfield(field.value,'Rho')
                field.value.Rho(ind_box_t,ind_box_x) = ...
                    field.value.Rho(ind_box_t,ind_box_x)+traj_dt(i)*G;
            end
            if isfield(field.value,'Q')
                field.value.Q(ind_box_t,ind_box_x) = ...
                    field.value.Q(ind_box_t,ind_box_x)+traj_dt(i)*traj_v(i)*G;
            end
            if isfield(field.value,'F')
                field.value.F(ind_box_t,ind_box_x) = ...
                    field.value.F(ind_box_t,ind_box_x)+traj_dt(i)*traj_f(i)*G;
            end
            if isfield(field.value,'F0')
                field.value.F0(ind_box_t,ind_box_x) = ...
                    field.value.F0(ind_box_t,ind_box_x)+traj_dt(i)*traj_f0(i)*G;
            end
        end
    end
    fprintf(' Done (%0.0fsec).\n',toc)
end
% Normalize fields by kernal factors
if isfield(field.value,'Rho'), field.value.Rho = field.value.Rho*fac; end
if isfield(field.value,'Q'), field.value.Q = field.value.Q*fac; end
if isfield(field.value,'F'), field.value.F = field.value.F*fac; end
if isfield(field.value,'F0'), field.value.F0 = field.value.F0*fac; end

%========================================================================
% Derived fields
%========================================================================
if any(strcmp(created_fields,'U'))
    field.name.U = 'bulk velocity'; field.unit.U = 'm/s';
    field.factor.U = 1; field.value.U = field.value.Q./field.value.Rho;
    field.value.U(field.value.Rho<1e-3) = nan;
end
if any(strcmp(created_fields,'Phi'))
    field.name.Phi = 'bulk fuel rate'; field.unit.Phi = 'g/s';
    field.factor.Phi = 1; field.value.Phi = field.value.F./field.value.Rho;
    field.value.Phi(field.value.Rho<1e-3) = nan;
end
if any(strcmp(created_fields,'Phi0'))
    field.name.Phi0 = 'bulk fuel rate (flat road)'; field.unit.Phi0 = 'g/s';
    field.factor.Phi0 = 1; field.value.Phi0 = field.value.F0./field.value.Rho;
    field.value.Phi0(field.value.Rho<1e-3) = nan;
end
if any(strcmp(created_fields,'Psi'))
    field.name.Psi = 'bulk fuel consumption'; field.unit.Psi = 'g/m';
    field.factor.Psi = 1; field.value.Psi = field.value.F./field.value.Q;
    field.value.Psi(field.value.Q<1e-2) = nan;
end
if any(strcmp(created_fields,'Psi0'))
    field.name.Psi0 = 'bulk fuel consumption (flat road)'; field.unit.Psi0 = 'g/m';
    field.factor.Psi0 = 1; field.value.Psi0 = field.value.F0./field.value.Q;
    field.value.Psi0(field.value.Q<1e-2) = nan;
end

%========================================================================
% Save fields
%========================================================================
filename = [parentDirectory '\Data_Macroscopic_Fields\fields_motion_2022-11-'...
    num2str(DAY_TO_PROCESS)  '.mat'];
fprintf('Saving file %s ...',filename), tic
save(filename,'field','t','x','direction','lane')
fprintf(' Done (%0.0fsec).\n',toc)

%========================================================================
% Plot results
%========================================================================
field_plot = 'Rho';
n_xticks = 20; % number of time ticks

if flag_plot_field
    if direction<0
        direction_name = 'Westbound';
    else
        direction_name = 'Eastound';
    end
    if lane>0
        lane_name = sprintf('lane %d',lane);
    else
        lane_name = 'all lanes';
    end
    fprintf('Produce figure ...'), tic
    t_local = datetime(t,'ConvertFrom','posixtime','TimeZone','America/Chicago');
    clf
    imagesc(t,x/1000,field.value.(field_plot)'*field.factor.(field_plot))
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
    if flag_plot_trajectories
        hold on
        for j = ind % add trajectories
            traj_t = data(j).timestamp(1:skip_t_plot:end);
            traj_x = data(j).(subfield_name_x)(1:skip_t_plot:end);
            plot(traj_t,traj_x/1000,'w-')
        end
        hold off
    end
    xlabel('time'), ylabel('position / km')
    fprintf(' Done (%0.0fsec).\n',toc)
end
