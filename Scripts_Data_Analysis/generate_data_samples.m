% Script to plot graphs comparing statistics relative to distance from AVs
% by Sulaiman Almatrudi for CIRCLES energy team

clearvars -except DAY_TO_PROCESS
if ~exist('DAY_TO_PROCESS', 'var')
    global DAY_TO_PROCESS
    DAY_TO_PROCESS = input('Enter the day of Nov. 2022 MVT to generate AVs GPS data files (from 16 to 18): ');
end
fprintf('Starting the assembly of AV GPS data for %dth Nov. 2022\n', DAY_TO_PROCESS);

%========================================================================
% Locating data files
%========================================================================

fprintf('\nGenerating data relative to active AVs\n')

[parentDirectory, ~, ~] = fileparts(pwd);
folderPath = fullfile(parentDirectory,[ 'Data_2022-11-' char(num2str(DAY_TO_PROCESS)) '__MVT_Slim']);
if ~isfolder(folderPath)
    error('Folder %s does not exist.\n',folderPath)
end
I24files = dir([ folderPath '\I-24MOTION_slim_2022-11-' char(num2str(day_to_process)) '*.json']);
NRfiles = length(I24files);
if NRfiles < 24
    warning('Some I24 slim files for the day: %d, Nov. 2022 are missing.',DAY_TO_PROCESS)
end

%========================================================================
% Assemble data from JSON files
%========================================================================
% script outputs
samples_dist = []; %distance relative to engaged av (positive = behind an av, negative = ahead of an av)
samples_speed = []; %speed samples with downstream/upstream av
samples_fr = []; %fuel rate (g/s) samples with downstream/upstream av
samples_fcons = []; %fuel consumption (g/m) samples with downstream/upstream av
samples_class = []; %vehicle class for samples with downstream/upstream  av
samples_xpos = []; %x position for samples with downstream/upstream  av
samples_lane = []; %lane number for samples with downstream/upstream  av
samples_t = []; %timestamp for samples with downstream/upstream  av (seconds after 6am)

% data collection
for filInd =1:NRfiles
    fprintf('\n Loading original data file %d / %d ... ', filInd,NRfiles)
    tic
    filenameLoad = [folderPath '\' I24files(filInd).name];
    dataTemp = jsondecode(fileread(filenameLoad));
    fprintf('Done (%0.0fsec).\n',toc)
    tic
    fprintf('\n Processing data ... ')
    
    [dist, speed, class, fr, xpos, t, lane, fcons] = ...
        stats_to_av_dist(dataTemp,day_to_process);
    
    fprintf('Done (%0.0fsec).\n',toc)
    
    tic
    fprintf('\n appending data ... ')
    samples_dist = [samples_dist;dist];
    samples_speed = [samples_speed;speed];
    samples_class = [samples_class;class];
    samples_fr = [samples_fr;fr];
    samples_xpos = [samples_xpos;xpos];
    samples_t = [samples_t;t];
    samples_lane = [samples_lane;lane];
    samples_fcons = [samples_fcons;fcons];
    fprintf('Done (%0.0fsec).\n',toc)
end

%========================================================================
% Generate samples_for_distance_analysis
%========================================================================

%% save the output file to .mat data file
fprintf('\n Saving output file ... ');tic
save([parentDirectory '\Data_for_Data_Analysis\' 'samples_for_distance_analysis_' char(num2str(day_to_process)) '.mat'],'samples_*','-v7.3')
fprintf('Done (%0.0fsec).\n',toc)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% function definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [all_d, all_v, all_vc, all_fr, all_x, all_t, all_lane, all_fcons] = stats_to_av_dist(dataTemp,dayDate)

Max_Dist = 1000; %maximum distance from an engaged AV for samples to be collected

% initialize with appropriate data type to optimize for memory use
all_v = inf*ones(size(vertcat(dataTemp.timestamp)));
all_fr = inf*ones(size(vertcat(dataTemp.timestamp)));
all_d = inf*ones(size(vertcat(dataTemp.timestamp)));
all_x = int16(inf*ones(size(vertcat(dataTemp.timestamp))));
all_vc = uint8(all_x);
all_lane = uint8(all_x);
all_t = uint16(all_x);


jj1 = 1;

for i=1:length(dataTemp)
    veh = dataTemp(i);
    if veh.direction>0; continue;end % only aggregate for westbound traffic
    
    nemp_p = ~isempty(veh.distance_to_downstream_engaged_av_meters);
    nemp_n = ~isempty(veh.distance_to_upstream_engaged_av_meters);
    
    if ((nemp_p)||(nemp_n))
        
        for ii =1:length(veh.timestamp)
            
            %aggregate data for vehicles that have an AV downstream within Max_Dist
            if nemp_p &&(~isnan(veh.distance_to_downstream_engaged_av_meters(ii))) && ...
                    (veh.distance_to_downstream_engaged_av_meters(ii) <= Max_Dist)
                
                temp_d = veh.distance_to_downstream_engaged_av_meters(ii);
                
                all_v(jj1) = veh.speed_meters_per_second(ii);
                all_fr(jj1) = veh.fuel_rate_grams_per_second(ii);
                all_d(jj1) =   temp_d;
                
                all_vc(jj1) = uint8(veh.coarse_vehicle_class);
                all_x(jj1) = int16(veh.x_position_meters(ii));
                all_lane(jj1) = uint8(veh.lane_number);
                %record time in (s) after 6am of each test day
                all_t(jj1) = uint16(veh.timestamp(ii)-(1668772800-(18-dayDate)*24*60*60)); 
                jj1 =jj1+1;
            end
            %aggregate data for vehicles that have an AV upstream within Max_Dist
            if nemp_n && (~isnan(veh.distance_to_upstream_engaged_av_meters(ii))) && ...  
                    (veh.distance_to_upstream_engaged_av_meters(ii) >= -Max_Dist)
                
                temp_d = veh.distance_to_upstream_engaged_av_meters(ii);
                
                all_v(jj1) = veh.speed_meters_per_second(ii);
                all_fr(jj1) = veh.fuel_rate_grams_per_second(ii);
                all_d(jj1) =   temp_d;
                
                all_vc(jj1) = uint8(veh.coarse_vehicle_class);
                all_x(jj1) = int16(veh.x_position_meters(ii));
                all_lane(jj1) = uint8(veh.lane_number);
                % record time in (s) after 6am of each test day
                all_t(jj1) = uint16(veh.timestamp(ii)-(1668772800-(18-dayDate)*24*60*60)); 
                jj1 =jj1+1;
                
            end
        end
    end
end



all_d(jj1:end) = [];
all_v(jj1:end) = [];
all_vc(jj1:end) = [];
all_fr(jj1:end) = [];
all_x(jj1:end) = [];
all_t(jj1:end) = [];
all_lane(jj1:end) = [];
all_fcons = all_fr./ max(1e-6,all_v);

end