% Run to generate data samples in .mat format for generating fuel consumption figures.
% (C) 11/8/2024 by Sulaiman Almatrudi for CIRCLES energy team; further adapted by Sean McQuade for CIRCLES scenario team.
% This licensed under BSD-3 clause license: https://opensource.org/license/bsd-3-clause
tic
DAYS_TO_PROCESS = [16,17,18]; % /NOV/2022 Enter 16, 17, or 18 here.

%========================================================================
% Generate samples_for_distance_analysis
%========================================================================

for day = DAYS_TO_PROCESS

    %% data preprocessing to generate data samples in .mat
    av_file_name = ['samples_for_distance_analysis_' char(num2str(day)) '.mat'];
    fprintf('\nGenerating data relative to active AVs\n')
    generate_data_samples_v10(day);

end

%% Functions used in this script 
%%%%%%%%%%%%%%%%%%%%%% auxiliary function 1

function [] = generate_data_samples_v10(day_to_process)

folder = [ '../Data/Data_2022-11-' char(num2str(day_to_process)) '__MVT_Data_Slim'];
I24_files_in_dir = dir([ folder '\I-24*' char(num2str(day_to_process)) '*.json']);
%I24_files_in_dir = dir([ 'I-24*' char(num2str(day_to_process)) '*.json']);
NR_files = length(I24_files_in_dir);

%% script outputs
samples_dist = []; %distance relative to engaged av (positive = behind an av, negative = ahead of an av)
samples_speed = []; %speed samples with downstream/upstream av
samples_fr = []; %fuel rate (g/s) samples with downstream/upstream av
samples_fcons = []; %fuel consumption (g/m) samples with downstream/upstream av
samples_class = []; %vehicle class for samples with downstream/upstream  av
samples_xpos = []; %x position for samples with downstream/upstream  av
samples_lane = []; %lane number for samples with downstream/upstream  av
samples_t = []; %timestamp for samples with downstream/upstream  av (seconds after 6am)

%% data collection
for file_nr =1:NR_files
    fprintf(sprintf('\n Loading original data file %d / %d \n', file_nr,NR_files ))
    tic
    filename_load = [folder '\' I24_files_in_dir(file_nr).name];
    data_v22 = jsondecode(fileread(filename_load));
    toc
    if ~isfield(data_v22,'distance_to_upstream_engaged_av_meters');fprintf('no AVs on the road');continue;end
    tic
    fprintf('\n Processing data\n')

    [dist, speed, class, fr, xpos, t, lane, fcons] = stats_to_av_dist(data_v22,day_to_process);

    toc
    
    tic
    fprintf('\n appending data\n')
    samples_dist = [samples_dist;dist];
    samples_speed = [samples_speed;speed];
    samples_class = [samples_class;class];
    samples_fr = [samples_fr;fr];
    samples_xpos = [samples_xpos;xpos];
    samples_t = [samples_t;t];
    samples_lane = [samples_lane;lane];
    samples_fcons = [samples_fcons;fcons];
    toc
end

%% save the output file to .mat data file
savefolder = '../Data/Data_for_Figures/';
save([savefolder 'samples_for_distance_analysis_' char(num2str(day_to_process)) '.mat'],'samples_*','-v7.3')

end

%%%%%%%%%%%%%%%%%%%%%% auxiliary function 2
function [all_d, all_v, all_vc, all_fr, all_x, all_t, all_lane, all_fcons] = stats_to_av_dist(data_v22,dayDate)

    Max_Dist = 1000; %maximum distance from an engaged AV for samples to be collected

    % initialize with appropriate data type to optimize for memory use  
    all_v = inf*ones(size(vertcat(data_v22.timestamp)));
    all_fr = inf*ones(size(vertcat(data_v22.timestamp)));
    all_d = inf*ones(size(vertcat(data_v22.timestamp)));

    all_x = int16(inf*ones(size(vertcat(data_v22.timestamp))));
    
    
    
    all_vc = uint8(all_x);
    all_lane = uint8(all_x);

    all_t = uint16(all_x);
    

    jj1 = 1;

    for i=1:length(data_v22)
        veh = data_v22(i);
        if veh.direction>0; continue;end % only aggregate for westbound traffic 

        nemp_p = ~isempty(veh.distance_to_downstream_engaged_av_meters);
        nemp_n = ~isempty(veh.distance_to_upstream_engaged_av_meters);

        if ((nemp_p)||(nemp_n))

            for ii =1:length(veh.timestamp)
                
               
                if nemp_p &&(~isnan(veh.distance_to_downstream_engaged_av_meters(ii))) && ...   %aggregate data for vehicles that have an AV downstream within Max_Dist
                        (veh.distance_to_downstream_engaged_av_meters(ii) <= Max_Dist) 

                    temp_d = veh.distance_to_downstream_engaged_av_meters(ii);

                    all_v(jj1) = veh.speed_meters_per_second(ii);
                    all_fr(jj1) = veh.fuel_rate_grams_per_second(ii);

                    all_vc(jj1) = uint8(veh.coarse_vehicle_class);
                    all_x(jj1) = int16(veh.x_position_meters(ii));
                    all_d(jj1) =   temp_d;
                    all_lane(jj1) = veh.lane_number;
                    all_t(jj1) = uint16(veh.timestamp(ii)-(1668772800-(18-dayDate)*24*60*60)); %record time in (s) after 6am of each test day
                    jj1 =jj1+1;
                end
                
                if nemp_n && (~isnan(veh.distance_to_upstream_engaged_av_meters(ii))) && ...  %aggregate data for vehicles that have an AV upstream within Max_Dist
                        (veh.distance_to_upstream_engaged_av_meters(ii) >= -Max_Dist) 

                    temp_d = veh.distance_to_upstream_engaged_av_meters(ii);

                        all_v(jj1) = veh.speed_meters_per_second(ii);
                        all_fr(jj1) = veh.fuel_rate_grams_per_second(ii);

                        all_vc(jj1) = uint8(veh.coarse_vehicle_class);
                        all_x(jj1) = int16(veh.x_position_meters(ii));
                        all_d(jj1) =   temp_d;
                        all_lane(jj1) = veh.lane_number;

                        all_t(jj1) = uint16(veh.timestamp(ii)-(1668772800-(18-dayDate)*24*60*60)); %record time in (s) after 6am of each test day
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
    all_fcons = all_fr./ (1e-6+all_v);
    
end

toc