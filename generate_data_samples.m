% Script to plot graphs comparing statistics relative to distance from different ego vehicles (Active AVs, Ghost, MOTION counterfactuals)
% (C) 11/8/2024 by Sulaiman Almatrudi for CIRCLES energy team
tic

TestDAY = 16; % /NOV/2022 Enter 16, 17, or 18 here.

%========================================================================
% Generate samples_for_distance_analysis
%========================================================================

%% data preprocessing
%Vehicles near active AVs 
av_file_name = ['samples_for_distance_analysis_' char(num2str(TestDAY)) '.mat'];
% if isfile(av_file_name)
%     fprintf('\nLoading data relative to active AVs\n')
%     load(av_file_name)
% else
%     fprintf('\nGenerating data relative to active AVs\n')
%     generate_data_samples_v10(TestDAY);
%     load(av_file_name)
% end

fprintf('\nGenerating data relative to active AVs\n')
generate_data_samples_v10(TestDAY);


% load(av_file_name)
% 
% fprintf('\nCalculating stats to plot\n')
% [mean_v,perc_v,mean_fc,perc_fc,mean_fov,perc_fov,count_d,eff_fov] = ...
%     calculate_stats_to_plot(dist_to_a,percentiles,samples_dist,samples_speed,samples_fr,...
%     samples_xpos, samples_t,XWINDOW,TWINDOW_s,min_spd_fov);
% 
% stats_tp.AV.mean_v = mean_v;
% stats_tp.AV.perc_v = perc_v;
% stats_tp.AV.mean_fc = mean_fc;
% stats_tp.AV.perc_fc = perc_fc;
% stats_tp.AV.mean_fov = mean_fov;
% stats_tp.AV.perc_fov = perc_fov;
% stats_tp.AV.count_d = count_d;
% 
% AV_mean_fc = stats_tp.AV.mean_fc;

% scripts to aggregate samples form I24MOTION data relative to position
% from nearest engaged AV (Sulaiman Almatrudi) 

%% Functions used in this script 
%%%%%%%%%%%%%%%%%%%%%% auxiliary function 1

function [] = generate_data_samples_v10(day_to_process)

folder = [ '2022-11-' char(num2str(day_to_process)) '__MVT_Data_Slim'];
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
save(['samples_for_distance_analysis_' char(num2str(day_to_process)) '.mat'],'samples_*','-v7.3')

end

%%%%%%%%%%%%%%%%%%%%%% auxiliary function 2
% function   [mean_v,perc_v,mean_fc,perc_fc,mean_fov,perc_fov,count_d,eff_fov] = ...
%     calculate_stats_to_plot(dist_to_a,percentiles,all_d,all_v,all_fc,...
%     all_x, all_t, x_win, t_win, min_spd_fov)
% 
%     % Function to calculate the different stats to plot 
%     % (C) 2024/8/11 by Sulaiman Almatrudi for CIRCLES energy team
% 
%      fprintf('\n_______20%%_______40%%_______60%%_______80%%______100%%\n')
%     for i = 1: length(dist_to_a)-1
% 
%         if floor(i/length(dist_to_a)*50)>floor((i-1)/length(dist_to_a)*50) % progress marker
%            fprintf('^')
%         end
%         d = dist_to_a(i);
%         d_nxt =  dist_to_a(i+1);
% 
%         indx_d = all_d>=d & all_d<d_nxt &  all_v>= min_spd_fov;
% 
%         if length(x_win)>1
%             indx_d = indx_d & all_x >= x_win(1) & all_x <= x_win(2);
%         end
%         if length(t_win)>1
%             indx_d = indx_d & all_t >= t_win(1) & all_t <= t_win(2);
%         end
% 
% 
%         count_d(i) = sum(indx_d);
% 
%         v_d = all_v(indx_d);
%         fc_d = all_fc(indx_d);
% 
%         mean_v(i) = mean(v_d);
%         perc_v(i,:) = prctile(v_d,percentiles);
% 
%         mean_fc(i) = mean(fc_d);
%         perc_fc(i,:) = prctile(fc_d,percentiles);
% 
% 
%         fov = fc_d./(1e-6+v_d);
% 
%         min_spd_mean_fov_only = 5;
%         % here I can modify ... sum(fc_d)/sum(v_d) total fuel consumed / total dist;  sum(r*dt)/sum(v*dt); mean(fov)!!!!
%         mean_fov(i) = mean(fov(v_d >= min_spd_mean_fov_only)); %augmented by Benni
%         eff_fov(i) = sum(fc_d)/max(sum(v_d),1e-6); %Total fuel consumption, also grams per meter
%         % mean_fov(i) = mean(fov(v_d >= min_spd_fov)); %augmented by Benni
%         % mean_fov(i) = mean(fov); 
%         perc_fov(i,:) = prctile(fov,percentiles);
% 
%     end
% end
% % legend_mean_min_speed = sprintf('Rel. to AVs (GPS) minimum speed (mean only) = %d',min_spd_mean_fov_only);
% dist_to_AV_bins = [-dist_x(end:1)]
% %extra plotting
% dist_to_AV_bins = [-flip(dist_x) 0 dist_x];
% figure(1)
% hold on
% plot(dist_to_AV_bins,eff_fov,'.:','Color','r','MarkerSize',6,'LineWidth',2)
% legend('eff_fov')

%%%%%%%%%%%%%%%%%%%%%% auxiliary function 3
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
