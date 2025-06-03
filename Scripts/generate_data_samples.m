% Run to generate data samples in .mat format for generating fuel consumption figures.
% (C) 2025 by Sulaiman Almatrudi for CIRCLES energy team; 
% further adapted by Sean McQuade for CIRCLES scenario team.
% This is licensed under BSD-3 clause license: https://opensource.org/license/bsd-3-clause

clearvars -except DAY_TO_PROCESS 
if ~exist('DAY_TO_PROCESS', 'var') || isempty(DAY_TO_PROCESS) || ~ismember( DAY_TO_PROCESS,[16:18])
    global DAY_TO_PROCESS
    DAY_TO_PROCESS = input(['Enter the day of Nov. 2022 to collect samples from '...
        'MVT data files (from 16 to 18): ']);
end
%========================================================================
% Generate samples_for_distance_analysis
%========================================================================
[parentDirectory, ~, ~] = fileparts(pwd);
dataFolder = fullfile(parentDirectory,'Data',['Data_2022-11-' char(num2str(DAY_TO_PROCESS)) '__MVT_Slim']);
I24FilesInDir = dir(fullfile(dataFolder, ['I-24*' char(num2str(DAY_TO_PROCESS)) '*.json']));
nrFiles = length(I24FilesInDir);
if nrFiles < 24
    dataFolder = fullfile(parentDirectory,'Data',['Data_2022-11-' char(num2str(DAY_TO_PROCESS)) '__MVT_Full']);
    I24FilesInDir = dir(fullfile(dataFolder, ['I-24*' char(num2str(DAY_TO_PROCESS)) '*.json']));
    nrFiles = length(I24FilesInDir);
    if nrFiles < 24
        error(['Processed I-24 data files not found. Please generate files'...
        ' by running generate_data_mvt_slim.m first.'])
    end
end

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
for file_nr =1:nrFiles
    fprintf(sprintf('Loading original data file %d / %d ...', file_nr,nrFiles ))
    tic
    filenameLoad = fullfile(dataFolder , I24FilesInDir(file_nr).name);
    data = jsondecode(fileread(filenameLoad));
    fprintf('Done (%0.0fsec).\n',toc)
    
    if ~isfield(data,'distance_to_upstream_engaged_av_meters');fprintf('no AVs on the road\n');continue;end
    
    fprintf('Processing data...')
    tic

    [dist, speed, class, fr, xpos, t, lane, fcons] = stats_to_av_dist(data,DAY_TO_PROCESS);
    fprintf('Done (%0.0fsec).\n',toc)

    fprintf('Appending data...');tic
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

%% save the output file to .mat data file
saveFolder = fullfile(parentDirectory,'Data','Data_for_Figures');
save(fullfile(saveFolder, ['samples_for_distance_analysis_' char(num2str(DAY_TO_PROCESS)) ...
    '.mat']),'samples_*','-v7.3')


%%%%%%%%%%%%%%%%%%%%%% auxiliary function %%%%%%%%%%%%%%%%%%%%%%
function [all_d, all_v, all_vc, all_fr, all_x, all_t, all_lane, all_fcons] = ...
    stats_to_av_dist(dataTemp,dayDate)

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
        if isfield(veh,'direction') % only aggregate for westbound traffic 
            if veh.direction>0 
                continue
            end 
        end

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

                    all_vc(jj1) = uint8(veh.coarse_vehicle_class);
                    all_x(jj1) = int16(veh.x_position_meters(ii));
                    all_d(jj1) =   temp_d;
                    all_lane(jj1) = veh.lane_number;
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

                        all_vc(jj1) = uint8(veh.coarse_vehicle_class);
                        all_x(jj1) = int16(veh.x_position_meters(ii));
                        all_d(jj1) =   temp_d;
                        all_lane(jj1) = veh.lane_number;
                        %record time in (s) after 6am of each test day
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
    all_fcons = all_fr./ (1e-6+all_v);
    
end
