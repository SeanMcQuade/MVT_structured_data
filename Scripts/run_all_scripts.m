% Run all scripts to generate a slim version of the MVT data for the three 
% test days (16-18/11/2022), aggregate samples from the data and generate
% the results figures. 
% (C) 2025 CIRCLES - Sulaiman Almatrudi

clear
global DAY_TO_PROCESS
% Run through each of the test days

for DAY_TO_PROCESS = 16:18
    
    % Assemble GPS data into run trajectories
    assemble_data_GPS
    
    % Generate slim version of the data for a full test day
    generate_data_mvt_slim
    
    % Generate samples for analysis
    generate_data_samples
    
    % Generate macroscopic fields from MOTION data
    generate_macroscopic_fields_motion
    
    % Plot macroscopic fields 
    plot_macroscopic_fields
    
    % Plot stats of contorlled AVs and sampled leaders/followers
    plot_AV_stats



end

