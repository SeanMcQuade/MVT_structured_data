% Run all scripts to generate a slim version of the MVT data for the three 
% test days (16-18/11/2022), aggregate samples from the data and generate
% the results figures. 
% (C) 2025 CIRCLES - Sulaiman Almatrudi

clear
global DAY_TO_PROCESS
% Run through each of the test days

for DAY_TO_PROCESS = 16:18
    
    % Assemble GPS data into run trajectories
    fprintf(['running script assemble_data_GPS.m for day ' num2str(DAY_TO_PROCESS)...
        '-11-2022\n'])
    assemble_data_GPS
    
    % Generate slim version of the data for a full test day
    fprintf(['running script generate_data_mvt_slim.m for day ' num2str(DAY_TO_PROCESS)...
        '-11-2022\n'])
    generate_data_mvt_slim
    
    % Generate samples for analysis
    fprintf(['running script generate_data_samples.m for day ' num2str(DAY_TO_PROCESS)...
        '-11-2022\n'])
    generate_data_samples
    
    % Generate macroscopic fields from MOTION data
    fprintf(['running script generate_macroscopic_fields_motion.m for day ' num2str(DAY_TO_PROCESS)...
        '-11-2022\n'])
    generate_macroscopic_fields_motion
    
    % Plot macroscopic fields 
    fprintf(['running script plot_macroscopic_fields.m for day ' num2str(DAY_TO_PROCESS)...
        '-11-2022\n'])
    plot_macroscopic_fields
    
end

% Plot stats of contorlled AVs and sampled leaders/followers
fprintf('running script plot_AV_stats.m for \n')
plot_AV_stats

