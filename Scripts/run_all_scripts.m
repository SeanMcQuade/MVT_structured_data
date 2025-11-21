% Run all scripts to generate a slim and full version of the MVT data for the three 
% test days (16-18/11/2022), aggregate samples from the data and generate
% the results figures. 
% (C) 2025 CIRCLES - Sulaiman Almatrudi
clear
% Run through each of the test days
for processingDay = 16:18
    % Assemble GPS data into run trajectories
    fprintf(['running function assemble_data_GPS.m for day ' num2str(processingDay)...
        '-11-2022\n'])
    assemble_data_GPS(processingDay)
    % Generate full version of the data for a full test day
    fprintf(['running function generate_data_mvt_full.m for day ' num2str(processingDay)...
        '-11-2022\n'])
    generate_data_mvt_full(processingDay)
    % Generate slim version of the data for a full test day
    fprintf(['running function generate_data_mvt_slim.m for day ' num2str(processingDay)...
        '-11-2022\n'])
    generate_data_mvt_slim(processingDay)
    % Generate samples for analysis
    fprintf(['running function generate_data_samples.m for day ' num2str(processingDay)...
        '-11-2022\n'])
    generate_data_samples(processingDay)
    % Generate macroscopic fields from MOTION data
    fprintf(['running function generate_macroscopic_fields.m for day ' num2str(processingDay)...
        '-11-2022\n'])
    generate_macroscopic_fields(processingDay)
    % Plot macroscopic fields 
    fprintf(['running function plot_macroscopic_fields.m for day ' num2str(processingDay)...
        '-11-2022\n'])
    plot_macroscopic_fields(processingDay)
    % Plot microscopic trajectories 
    fprintf(['running function plot_microscopic_trajectories.m for day ' num2str(processingDay)...
        '-11-2022\n'])
    plot_microscopic_trajectories(processingDay)
end
% Plot stats of contorlled AVs and sampled leaders/followers and compare
% days
fprintf('running function plot_AV_stats.m for \n')
plot_AV_analysis()
