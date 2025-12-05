% Produces macroscopic fields and fuel results for a subset of data, to
% validate processes are part of peer review.
% (C) 2025 CIRCLES - Sulaiman Almatrudi, Jonathan Sprinkle
clear
% Run through each of the test days
% the peer-review example processes only 11-17-2022 (Thursday)
for processingDay = 17
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
% close other plots
close all
% Plot stats of contorlled AVs and sampled leaders/followers and compare
% days
fprintf('running function plot_AV_stats.m for \n')
plot_AV_analysis(17)
