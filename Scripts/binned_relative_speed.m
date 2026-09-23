function [] = binned_relative_speed(filtered_dist_all_files, ...
    filtered_speed_all_files, day, outputFile)
% BINNED_RELATIVE_SPEED  Mean/median relative speed in distance bins behind an AV.
%
% Purpose
%   Produces the "Relative speeds behind AV" panel of the paper from the
%   relative speeds relative_speed_histogram has already computed for one day.
%
% Inputs
%   filtered_dist_all_files   relative distance to the nearest engaged AV (m),
%                             pooled over the day's segments
%   filtered_speed_all_files  matching relative speeds (m/s)
%   day                       16, 17, or 18 (November 2022); titles only
%   outputFile                absolute path of the .pdf to write
%
% Outputs
%   (none; writes outputFile)
%
% Notes
%   This was a script that ran inside relative_speed_histogram's workspace and
%   silently inherited d, filtered_*, day, j_start and j_end from it. Taking
%   them as arguments is what lets the caller be a function, and lets
%   mvt.sources see the dependency.
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

MAXDIST = 350; % (m) maximum relative distance to an ego vehicle to be considered
XWINDOW = [MAXDIST-122 7000]; % (m) , [West to East] set to -1 to include the full testbed.
AVLOCBUFFER = 1; % (m) buffer around 0 distance from av within which data is not considered.
TWINDOW = [0645, 0915]; % time window of interest in military time
% Output figure canvas in pixels. exportgraphics crops the PDF to the drawn
% content, and the content extent follows the figure size - which MATLAB
% otherwise derives from the screen, so the same code produced a different page
% under `matlab -batch` than in an interactive session. Pinning it makes the
% figure reproducible; the value is the batch default this was rendering at.
figRes = [700 420];
figScale = '-r384';
plotCol = [0,0,1;1,0,0;0,.6,0]; %Plotting colors
% initiate grid of bin edges
distToA = round(linspace(-MAXDIST,MAXDIST,ceil(MAXDIST*2/10)+1));
distToA = [distToA(1:floor(length(distToA)/2)),-AVLOCBUFFER,AVLOCBUFFER,...
    distToA(floor(length(distToA)/2)+2:end)];
binCenters = (distToA(1:end-1)+distToA(2:end))/2;

%used to put data into bins by distance;
d = filtered_dist_all_files; 
for binInd = 1:length(binCenters) %bin the indices by distance from AV
    binnedInd{binInd}   = find(distToA(binInd) <= d & d < distToA(binInd+1));
    binnedSpeed{binInd} = filtered_speed_all_files(binnedInd{binInd});
    %calculate mean speed in bin
    mean_binned_speed(binInd) = mean(binnedSpeed{binInd},"omitnan");
    median_binned_speed(binInd) = prctile(binnedSpeed{binInd},50);
end

fontSize = 12;

%plot histogram
axX = 350;
axY = [0.06,0.13];
fig = figure('Visible', 'off', 'Position', [10 50 figRes], ...
    'PaperPositionMode', 'auto');
closeFigure = onCleanup(@() close(fig));
hold on
plot(binCenters, mean_binned_speed,"LineWidth", 2);
plot(binCenters, median_binned_speed,"LineWidth", 2,"LineStyle", "--");
plot([min(binCenters), max(binCenters)],[0,0],"LineStyle",":")
plot([0,0],[-0.1,0.5],"LineStyle",":")
text(axX*.03,axY(1)+diff(axY)*.90,'behind AV','FontSize',fontSize,...
    'HorizontalAlignment','left','VerticalAlignment','bottom')
% text(-axX*.03,axY(1)+diff(axY)*.90,'ahead of AV','FontSize',fontSize,...
%     'HorizontalAlignment','right','VerticalAlignment','bottom')
legend("mean relative speed (m/s)", "median relative speed (m/s)", "","","FontSize", fontSize,"Location","northwest");
if day == 16 %Write the day in the title
    formatSpec = "Relative speed histogram Wednesday 16-Nov-2022";
    title_string = sprintf(formatSpec);
elseif day == 17
    formatSpec = "Relative speed histogram Thursday 17-Nov-2022";
    title_string = sprintf(formatSpec);
elseif day == 18
    formatSpec = "Relative speed histogram Friday 18-Nov-2022";
    title_string = sprintf(formatSpec);
end

title(title_string,"FontSize", fontSize)
xlabel('Distance behind nearest AV in same lane (m)','FontSize',fontSize)
ylabel('Statistics calculated per bin','FontSize',fontSize)
xlim([0 max(binCenters)]);
% xlim([min(binCenters) max(binCenters)]);

fprintf('Save figure in %s ...', outputFile), tic
exportgraphics(fig, outputFile, 'ContentType', 'vector');
fprintf(' Done (%0.0fsec).\n', toc)
end
