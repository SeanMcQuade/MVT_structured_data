%Plot average vehicle speed according to relative distance bins.
%figure;
%scatter(all_av_dist,rel_speed_all_files,0.2)

MAXDIST = 350; % (m) maximum relative distance to an ego vehicle to be considered
XWINDOW = [MAXDIST-122 7000]; % (m) , [West to East] set to -1 to include the full testbed.
AVLOCBUFFER = 1; % (m) buffer around 0 distance from av within which data is not considered.
TWINDOW = [0645, 0915]; % time window of interest in military time
% Output figure resolution  and scale
figRes = [1600 1000]; 
figScale = '-r384';
plotCol = [0,0,1;1,0,0;0,.6,0]; %Plotting colors
% Locate data to load
[parentDirectory, ~, ~] = fileparts(pwd);
[dataRootDirectory, ~, ~] = fileparts(parentDirectory);
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

fontSize = 18;

%plot histogram
axX = 350;
axY = [0.06,0.13];
figure;
hold on
plot(binCenters, mean_binned_speed,"LineWidth", 2);
plot(binCenters, median_binned_speed,"LineWidth", 2,"LineStyle", "--");
plot([min(binCenters), max(binCenters)],[0,0],"LineStyle",":")
plot([0,0],[-0.1,0.5],"LineStyle",":")
text(axX*.03,axY(1)+diff(axY)*.90,'behind AV','FontSize',fontSize,...
    'HorizontalAlignment','left','VerticalAlignment','bottom')
text(-axX*.03,axY(1)+diff(axY)*.90,'ahead of AV','FontSize',fontSize,...
    'HorizontalAlignment','right','VerticalAlignment','bottom')
legend("mean relative speed (m/s)", "median relative speed (m/s)", "","","FontSize", fontSize,"Location","northwest");
if day == 16 %Write the day in the title
    formatSpec = "Relative speed histogram Wed";
    title_string = sprintf(formatSpec);
elseif day == 17
    formatSpec = "Relative speed histogram Thurs";
    title_string = sprintf(formatSpec);
elseif day == 18
    formatSpec = "Relative speed histogram Fri";
    title_string = sprintf(formatSpec);
end

titlestring = sprintf(title_string,day,j_start, j_end);
title(titlestring,"FontSize", fontSize)
xlabel('Distance to nearest AV in same lane (m)','FontSize',fontSize)
ylabel('Statistics calculated per bin','FontSize',fontSize)
xlim([min(binCenters) max(binCenters)]);

formatSpecSave = "../../results/figures/Relative speeds behind AV, day %d.fig";
savename = sprintf(formatSpecSave,day);
saveas(gcf,savename)