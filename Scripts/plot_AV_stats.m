% script to plot multiple statistics for AV fuel consumption for experiments
% run on the following dates: 11/16/2022, 11/17/2022, and 11/18/2022
% (C) 2025 Sean McQuade and Sulaiman Almatrudi
% This script uses data in ../Data/Data_for_Figures.
% This licensed under BSD-3 clause license:
% https://opensource.org/license/bsd-3-clause

clear
statsToPlotChoices = ["effective","mean","median"];
MAXDIST = 350; % (m) maximum relative distance to an ego vehicle to be considered
XWINDOW = [MAXDIST-122 7000]; % (m) , [West to East] set to -1 to include the full testbed.
AVLOCBUFFER = 1; % (m) buffer around 0 distance from av within which data is not considered.
TWINDOW = [0645, 0915]; % time window of intrest in military time
flagSave = true; % save figures
minSpdRaw = 0; %(m/s) minimum speed to filter for all samples
minSpdMean = 1; % (m/s) minimum speed for samples when calculating the mean fuel consumption

figRes = [1600 1000];
figScale = '-r384';
plotCol = [0,0,1;1,0,0;0,.6,0]; %Plotting colors

[parentDirectory, ~, ~] = fileparts(pwd);

% intiate grid of bin edges
distToA = round(linspace(-MAXDIST,MAXDIST,ceil(MAXDIST*2/10)+1));
distToA = [distToA(1:floor(length(distToA)/2)),-AVLOCBUFFER,AVLOCBUFFER,...
    distToA(floor(length(distToA)/2)+2:end)];
binCenters = (distToA(1:end-1)+distToA(2:end))/2;

maskingDist = 30; % bins closer than this num meters are masked
indUnmaskedLeft = find(distToA(2:end)<=-maskingDist);
indMaskedLeft = [indUnmaskedLeft(end) find(-maskingDist<distToA(2:end) &...
    distToA(2:end)<0)];
indUnmaskedRight = find( maskingDist <= distToA(1:end-1));
indMaskedRight = [find(0<distToA(1:end-1) & distToA(1:end-1) < maskingDist),...
    indUnmaskedRight(1)];

indicesPartitions{1} = indUnmaskedLeft;
indicesPartitions{2} = indMaskedLeft;
indicesPartitions{3} = indMaskedRight;
indicesPartitions{4} = indUnmaskedRight;

testDays = [16, 17, 18];
fprintf('Loading aggregated data samples... ') ; tic
for i=1:length(testDays) % three days, Wed thrus fri
    loadedData = load(fullfile(parentDirectory , 'Data', 'Data_for_Figures',...
        ['samples_for_distance_analysis_' char(num2str(testDays(i))) '.mat']));
    
    % filter data based on x position, time, and speed
    ind = loadedData.samples_speed >= minSpdRaw; %set minSpdRaw to '0' to keep all speed samples
    
    if length(XWINDOW)>1
        ind = loadedData.samples_xpos >= XWINDOW(1) & loadedData.samples_xpos <= XWINDOW(2);
    end
    
    if length(TWINDOW)>1
        tWinHr = floor(TWINDOW/100) ;
        tWinMin = mod(TWINDOW,100) ;
        tWinSec = ((tWinHr-6)*60 + tWinMin)*60;
        ind = ind & loadedData.samples_t >= tWinSec(1) & ...
            loadedData.samples_t <= tWinSec(2);
    end
    
    samplesClassCell{i} = loadedData.samples_class(ind);
    loadedData = rmfield(loadedData,'samples_class');
    samplesDistCell{i} = loadedData.samples_dist(ind);
    loadedData = rmfield(loadedData,'samples_dist');
    samplesFconsCell{i} = loadedData.samples_fcons(ind);
    loadedData = rmfield(loadedData,'samples_fcons');
    samplesFrCell{i} = loadedData.samples_fr(ind);
    loadedData = rmfield(loadedData,'samples_fr');
    samplesLaneCell{i} = loadedData.samples_lane(ind);
    loadedData = rmfield(loadedData,'samples_lane');
    samplesSpeedCell{i} = loadedData.samples_speed(ind);
    loadedData = rmfield(loadedData,'samples_speed');
    samplesTCell{i} = loadedData.samples_t(ind);
    loadedData = rmfield(loadedData,'samples_t');
    samplesXposCell{i} = loadedData.samples_xpos(ind);
    loadedData = rmfield(loadedData,'samples_xpos');
end
clear loadedData
fprintf('Done (%0.0fsec).\n',toc)

%%%%%%%%%%%% preallocate statistics %%%%%%%%%%%%%%%
fprintf('Binning data to distance from active AVs... ') ; tic
meanFcons   = zeros(length(testDays),length(binCenters));
binnedFcons = cell(length(testDays),length(binCenters));
binnedFr    = cell(length(testDays),length(binCenters));
binnedSpeed = cell(length(testDays),length(binCenters));
effectiveFc = zeros(length(testDays),length(binCenters));
medianFcons = zeros(length(testDays),length(binCenters));
binCounts   = zeros(length(testDays),length(binCenters));
binnedInd   = cell(1,length(binCenters));

for i=1:length(testDays)
    d = samplesDistCell{i}; %used to put data into bins by distance
    for j = 1:length(binCenters) %bin the indices by distance from AV
        binnedInd{j}     = find(distToA(j) <= d & d < distToA(j+1));
        binnedSpeed{i,j} = samplesSpeedCell{i}(binnedInd{j});
        binnedFr{i,j}    = samplesFrCell{i}(binnedInd{j});
        binnedFcons{i,j} = samplesFconsCell{i}(binnedInd{j});
        effectiveFc(i,j) = sum(binnedFr{i,j})/max(sum(binnedSpeed{i,j}),1e-6);
        meanFcons(i,j)   = mean(binnedFcons{i,j}(binnedSpeed{i,j}>=minSpdMean));
        medianFcons(i,j) = median(binnedFcons{i,j});
        binCounts(i,j)   = length(binnedInd{j});
    end
end
fprintf('Done (%0.0fsec).\n',toc)

% Create triplets of mean median effictive
% Partition them into Left, Masked_left, Masked_right, Right
% Calculate regression

fprintf('Processing data samples... ') ; tic
statsToPlot         = cell(length(testDays),1);
partitionedStats    = cell(length(testDays),4);
dLeft               = zeros(length(testDays),length(indUnmaskedLeft));
dRight              = zeros(length(testDays),length(indUnmaskedRight));
dMaskedLeft         = zeros(length(testDays),length(indMaskedLeft));
dMaskedRight        = zeros(length(testDays),length(indMaskedRight));
valLeft             = zeros(length(testDays),length(indUnmaskedLeft));
valRight            = zeros(length(testDays),length(indUnmaskedRight));
valMaskedLeft       = zeros(length(testDays),length(indMaskedLeft));
valMaskedRight      = zeros(length(testDays),length(indMaskedRight));
wRegLeft            = cell(length(testDays),1);
wRegRight           = cell(length(testDays),1);
sRegression         = cell(length(testDays),1);

for i=1:length(testDays)
    statsToPlot{i} = [meanFcons(i,:); effectiveFc(i,:); medianFcons(i,:)];
    % partitioned_stats{1}left, {2}masked left, {3}masked right, {4}right.
    % Each has three rows [mean_fcons; effective_fc; median_fcons]
    [partitionedSingleDay] = mask_nearby_points(statsToPlot{i},indicesPartitions);
    for k = 1:3
        wRegLeft{i}(k,:) = find_regression_line(binCenters(indicesPartitions{1}).', ...
            partitionedSingleDay{1}(k,:).', diag(binCounts(i,indicesPartitions{1})));
        wRegRight{i}(k,:) = find_regression_line(binCenters(indicesPartitions{4}).',...
            partitionedSingleDay{4}(k,:).', diag(binCounts(i,indicesPartitions{4})));
        yValsLeft  = partitionedSingleDay{1}(k,:).';
        yValsRight = partitionedSingleDay{4}(k,:).';
        ratioBehindVsAhead =  -(yValsRight./yValsLeft(end:-1:1,:)-1)*100;
        sRegression{i}(k,:) = find_regression_line(binCenters(indicesPartitions{4}).', ...
            ratioBehindVsAhead);
    end
    [partitionedStats{i,1},partitionedStats{i,2},partitionedStats{i,3},...
        partitionedStats{i,4}] = partitionedSingleDay{1:4};
end
fprintf('Done (%0.0fsec).\n',toc)
% Calculate Mean Bulk speed and AVs for Wed. Thurs. and Fri.
strArray = ["fields_motion_2022-11-16" "fields_motion_2022-11-17" ...
    "fields_motion_2022-11-18"];

if length(TWINDOW) > 1
    TWINDOWSTR = num2str(TWINDOW);
else
    TWINDOW = [0600,1000];
    TWINDOWSTR = num2str(TWINDOW);
end
%%

fprintf('Loading and processing GPS and field data... ') ; tic
uTempMean = [];
nAVsActive = [];
for i=1:3 %wed. thurs. fri.
    fieldFilename = fullfile(parentDirectory,'Data' ,'Data_for_Figures',...
        char(strArray(i)));
    load(fieldFilename)
    GPSFilename = fullfile(parentDirectory, 'Data',...
        'Data_GPS',['CIRCLES_GPS_10Hz_2022-11-' num2str(15+i) '.json']);
    data = jsondecode(fileread(GPSFilename));
    tStep = floor(mean(diff(t)));
    
    for ii = 1:length(t)
        tTemp = t(ii);
        uTempV = field.value.U(ii,:)*field.factor.U;
        rhoTempV = field.value.Rho(ii,:)*field.factor.Rho;
        rhoTempV = rhoTempV(~isnan(uTempV));
        uTempV = uTempV(~isnan(uTempV));
        
        totalRhoTemp = sum(rhoTempV);
        uTempMean(i,ii) = sum(uTempV .* rhoTempV)/totalRhoTemp;
        
        avsOnTheRoad = data([data.direction] == -1 & [data.first_timestamp] ...
            <= tTemp & [data.last_timestamp] >= tTemp + tStep);
        nAvsActiveTemp = 0;
        nAvsTemp = 0;
        for ii2 = 1:length(avsOnTheRoad)
            av = avsOnTheRoad(ii2);
            xMin = -122;
            xMax = 7000;
            % count only AVs in xwindow
            avInCellInd = av.timestamp >= tTemp & av.timestamp <= tTemp...
                + tStep & av.x_position >= xMin & av.x_position <= xMax;
            if isfield(av,"control_car") && ~isempty(av.control_car)
                isAvControl = av.control_car(avInCellInd);
            else
                isAvControl = av.controller_engaged(avInCellInd);
            end
            % the AV is considered active if it has an active signal
            % for more than half the time bin (5s default)
            isAvControl = (sum(isAvControl)/length(isAvControl))>= 0.5;
            if isAvControl
                nAvsActiveTemp = nAvsActiveTemp + 1;
            end
            nAvsTemp = nAvsTemp +1;
        end
        nAVsActive(i,ii) = nAvsActiveTemp;
        nAVs(i,ii) = nAvsTemp;
    end
    % Plotting bulk speed
    LIGHTBLUE = 1/2*[0,0,1]+1/2*[1,1,1];
    LIGHTRED = 1/2*[1,0,0]+1/2*[1,1,1];
    FONTSIZE = 14;
    tCST(i,:) = datetime(t, 'convertfrom', 'posixtime', 'Format', 'HH:mm:ss.SSS',...
        'TimeZone' ,'America/Chicago') ;
    
    
end
fprintf('Done (%0.0fsec).\n',toc)
%%
% Plot in 3x3 grid
fprintf('Plotting stats... ') ; tic
figure(1);
clf
timeString = [num2str(firstHour), ':', num2str(firstMinutes), ' to ',...
    num2str(secondHour), ':', num2str(secondMinutes)];
dateStrAry = {['Wed Nov 16, time ' timeString],['Thu Nov 17, time ' timeString],...
    ['Fri Nov 18, time ' timeString]};
if length(statsToPlotChoices) == 3
    figsToProduce = [2,3];
else
    figsToProduce = 2 ;
end
for figInd = figsToProduce
    if figInd == 2
        statsToPlotChoice = statsToPlotChoices(1);
    else
        statsToPlotChoice = statsToPlotChoices;
    end
    
    for i=1:length(testDays)
        firstHour = floor(TWINDOW(1)/100); firstMinutes = mod(TWINDOW(1),100);
        secondHour = floor(TWINDOW(2)/100); secondMinutes = mod(TWINDOW(2),100);
        
        plotStartT = datetime(2022,11,testDays(i), firstHour, firstMinutes ,0,...
            'TimeZone' ,'America/Chicago');
        plotEndT = datetime(2022,11,testDays(i), secondHour, secondMinutes,0, ...
            'TimeZone' ,'America/Chicago');
        tWindowNAvsActive = nAVsActive(i,tCST(i,:)>plotStartT);
        
        subplot(3,3,3*i-2)
        hold on
        yyaxis right
        hold on %plot number of AVs
        plot(tCST(i,tCST(i,:) <= plotStartT), nAVsActive(i,tCST(i,:) <= plotStartT),...
            'Color',LIGHTRED,'LineStyle','-','LineWidth',1.5)
        plot(tCST(i,tCST(i,:) >= plotEndT), nAVsActive(i,tCST(i,:) >= plotEndT),...
            'Color',LIGHTRED,'LineStyle','-','LineWidth',1.5)
        plot(tCST(i,tCST(i,:) >= plotStartT & tCST(i,:) <= plotEndT), ...
            nAVsActive(i,tCST(i,:) >= plotStartT & tCST(i,:) <= plotEndT),'r-','LineWidth',1.5)
        plot([plotStartT plotStartT],[0,55],'k:','LineWidth',1.5) %vertical lines
        plot([plotEndT plotEndT],[0,55],'k:','LineWidth',1.5)
        ylabel('Number of active AVs')
        set(gca,'YColor','r');
        grid on
        axx = [plotStartT - minutes(30), plotEndT + minutes(30)];
        axy = [0,55];
        title('Traffic state and AV deployment vs. time')
        text(plotStartT + minutes(2),(axy(2) - axy(1))*0.9,'start time','FontSize',FONTSIZE)
        text(plotEndT - minutes(47),(axy(2) - axy(1))*0.9,'end time','FontSize',FONTSIZE)
        set(gca,'FontSize',FONTSIZE)
        xlim(axx)
        ylim(axy)

        yyaxis left
        plot(tCST(i,tCST(i,:) <= plotStartT), uTempMean(i,tCST(i,:) <= plotStartT),'Color',...
            LIGHTBLUE,'LineStyle','-','LineWidth',2)
        plot(tCST(i,tCST(i,:) >= plotEndT), uTempMean(i,tCST(i,:) >= plotEndT),'Color',...
            LIGHTBLUE,'LineStyle','-','LineWidth',2)
        plot(tCST(i,tCST(i,:) >= plotStartT & tCST(i,:) <= plotEndT),...
            uTempMean(i,tCST(i,:) >= plotStartT & tCST(i,:) <= plotEndT),'Color','b',...
            'LineStyle','-','LineWidth',2)
        ylabel('Mean speed (m/s)')
        ylim([0,32])
        meanBulkSpeedInWindow(i) = mean(uTempMean(i,tCST(i,:) >= plotStartT & tCST(i,:) <= plotEndT));
        hold off
        hold on
        dateStr = dateStrAry{i};
        [partitionedStatsPerDay{1},partitionedStatsPerDay{2},partitionedStatsPerDay{3},...
            partitionedStatsPerDay{4}] = partitionedStats{i,:};
        subplot(3,3,3*i-1)
        cla
        plot_two_sided(dateStr,plotCol,partitionedStatsPerDay,wRegLeft{i},...
            wRegRight{i}, binCenters, indicesPartitions,statsToPlotChoice)
        subplot(3,3,3*i)
        cla
        plot_one_sided(dateStr,plotCol,partitionedStatsPerDay,sRegression{i},...
            binCenters, indicesPartitions,statsToPlotChoice)
    end
    
    
    % Adjust subfig placements and size
    posX = [.05,.355,.76];
    posY = [.74,0,.07]; posY(2) = (posY(1)+posY(3))/2;
    posS = [.21,.22]; posW = [.34, posS(2)];
    subplot(3,3,3), set(gca,'Position',[posX(3),posY(1),posS])
    subplot(3,3,6), set(gca,'Position',[posX(3),posY(2),posS])
    subplot(3,3,9), set(gca,'Position',[posX(3),posY(3),posS])
    
    subplot(3,3,1), set(gca,'Position',[posX(1),posY(1),posS])
    subplot(3,3,4), set(gca,'Position',[posX(1),posY(2),posS])
    subplot(3,3,7), set(gca,'Position',[posX(1),posY(3),posS])
    
    subplot(3,3,2), set(gca,'Position',[posX(2),posY(1),posW])
    subplot(3,3,5), set(gca,'Position',[posX(2),posY(2),posW])
    subplot(3,3,8), set(gca,'Position',[posX(2),posY(3),posW])
    fprintf('Done (%0.0fsec).\n',toc)
    
    if flagSave
        fname = fullfile(parentDirectory,'Figures', ['fig_' num2str(figInd)...
            '_fuel_results_',...
            char(strjoin(statsToPlotChoice(:),'_')), ...
            '_',char(num2str(TWINDOW(1))),'_',char(num2str(TWINDOW(2)))]);
        set(gcf,'Position',[10 40 figRes],'PaperPositionMode','auto')
        print(fname,'-dpng',figScale);
        savefig(fname)
    end
end
%%
figure(2)
clf
fontSize = 18;
axx = 350;
axy = [0, 1400000]; % [400000, 1400000]
equispacedBins = fix(binCenters);
equispacedBins(floor(length(equispacedBins)/2)+1) = [];
weekDaysStr = ["Wednesday","Thursday","Friday"];
for i =1:length(testDays)
    subplot(3,1,i)
    countsTemp = binCounts(i,:);
    countsTemp(floor(length(countsTemp)/2)+1) = [];
    bar(equispacedBins,countsTemp )
    hold on
    plot([0,0],axy,'k:','LineWidth',1.5)
    plot([-30,-30],axy,'r:','LineWidth',2.5)
    plot([30,30],axy,'r:','LineWidth',2.5)
    title(weekDaysStr(i),'FontSize',fontSize)
    if i ==2
        xlabel('Distance to nearest AV in same lane (m)','FontSize',fontSize)
        ylabel('number of data points','FontSize',fontSize)
    end
    
    axis([-axx,axx,axy])
    
    text(axx*.2,axy(1)+diff(axy)*.8,'behind AV','FontSize',fontSize,...
        'HorizontalAlignment','left','VerticalAlignment','bottom')
    text(-axx*.2,axy(1)+diff(axy)*.8,'ahead of AV','FontSize',fontSize,...
        'HorizontalAlignment','right','VerticalAlignment','bottom')
    hold off
    set(gca,'fontsize',18)
end
if flagSave
        fname = fullfile(parentDirectory,'Figures', ['fig_SM' num2str(2)...
            '_vehicle_samples_counts_',...
            char(strjoin(statsToPlotChoice(:),'_')), ...
            '_',char(num2str(TWINDOW(1))),'_',char(num2str(TWINDOW(2)))]);
        set(gcf,'Position',[10 40 figRes],'PaperPositionMode','auto')
        print(fname,'-dpng',figScale);
        savefig(fname)
end
fprintf('Done (%0.0fsec).\n',toc)


%========================================================================
function plot_two_sided(dateStr,col,partitionedStats,wRegLeft, wRegRight, ...
    binD, indicesPartitions,statsToPlotChoice)

axX = 350;
axY = [0.06,0.13];
FONTSIZE = 14;
colB = 1/5*col+4/5*[1,1,1]; % brightest version of color
colSB = 1/2*col+1/2*[1,1,1]; %slightly brighter version
colD = 1/2*col+1/2*[0,0,0]; % darker version of color
% color_ary = {colsb,col,cold};

avX1 = binD(indicesPartitions{1}); avY1 = partitionedStats{1};              %left_unmasked
avX2 = binD(:,indicesPartitions{2}); avY2 = partitionedStats{2};           %left_masked
avX3 = binD(:,indicesPartitions{3}); avY3 = partitionedStats{3};           %right_masked
avX4 = binD(:,indicesPartitions{4}); avY4 = partitionedStats{4};          %right_unmasked

plot([0,0],axY,'k:','LineWidth',1.5)
hold on
grid on
text(axX*.03,axY(1)+diff(axY)*.90,'behind AV','FontSize',FONTSIZE,...
    'HorizontalAlignment','left','VerticalAlignment','bottom')
text(-axX*.03,axY(1)+diff(axY)*.90,'ahead of AV','FontSize',FONTSIZE,...
    'HorizontalAlignment','right','VerticalAlignment','bottom')
statsAry = {'mean','effective','median'};
linestyleAry = {'-','-','-'};
linesizeAry = [2.5, 2.5, 2.5];

[ ~, indQuantityPlotted] = ismember(statsToPlotChoice,statsAry);% [mean, effective, median]
for i = indQuantityPlotted
    plot(avX1,wRegLeft(i,:),'--','Color',colD(i,:),'LineWidth',2.5);  %plot regression
    plot(avX4,wRegRight(i,:),'--','Color',colD(i,:),'LineWidth',2.5);
    plot(avX1,avY1(i,:) ,linestyleAry{i} ,'Color',col(i,:),'LineWidth',linesizeAry(i));
    plot(avX2,avY2(i,:) ,linestyleAry{i} ,'Color',colB(i,:),'LineWidth',2);
    plot(avX3,avY3(i,:) ,linestyleAry{i} ,'Color',colB(i,:),'LineWidth',2);
    hl(i) = plot(avX4,avY4(i,:) ,linestyleAry{i} ,'Color',col(i,:),'LineWidth',linesizeAry(i));
    title(['Fuel Consumption ',dateStr],'FontSize',FONTSIZE)
end
axis([-axX,axX,axY])
grid on
set(gca,'xtick',-300:100:300)
legend(hl(indQuantityPlotted),[statsAry(indQuantityPlotted) " FC"])
title(['Fuel consumption (FC). ',dateStr],'FontSize',FONTSIZE)
xlabel('Distance to nearest AV in same lane (m)','FontSize',FONTSIZE)
ylabel('Fuel consumption (g/m)','FontSize',FONTSIZE)
set(gca,'FontSize',FONTSIZE)
end

function plot_one_sided(dateStr,col,partitionedStats,sRegression, binD, ...
    indicesPartitions, statsToPlotChoice)
axX = 350;
axY = [-2,30];
FONTSIZE = 14;
colB = 1/4*col+3/4*[1,1,1]; % brightest version of color
colSB = 2/3*col+1/3*[1,1,1]; %slightly brighter version
colD = 1/2*col+1/2*[0,0,0]; % darker version of color
% color_ary = {colsb,col,cold};
linestyleAry = {'-','-','-'};
linesizeAry = [2.5, 2.5, 2.5];

avX1 = binD(indicesPartitions{1}); avY1 = partitionedStats{1};              %left_unmasked
avX2 = binD(:,indicesPartitions{2}); avY2 = partitionedStats{2};           %left_masked
avX3 = binD(:,indicesPartitions{3}); avY3 = partitionedStats{3};           %right_masked
avX4 = binD(:,indicesPartitions{4}); avY4 = partitionedStats{4};          %right_unmasked
plot([0,axX],[0,0],'k:','LineWidth',1.5)
hold on
grid on
statsAry = ["mean", "effective", "median"];
[~, indQuantityPlotted] = ismember(statsToPlotChoice,statsAry);% [mean, effective, median]
for i = indQuantityPlotted
    hl(i) = plot(avX4,-(avY4(i,:)./avY1(i,end:-1:1)-1)*100,linestyleAry{i},...
        'Color',col(i,:),'LineWidth',linesizeAry(i));
    plot(avX3,-(avY3(i,:)./avY2(i,end:-1:1)-1)*100,linestyleAry{i},'Color',...
        colB(i,:),'LineWidth',2)
    plot(avX4,sRegression(i,:),'--','Color',colD(i,:),'LineWidth',2.5);
end
axis([0,axX,axY])
axis ij
grid on
set(gca,'xtick',0:100:300)
ytickformat('percentage')
legend(hl(indQuantityPlotted),[statsAry(indQuantityPlotted) " FC"],'Location','southeast')
title(['Fuel reduction ahead vs. behind AV '],'FontSize',FONTSIZE) %changed from "behind vs ahead" 12/10/24
xlabel('Distance to nearest AV in same lane (m)','FontSize',FONTSIZE)
ylabel('Fuel consumption reduction','FontSize',FONTSIZE)
set(gca,'FontSize',FONTSIZE)
end

function [partitionedY] = mask_nearby_points(y,indicesPartitions)
partitionedY{1,1} = y(:,indicesPartitions{1});
partitionedY{1,2} = y(:,indicesPartitions{2});
partitionedY{1,3} = y(:,indicesPartitions{3});
partitionedY{1,4} = y(:,indicesPartitions{4});
end

function [regLine] = find_regression_line(xVals,yVals, weights)

aMatrix = [xVals(:) ones(size(xVals(:)))];
if nargin == 3
    solutionRatio = (transpose(aMatrix)*weights*aMatrix)\(transpose(aMatrix)*weights*yVals);
else
    solutionRatio = (transpose(aMatrix)*aMatrix)\(transpose(aMatrix)*yVals);
end
regLine = transpose(aMatrix*solutionRatio);

end