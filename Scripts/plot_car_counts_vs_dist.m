% script to plot number of bulk traffic vehicles near AVs for experiments run on the following dates:
% 11/16/2022, 11/17/2022, and 11/18/2022 Authored by Sean McQuade. This script uses data in
% ../Data/Data_for_Figures. This licensed under BSD-3 clause license: https://opensource.org/license/bsd-3-clause

%% Load data from data_samples.mat
test_days = [16,17,18];

MAX_DIST = 350; % (m) maximum relative distance to an ego vehicle to be considered 
XWINDOW = [MAX_DIST-122 7000]; %[MAX_DIST, 6700-MAX_DIST]; % (m) , [West to East] set to -1 to include the full testbed. Find max and min (samples.xpos) and subtract 400m
TWINDOW = [0645, 0915];%; 0700, 0730; 0730, 0800]; %windows of intrest 12/5/2024

dist_to_a = round(linspace(-MAX_DIST,MAX_DIST,MAX_DIST*2/10+1));
dist_to_a = [dist_to_a(1:floor(length(dist_to_a)/2)),-1,1,dist_to_a(floor(length(dist_to_a)/2)+2:end)];
bin_d = (dist_to_a(1:end-1)+dist_to_a(2:end))/2;
flag_save = 1;
fig_res = [1600 1000];
%fig_scale = '-r384';
fig_scale = '-r192';
col = [0,0,1;1,0,0;0,.6,0]; %Plotting colors

filenames = {'samples_for_distance_analysis_16','samples_for_distance_analysis_17', 'samples_for_distance_analysis_18'};
for i=1:length(test_days) % three days, Wed thrus fri
    data_path = ['../Data/Data_for_Figures/', filenames{i}];
    S = load(data_path);


    % filter data
    if length(XWINDOW)>1  %!!!!! uncomment for filtering on testbed 
        ind = S.samples_xpos >= XWINDOW(1) & S.samples_xpos <= XWINDOW(2);
    end
    if length(TWINDOW)>1
        t_win_hr = floor(TWINDOW/100) ;
        t_win_min = mod(TWINDOW,100) ;
        t_win_s = ((t_win_hr-6)*60 + t_win_min)*60;
        ind = ind & S.samples_t >= t_win_s(1) & S.samples_t <= t_win_s(2);
    end
    samples_class_cell{i} = S.samples_class(ind);
    samples_dist_cell{i} = S.samples_dist(ind);
    samples_fcons_cell{i} = S.samples_fcons(ind);
    samples_fr_cell{i} = S.samples_fr(ind);
    samples_lane_cell{i} = S.samples_lane(ind);
    samples_speed_cell{i} = S.samples_speed(ind);
    samples_t_cell{i} = S.samples_t(ind);
    samples_xpos_cell{i} = S.samples_xpos(ind);
end  

%%%%%%%%%%%% preallocate statistics %%%%%%%%%%%%%%%
bin_counts   = zeros(length(test_days),length(bin_d));
binned_ind   = cell(1,length(bin_d));

for i=1:length(test_days)
    d=samples_dist_cell{i}; %used to put data into bins by distance
    for j = 1:length(bin_d) %bin the indices by distance from AV
        binned_ind{j}   = find(dist_to_a(j) <= d&d < dist_to_a(j+1));
        bin_counts(i,j)   = length(binned_ind{j});
    end
end

%% Plot data from bin counts.
axx = 350;
axy = [0, 1400000]; % [400000, 1400000]
equispaced_bins = [bin_d(1:34),-5,5,bin_d(38:end)];
fontsize = 18;
bin_counts_no_center = bin_counts(:,[1:35,37:end]);
wed_counts = bin_counts_no_center(1,:);
thurs_counts = bin_counts_no_center(2,:);
fri_counts = bin_counts_no_center(3,:);
%plot(bin_d, wed_counts)

figure(1)
subplot(3,1,1)
bar(equispaced_bins, wed_counts)
hold on
plot([0,0],axy,'k:','LineWidth',1.5)
plot([-30,-30],axy,'r:','LineWidth',2.5)
plot([30,30],axy,'r:','LineWidth',2.5)
title('Wednesday','FontSize',fontsize)
axis([-axx,axx,axy])
ax.FontSize = fontsize;
text(axx*.2,axy(1)+diff(axy)*.8,'behind AV','FontSize',fontsize,...
    'HorizontalAlignment','left','VerticalAlignment','bottom')
text(-axx*.2,axy(1)+diff(axy)*.8,'ahead of AV','FontSize',fontsize,...
    'HorizontalAlignment','right','VerticalAlignment','bottom')
hold off
set(gca,'fontsize',18)

figure(1)
subplot(3,1,2)
bar(equispaced_bins, thurs_counts)
hold on
plot([0,0],axy,'k:','LineWidth',1.5)
plot([-30,-30],axy,'r:','LineWidth',2.5)
plot([30,30],axy,'r:','LineWidth',2.5)
title('Thursday','FontSize',fontsize)
xlabel('Distance to nearest AV in same lane (m)','FontSize',fontsize)
ylabel('number of data points','FontSize',fontsize)
axis([-axx,axx,axy])

text(axx*.2,axy(1)+diff(axy)*.8,'behind AV','FontSize',fontsize,...
    'HorizontalAlignment','left','VerticalAlignment','bottom')
text(-axx*.2,axy(1)+diff(axy)*.8,'ahead of AV','FontSize',fontsize,...
    'HorizontalAlignment','right','VerticalAlignment','bottom')
hold off
set(gca,'fontsize',18)

figure(1)
subplot(3,1,3)
bar(equispaced_bins, fri_counts)
hold on
plot([0,0],axy,'k:','LineWidth',1.5)
plot([-30,-30],axy,'r:','LineWidth',2.5)
plot([30,30],axy,'r:','LineWidth',2.5)
title('Friday','FontSize',fontsize)
axis([-axx,axx,axy])
ax.FontSize = fontsize;
ax.YAxis.FontSize = 24;
text(axx*.2,axy(1)+diff(axy)*.8,'behind AV','FontSize',fontsize,...
    'HorizontalAlignment','left','VerticalAlignment','bottom')
text(-axx*.2,axy(1)+diff(axy)*.8,'ahead of AV','FontSize',fontsize,...
    'HorizontalAlignment','right','VerticalAlignment','bottom')
hold off
set(gca,'fontsize',18)