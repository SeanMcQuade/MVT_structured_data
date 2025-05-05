% script to plot AV fuel consumption for for 16 17 and 18 with additions16
% from Benni Seibold, modified by Sean McQuade to loop through different
% paramters, and save each plot. check for largest XWINDOW: x_in =
% min(x^I24_in,x^GPS_in-350m) or +350m, depending on direction of x-coord

tic

MAX_DIST = 350; % (m) maximum relative distance to an ego vehicle to be considered 
XWINDOW = [MAX_DIST-122 7000]; % (m) , set to -1 to include the full testbed. Find max and min (samples.xpos) and subtract 400m

%TWINDOWS = [0645, 0745; 0700, 0800; 0715, 0815; 0730, 0830; 0645, 0915; 0645, 0930];
% 6:45 to 9:15 is a multiple of 15min that shows mean speed <20 m/s for all days with AV present
TWINDOWS = [0645, 0915];%; 0700, 0730; 0730, 0800]; %windows of intrest 12/5/2024

for i_t = 1:size(TWINDOWS,1)
TWINDOW = TWINDOWS(i_t,:);
min_spd_mean_only = 1; % (m/s) % minimum speed for samples (to remove ideling vehicles)
min_spd_raw = 0;
dist_to_a = round(linspace(-MAX_DIST,MAX_DIST,MAX_DIST*2/10+1));
dist_to_a = [dist_to_a(1:floor(length(dist_to_a)/2)),-1,1,dist_to_a(floor(length(dist_to_a)/2)+2:end)];
bin_d = (dist_to_a(1:end-1)+dist_to_a(2:end))/2;
flag_save = 1;
fig_res = [1600 1000];
%fig_scale = '-r384';
fig_scale = '-r192';
col = [0,0,1;1,0,0;0,.6,0]; %Plotting colors

% binindices 
% ind_unmasked_left = find(dist_to_a<-40);
% ind_unmasked_right = find(dist_to_a>40);

masking_dist = 30; % bins closer than this num meters are masked
ind_unmasked_left = find(dist_to_a(2:end)<=-masking_dist);
ind_masked_left = [ind_unmasked_left(end) find(-masking_dist<dist_to_a(2:end) & dist_to_a(2:end)<0)];
ind_unmasked_right = find( masking_dist <= dist_to_a(1:end-1)); 
ind_masked_right = [find(0<dist_to_a(1:end-1) & dist_to_a(1:end-1) < masking_dist) ind_unmasked_right(1)];
indices_partitions{1} = ind_unmasked_left;
indices_partitions{2} = ind_masked_left;
indices_partitions{3} = ind_masked_right;
indices_partitions{4} = ind_unmasked_right;

%% load semi-filtered raw data      Uncomment line 26 to line 76 to create data files, 
%   uncomment 79 to load data.
test_days = [16,17,18];
filenames = {'samples_for_distance_analysis_16','samples_for_distance_analysis_17', 'samples_for_distance_analysis_18'};
for i=1:length(test_days) % three days, Wed thrus fri
    S = load(filenames{i});

    % filter data
    if length(XWINDOW)>1
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
mean_fcons   = zeros(length(test_days),length(bin_d));
binned_fcons = cell(length(test_days),length(bin_d));
binned_fr    = cell(length(test_days),length(bin_d));
binned_speed = cell(length(test_days),length(bin_d));
effective_fc = zeros(length(test_days),length(bin_d));
median_fcons = zeros(length(test_days),length(bin_d));
bin_counts   = zeros(length(test_days),length(bin_d));
binned_ind   = cell(1,length(bin_d));

for i=1:length(test_days)
    d=samples_dist_cell{i}; %used to put data into bins by distance
    for j = 1:length(bin_d) %bin the indices by distance from AV
        binned_ind{j}   = find(dist_to_a(j) <= d&d < dist_to_a(j+1));
        binned_speed{i,j} = samples_speed_cell{i}(binned_ind{j});
        binned_fr{i,j}    = samples_fr_cell{i}(binned_ind{j});
        binned_fcons{i,j} = samples_fcons_cell{i}(binned_ind{j});
        effective_fc(i,j) = sum(binned_fr{i,j})/max(sum(binned_speed{i,j}),1e-6);
        mean_fcons(i,j)    = mean(binned_fcons{i,j}(binned_speed{i,j}>=min_spd_mean_only));
        median_fcons(i,j) = prctile(binned_fcons{i,j},50);
        bin_counts(i,j)   = length(binned_ind{j});
    end
end

%% Create triplets of mean med eff, partition them into Left, Masked_left, Masked_right, Right; Calculate regression
stats_to_plot         = cell(length(test_days),1);
partitioned_stats     = cell(length(test_days),4);
d_left                = zeros(length(test_days),length(ind_unmasked_left));
d_right               = zeros(length(test_days),length(ind_unmasked_right));
d_masked_left         = zeros(length(test_days),length(ind_masked_left)); 
d_masked_right        = zeros(length(test_days),length(ind_masked_right));
val_left              = zeros(length(test_days),length(ind_unmasked_left));
val_right             = zeros(length(test_days),length(ind_unmasked_right));
val_masked_left       = zeros(length(test_days),length(ind_masked_left)); 
val_masked_right      = zeros(length(test_days),length(ind_masked_right)); 
W_reg_left            = cell(length(test_days),1);
W_reg_right           = cell(length(test_days),1);
S_regression          = cell(length(test_days),1);

for i=1:length(test_days)
    stats_to_plot{i} = [mean_fcons(i,:); effective_fc(i,:); median_fcons(i,:)];
    % partitioned_stats{1}left, {2}masked left, {3}masked right, {4}right. Each has three rows [mean_fcons; effective_fc; median_fcons]
    [partitioned_single_day] = mask_nearby_points(stats_to_plot{i},indices_partitions);
    [W_reg_left{i}, W_reg_right{i}] = find_weighted_regression_two_sided(partitioned_single_day, bin_counts(i,:), indices_partitions, bin_d);
    [S_regression{i}] = find_simple_regression_one_sided(partitioned_single_day, indices_partitions, bin_d);
    [partitioned_stats{i,1},partitioned_stats{i,2},partitioned_stats{i,3},partitioned_stats{i,4}] = partitioned_single_day{1:4};
end

%% Calculate Mean Bulk speed and AVs for Wed. Thurs. and Fri.
% %% Load Variables
load data_for_av_activation_wed_thurs_fri(2024).mat

t_cell_array = {t_16 t_17 t_18};
x_cell_array = {x_16 x_17 x_18};
data_cell_array = {data_16 data_17 data_18};
Str_array = ["fields_motion16" "fields_motion17" "fields_motion18"];
Fig_string = ["bulkSpeedAVs1" "bulkSpeedAVs2" "bulkSpeedAVs3"];
for i=1:3 %wed. thurs. fri.
    load(Str_array(i))
    data = data_cell_array{i};
    t = t_cell_array{i};
    x = x_cell_array{i};
for ii = 1:length(t)
    
    t_temp = t(ii);
    u_temp_v = field.value.U(ii,:);
    rho_temp_v = field.value.Rho(ii,:)*50;
    rho_temp_v = rho_temp_v(~isnan(u_temp_v));
    u_temp_v = u_temp_v(~isnan(u_temp_v));
    
    total_rho_temp = sum(rho_temp_v);
    u_temp_mean(ii) = sum(u_temp_v .* rho_temp_v)/total_rho_temp;
    u_avg(ii) = mean(u_temp_v);
    
    AVs_on_the_road = data([data.direction]==-1 & [data.first_timestamp]<=  t_temp & [data.last_timestamp]>=  t_temp + 5);
    nAVs_active_temp = 0;
    nAVs_temp = 0;
    for ii2 = 1:length(AVs_on_the_road)
        av = AVs_on_the_road(ii2);
        x_min = -122;
        x_max = 7000;
        av_in_cell_ind = av.timestamp >= t_temp & av.timestamp <= t_temp +5 & av.x_position >= x_min & av.x_position <=   x_max; %Sulaiman suggestion to count only AVs in xwindow.
        %av_in_cell_ind  = av.timestamp >= t_temp & av.timestamp <= t_temp +5 ;
        if isfield(av,"control_car") && ~isempty(av.control_car)
                is_av_control = av.control_car(av_in_cell_ind);
        else
                is_av_control = av.controller_engaged(av_in_cell_ind);
        end
        is_av_control = (sum(is_av_control)/length(is_av_control))>=0.5;  % the AV is considered active if it has an active signal for more than half the time bin of 5s
        if is_av_control
                nAVs_active_temp = nAVs_active_temp + 1;
        end
        nAVs_temp = nAVs_temp +1;
    end
    nAVs_active(ii) = nAVs_active_temp;
    nAVs(ii) = nAVs_temp;
end
%% Plotting bulk speed 
lightblue = 1/2*[0,0,1]+1/2*[1,1,1];
lightred = 1/2*[1,0,0]+1/2*[1,1,1];
fontsize = 14;
TWINDOWstr = num2str(TWINDOW); Firsthour=TWINDOWstr(1); Firstminutes=TWINDOWstr(2:3); Secondhour = TWINDOWstr(6); Secondminutes = TWINDOWstr(7:8);
TT = datetime(t, 'convertfrom', 'posixtime', 'Format', 'HH:mm:ss.SSS','TimeZone' ,'America/Chicago') ; %[7:00:00,9:00:00;7:20:00,8:50:00;7:00:00,8:30:00]
plot_start_t = datetime(2022,11,test_days(i), str2num(Firsthour), str2num(Firstminutes),0, 'TimeZone' ,'America/Chicago');
plot_end_t = datetime(2022,11,test_days(i), str2num(Secondhour), str2num(Secondminutes),0, 'TimeZone' ,'America/Chicago');
TT_before_ind = TT <= plot_start_t; TT_after_ind = TT <= plot_end_t; 
TT_inside_ind = TT >= plot_start_t & TT <= plot_end_t;

figure(1)
subplot(3,3,3*i-2)
hold on
yyaxis right
hold on %plot number of AVs
plot(TT(TT <= plot_start_t), nAVs_active(TT <= plot_start_t),'Color',lightred,'LineStyle','-','LineWidth',1.5)
plot(TT(TT >= plot_end_t), nAVs_active(TT >= plot_end_t),'Color',lightred,'LineStyle','-','LineWidth',1.5)
plot(TT(TT >= plot_start_t & TT <= plot_end_t), nAVs_active(TT >= plot_start_t & TT <= plot_end_t),'r-','LineWidth',1.5)
plot([plot_start_t plot_start_t],[0,55],'k:','LineWidth',1.5) %vertical lines
plot([plot_end_t plot_end_t],[0,55],'k:','LineWidth',1.5)
ylabel('Number of active AVs')
set(gca,'YColor','r');
grid on
axx=[plot_start_t-minutes(30),plot_end_t+minutes(30)];
axy = [0,55];
title('Traffic state and AV deployment vs. time')
text(plot_start_t+minutes(2),(axy(2) - axy(1))*0.9,'start time','FontSize',fontsize)
text(plot_end_t-minutes(47),(axy(2) - axy(1))*0.9,'end time','FontSize',fontsize)
set(gca,'FontSize',fontsize)
xlim(axx)
ylim(axy)

yyaxis left
plot(TT(TT <= plot_start_t), u_temp_mean(TT <= plot_start_t),'Color',lightblue,'LineStyle','-','LineWidth',2)
plot(TT(TT >= plot_end_t), u_temp_mean(TT >= plot_end_t),'Color',lightblue,'LineStyle','-','LineWidth',2)
plot(TT(TT >= plot_start_t & TT <= plot_end_t), u_temp_mean(TT >= plot_start_t & TT <= plot_end_t),'Color','b','LineStyle','-','LineWidth',2)
ylabel('Mean speed (m/s)')
ylim([0,32])
mean_bulk_speed_in_window(i) = mean(u_temp_mean(TT >= plot_start_t & TT <= plot_end_t));
% plot([plot_start_t, plot_end_t],[mean_bulk_speed_in_window(i), mean_bulk_speed_in_window(i)], 'b--','LineWidth',1.5); %plot mean bulk speed
hold off

end %for loop of days (Wed Thurs Fri)

%% Plot in 3x3 grid            
figure(1)
% for i = 1:length(test_days)
time_string = [TWINDOWstr(1), ':', TWINDOWstr(2:3), ' to ', TWINDOWstr(6), ':', TWINDOWstr(7:8)];
datestr_ary = {['Wed Nov 16, time ' time_string],['Thu Nov 17, time ' time_string] ,['Fri Nov 18, time ' time_string]};

for i=1:length(test_days)
    hold on
    datestr = datestr_ary{i};
    [partitioned_stats_per_day{1},partitioned_stats_per_day{2},partitioned_stats_per_day{3},partitioned_stats_per_day{4}]=partitioned_stats{i,:};
    subplot(3,3,3*i-1)
    plot_two_sided(datestr,col,partitioned_stats_per_day,W_reg_left{i}, W_reg_right{i}, bin_d, indices_partitions)
    subplot(3,3,3*i)
    plot_one_sided(datestr,col,partitioned_stats_per_day,S_regression{i}, bin_d, indices_partitions)
end


% Adjust subfig placements and size
posx = [.05,.355,.76];
posy = [.74,0,.07]; posy(2) = (posy(1)+posy(3))/2;
poss = [.21,.22]; posw = [.34, poss(2)];
subplot(3,3,3), set(gca,'Position',[posx(3),posy(1),poss])
subplot(3,3,6), set(gca,'Position',[posx(3),posy(2),poss])
subplot(3,3,9), set(gca,'Position',[posx(3),posy(3),poss])

subplot(3,3,1), set(gca,'Position',[posx(1),posy(1),poss])
subplot(3,3,4), set(gca,'Position',[posx(1),posy(2),poss])
subplot(3,3,7), set(gca,'Position',[posx(1),posy(3),poss])

subplot(3,3,2), set(gca,'Position',[posx(2),posy(1),posw])
subplot(3,3,5), set(gca,'Position',[posx(2),posy(2),posw])
subplot(3,3,8), set(gca,'Position',[posx(2),posy(3),posw])

if flag_save
    fname = ['fig_nature_fuel_results_',TWINDOWstr(1),TWINDOWstr(2),TWINDOWstr(3),'_',TWINDOWstr(6),TWINDOWstr(7),TWINDOWstr(8)];
    set(gcf,'Position',[10 40 fig_res],'PaperPositionMode','auto')
    print(fname,'-dpng',fig_scale);
    savefig(fname)
end
end
toc

%========================================================================
function plot_two_sided(datestr,col,partitioned_stats,W_reg_left, W_reg_right, bin_d, indices_partitions)

axx = 350;
axy = [0.06,0.13];
fontsize = 14;
colb = 1/5*col+4/5*[1,1,1]; % brightest version of color
colsb = 1/2*col+1/2*[1,1,1]; %slightly brighter version
cold = 1/2*col+1/2*[0,0,0]; % darker version of color
% color_ary = {colsb,col,cold};

AV_x1 = bin_d(indices_partitions{1}); AV_y1 = partitioned_stats{1};              %left_unmasked
AV_x2 = bin_d(:,indices_partitions{2}); AV_y2 = partitioned_stats{2};           %left_masked
AV_x3 = bin_d(:,indices_partitions{3}); AV_y3 = partitioned_stats{3};           %right_masked
AV_x4 = bin_d(:,indices_partitions{4}); AV_y4 = partitioned_stats{4};          %right_unmasked

plot([0,0],axy,'k:','LineWidth',1.5)
hold on
grid on
text(axx*.03,axy(1)+diff(axy)*.90,'behind AV','FontSize',fontsize,...
    'HorizontalAlignment','left','VerticalAlignment','bottom')
text(-axx*.03,axy(1)+diff(axy)*.90,'ahead of AV','FontSize',fontsize,...
    'HorizontalAlignment','right','VerticalAlignment','bottom')
stats_ary = {'mean FC','effective FC','median FC'};
linestyle_ary = {'-','-','-'};
linesize_ary = [2.5, 2.5, 2.5];

ind_quantity_plotted = [1,2,3]; %mean effective median

% for i = 1:length(stats_ary)% 
for i = ind_quantity_plotted
    plot(AV_x1,W_reg_left(i,:),'--','Color',cold(i,:),'LineWidth',2.5);  %plot regression
    plot(AV_x4,W_reg_right(i,:),'--','Color',cold(i,:),'LineWidth',2.5);
    plot(AV_x1,AV_y1(i,:) ,linestyle_ary{i} ,'Color',col(i,:),'LineWidth',linesize_ary(i)); 
    plot(AV_x2,AV_y2(i,:) ,linestyle_ary{i} ,'Color',colb(i,:),'LineWidth',2); 
    plot(AV_x3,AV_y3(i,:) ,linestyle_ary{i} ,'Color',colb(i,:),'LineWidth',2); 
    hl(i) = plot(AV_x4,AV_y4(i,:) ,linestyle_ary{i} ,'Color',col(i,:),'LineWidth',linesize_ary(i)); 
    title(['Fuel Consumption ',datestr],'FontSize',fontsize)
end
axis([-axx,axx,axy])
grid on
set(gca,'xtick',-300:100:300)
legend(hl,stats_ary{ind_quantity_plotted}) %%%check
title(['Fuel consumption (FC). ',datestr],'FontSize',fontsize)
xlabel('Distance to nearest AV in same lane (m)','FontSize',fontsize)
ylabel('Fuel consumption (g/m)','FontSize',fontsize)
set(gca,'FontSize',fontsize)
end
                      
function plot_one_sided(datestr,col,partitioned_stats,S_regression, bin_d, indices_partitions)
axx = 350;
fontsize = 14;
colb = 1/4*col+3/4*[1,1,1]; % brightest version of color
colsb = 2/3*col+1/3*[1,1,1]; %slightly brighter version
cold = 1/2*col+1/2*[0,0,0]; % darker version of color
% color_ary = {colsb,col,cold};
linestyle_ary = {'-','-','-'};
linesize_ary = [2.5, 2.5, 2.5];

AV_x1 = bin_d(indices_partitions{1}); AV_y1 = partitioned_stats{1};              %left_unmasked
AV_x2 = bin_d(:,indices_partitions{2}); AV_y2 = partitioned_stats{2};           %left_masked
AV_x3 = bin_d(:,indices_partitions{3}); AV_y3 = partitioned_stats{3};           %right_masked
AV_x4 = bin_d(:,indices_partitions{4}); AV_y4 = partitioned_stats{4};          %right_unmasked
plot([0,axx],[0,0],'k:','LineWidth',1.5)
hold on
grid on
stats_ary = {'mean FC','effective FC','median FC'};
for i = 1:length(stats_ary)
    hl(i) = plot(AV_x4,-(AV_y4(i,:)./AV_y1(i,end:-1:1)-1)*100,linestyle_ary{i},'Color',col(i,:),'LineWidth',linesize_ary(i));
    plot(AV_x3,-(AV_y3(i,:)./AV_y2(i,end:-1:1)-1)*100,linestyle_ary{i},'Color',colb(i,:),'LineWidth',2)
    plot(AV_x4,S_regression(i,:),'--','Color',cold(i,:),'LineWidth',2.5);
end
axis([0,axx,-2,30])
axis ij
grid on
set(gca,'xtick',0:100:300)
ytickformat('percentage')
legend(hl,stats_ary,'Location','southeast')
title('Fuel reduction ahead vs. behind AV','FontSize',fontsize) %changed from "behind vs ahead" 12/10/24
xlabel('Distance to nearest AV in same lane (m)','FontSize',fontsize)
ylabel('Fuel consumption reduction','FontSize',fontsize)
set(gca,'FontSize',fontsize)
end

function [partitioned_y] = mask_nearby_points(y,indices_partitions)
partitioned_y{1,1} = y(:,indices_partitions{1});
partitioned_y{1,2} = y(:,indices_partitions{2});
partitioned_y{1,3} = y(:,indices_partitions{3});
partitioned_y{1,4} = y(:,indices_partitions{4});
end

function [W_reg_left, W_reg_right] = find_weighted_regression_two_sided(partitioned_stats, Bin_Count_day, indices_partitions, bin_d)
x_vals(1,:) = transpose(bin_d(indices_partitions{1}));                     %left x_values
x_vals(2,:) = transpose(bin_d(indices_partitions{4}));                     %right x_values
weights{1} = diag(Bin_Count_day(indices_partitions{1}));                   %left weights
weights{2} = diag(Bin_Count_day(indices_partitions{4}));                   %right weights

unmasked_partitions = size(indices_partitions);                            % this chooses the first and last as "unmasked"
for i = 1:2 %four partitions LEFT, masked, masked, RIGHT
    y_vals = partitioned_stats{unmasked_partitions(i)};     
    A_matrix = [transpose(x_vals(i,:)) ones(size(transpose(indices_partitions{1})))];
    for j = 1:3
        Weighted_Solution = (transpose(A_matrix)*weights{i}*A_matrix)\(transpose(A_matrix)*weights{i}*transpose(y_vals(j,:)));
        W_regression(j,:) = A_matrix*Weighted_Solution;
    end
    partitioned_W_regression{i} = W_regression;
end
W_reg_left = partitioned_W_regression{1}; W_reg_right = partitioned_W_regression{2};
end

function [S_regression] = find_simple_regression_one_sided(partitioned_stats, indices_partitions, bin_d)
    x_vals = bin_d(indices_partitions{4});                                 %right x_values
    y_vals_left  = transpose(partitioned_stats{1}); 
    y_vals_right = transpose(partitioned_stats{4}); 
    ratio_behind_vs_ahead =  -(y_vals_right./y_vals_left(end:-1:1,:)-1)*100;
    % -(AV_y3(i,:)./AV_y2(i,end:-1:1)-1)
    A_matrix = [transpose(x_vals) ones(size(transpose(indices_partitions{1})))];
    Solution_ratio = (transpose(A_matrix)*A_matrix)\(transpose(A_matrix)*ratio_behind_vs_ahead);
    S_regression = transpose(A_matrix*Solution_ratio);
end

% function [W_reg_left, W_reg_right] = find_regression(partitioned_stats, Bin_Count_day, indices_partitions, bin_d)
% x_vals(1,:) = transpose(bin_d(indices_partitions{1}));                     %left x_values
% x_vals(2,:) = transpose(bin_d(indices_partitions{4}));                     %right x_values
% if length(Bin_Count_day) > 1
%     weights{1} = diag(Bin_Count_day(indices_partitions{1}));                   %left weights
%     weights{2} = diag(Bin_Count_day(indices_partitions{4}));                   %right weights
% 
%     unmasked_partitions = size(indices_partitions);                            % this chooses the first and last as "unmasked"
%     for i = 1:2 %four partitions LEFT, masked, masked, RIGHT
%         y_vals = partitioned_stats{unmasked_partitions(i)};
%         A_matrix = [transpose(x_vals(i,:)) ones(size(transpose(indices_partitions{1})))];
%         for j = 1:3
%             Weighted_Solution = (transpose(A_matrix)*weights{i}*A_matrix)\(transpose(A_matrix)*weights{i}*transpose(y_vals(j,:)));
%             W_regression(j,:) = A_matrix*Weighted_Solution;
%         end
%         partitioned_W_regression{i} = W_regression;
%     end
%     W_reg_left = partitioned_W_regression{1}; W_reg_right = partitioned_W_regression{2};
% else
% end
% end
