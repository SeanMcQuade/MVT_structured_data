% Specific for Nature
% (C) 2026/05/26 by Benjamin Seibold, added to by Sean McQuade
tic
%choose parameter: 16, 17, or 18 (16 = Wed, 17 = Thurs, or 18 = Fri). 
day = 16; 
%choose initial and terminal file, there are 24 files per day.
j_start = 1; j_end = 24;

%preallocate cell array to save rel speed strcat(datafolder,
filtered_dist_all_files = [];
filtered_speed_all_files = [];
all_av_dist = [];

%folder containing data
formatSpec = "../../results/slim/2022-11-%d/";
datafolder = sprintf(formatSpec, day);

for j=j_start:j_end
    %if 1 % activate upon first time, then deactivate
        data_files = dir(strcat(datafolder,'I-24MOTION_????-??-??_??-??-??.json'));
        % Load data file
        filename = data_files(j).name;
        fprintf('Loading %s ...',filename), tic
        fid = fopen(strcat(datafolder,filename));
        data = fread(fid,inf);
        fclose(fid);
        fprintf(' Done (%0.0fsec).\n',toc)
        fprintf('Decoding data structures ...'), tic
        data = jsondecode(char(data'));
        fprintf(' Done (%0.0fsec).\n',toc)

%histogram parameter
edges = -30:0.5:30;

% Process data
av_dist_all_secments = [];
rel_speed_all_segments = [];
for i = 1:length(data) % loop over all segments
    times = data(i).timestamp;
    av_id = data(i).downstream_engaged_av_id;
    av_dist = data(i).distance_to_downstream_engaged_av_meters; 
    
    ind = ~isnan(av_id); % only keep those instances where an AV is present
    if numel(ind)<2, continue, end % if segment too short, move to next segment
    times = times(ind); av_id = av_id(ind); av_dist = av_dist(ind);
    
    ind_av_switch = diff(av_id)~=0; % AV index switching
    ind_large_jump = abs(diff(av_dist)./diff(times))>70; % excessive jumps
    ind_switch = find(ind_av_switch|ind_large_jump); % indices of switching
    
    n_times = length(times);
    % finite difference stencil, usually centered, but one-sided at switchings
    stencil_forward = [2:n_times,n_times];
    stencil_forward(ind_switch) = stencil_forward(ind_switch)-1;
    stencil_backward = [1,1:n_times-1];
    stencil_backward(ind_switch+1) = stencil_backward(ind_switch+1)+1;
    
    rel_speed_raw = (av_dist(stencil_forward)-av_dist(stencil_backward))./...
        (times(stencil_forward)-times(stencil_backward)); % finite differences
    
    rel_speed = rel_speed_raw;
    for k = 1:40 % apply local averaging to get smooth relative speed
        rel_speed = (rel_speed(stencil_backward)+...
            rel_speed+rel_speed(stencil_forward))/3;
    end

    % save rel_speed and av_dist for all segments
    av_dist_all_secments = [av_dist_all_secments; av_dist];
    rel_speed_all_segments = [rel_speed_all_segments; rel_speed];
end

%only keep
clear filtered_ind filtered_dist filtered_speed;
lower_Bnd = 35;
upper_Bnd = 350;
filtered_ind_low = find(lower_Bnd < av_dist_all_secments);
filtered_ind_up = find(av_dist_all_secments < upper_Bnd);
filtered_ind = intersect(filtered_ind_low,filtered_ind_up);
filtered_dist = av_dist_all_secments(filtered_ind);
filtered_speed = rel_speed_all_segments(filtered_ind);

number_of_data = length(filtered_dist);
mean_rel_speed(j) = mean(filtered_speed,'omitnan');
median_rel_speed(j) = prctile(filtered_speed, 50);
stddev_rel_speed(j) = std(filtered_speed, 'omitnan');
first_quartile(j) = prctile(filtered_speed, 25);
third_quartile(j) = prctile(filtered_speed,75);
interquartile = [first_quartile(j), third_quartile(j)];

%save rel distance and speed from each iteration
filtered_dist_all_files = [filtered_dist_all_files; filtered_dist];
filtered_speed_all_files = [filtered_speed_all_files; filtered_speed];

%plot relative speed histogram for jth file
figure;
hold on
H = histogram(filtered_speed,edges);
plot(mean_rel_speed(j)*ones(2,1), [0,max(H.Values)],"LineWidth",2)
plot(median_rel_speed(j)*ones(2,1),[0,0.5*max(H.Values)],"LineWidth",2)
std_dev = [mean_rel_speed(j)-stddev_rel_speed(j),...
    mean_rel_speed(j) + stddev_rel_speed(j)];
plot(std_dev, [0.5,0.5],"LineWidth",3)
plot(interquartile, [0.5,0.5],"LineWidth",4)
fontsize = 24;
legend("Histogram", "Mean speed", "Median Speed", ...
                  "Standard Dev", "Interquartile","FontSize",fontsize)
xlabel("Relative Speed m/s","FontSize",fontsize);
ylabel_formatSpec = "Frequency of speeds between %d and %d m/s";
ylabel_string = sprintf(ylabel_formatSpec, lower_Bnd, upper_Bnd);
ylabel(ylabel_string,"FontSize",fontsize);


if day == 16 %Write the day in the title
    formatSpec = "Relative speed histogram Wed," + ...
        "for file %d, total data points = %d";
    title_string = sprintf(formatSpec,j,number_of_data);
elseif day == 17
    formatSpec = "Relative speed histogram Thurs, " + ...
        "for file %d, total data points = %d";
    title_string = sprintf(formatSpec,j,number_of_data);
elseif day == 18
    formatSpec = "Relative speed histogram Fri," + ...
        "for file %d, total data points = %d";
    title_string = sprintf(formatSpec,j,number_of_data);
end

title(title_string,"FontSize",fontsize);
formatSpecSave= "../../results/figures/Relative speed histogram, day %d for file j = %d.png";
savename = sprintf(formatSpecSave,day, j);
saveas(gcf,savename)
end

number_of_data_all_files = length(filtered_speed_all_files);
mean_all_fil_rel_speed = mean(filtered_speed_all_files,'omitnan');
median_fil_all_files = prctile(filtered_speed_all_files,50);
stddev_fil_all_rel_speed = std(filtered_speed_all_files, 'omitnan');
first_quartile_fil_all = prctile(filtered_speed_all_files, 25);
third_quartile_fil_all = prctile(filtered_speed_all_files,75);

%plot relative speed histogram for all files
figure;
hold on
H = histogram(filtered_speed_all_files,edges);
plot(mean_all_fil_rel_speed*ones(2,1), [0,max(H.Values)],"LineWidth",2)
plot(median_fil_all_files*ones(2,1),[0,0.5*max(H.Values)],"LineWidth",2)
std_dev = [mean_all_fil_rel_speed-stddev_fil_all_rel_speed,...
           mean_all_fil_rel_speed + stddev_fil_all_rel_speed];
interquartile_all = [first_quartile_fil_all, third_quartile_fil_all];
plot(std_dev, [10,10],"LineWidth",3)
plot(interquartile_all, [0.5,0.5],"LineWidth",4)
fontsize = 24;
legendformatSpec_mean = "mean speed = %3.3f";
mean_legend = sprintf(legendformatSpec_mean,mean_all_fil_rel_speed);
legendformatSpec_med = "median speed = %3.3f";
med_legend = sprintf(legendformatSpec_med,median_fil_all_files);
legend("histogram", mean_legend,med_legend,"Standard Dev", "Interquartile", ...
                                                       "FontSize",fontsize)

if day == 16 %Write the day in the title
    formatSpec = "Relative speed histogram Wed," + ...
        "for files %d to %d, total data points = %d";
    title_string = sprintf(formatSpec,j_start,j_end,number_of_data_all_files);
elseif day == 17
    formatSpec = "Relative speed histogram Thurs," + ...
        "for files %d to %d, total data points = %d";
    title_string = sprintf(formatSpec,j_start,j_end,number_of_data_all_files);
elseif day == 18
    formatSpec = "Relative speed histogram Fri," + ...
        "for file %d to %d, total data points = %d";
    title_string = sprintf(formatSpec, j_start, j_end,number_of_data_all_files);
end

title(title_string,"FontSize",fontsize);
xlabel("Relative Speed m/s","FontSize",fontsize);
ylabel_formatSpec = "Frequency of speeds between %d and %d m/s";
ylabel_string = sprintf(ylabel_formatSpec, lower_Bnd, upper_Bnd);
ylabel(ylabel_string,"FontSize",fontsize);
formatSpecSave= "../../results/figures/Relative speed histogram, day %d for files j = %d to %d.png";
savename = sprintf(formatSpecSave,day, j_start, j_end);
saveas(gcf,savename)

%run the script to show average speed in standard distance bins from AV
binned_relative_speed

% Returns h = 1 if the data is NOT Gaussian, h = 0 if it IS Gaussian
[h,p] = adtest(filtered_speed_all_files)
toc