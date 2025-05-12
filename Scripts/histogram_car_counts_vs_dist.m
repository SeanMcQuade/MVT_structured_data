%this program counts the number of cars at dist from AV

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
