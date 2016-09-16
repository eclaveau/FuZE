function [ppp, range, timecursor_ip, timecursor_m1m0, fontsize, htitle, h, h5, h10, h15, h20, h25, h30, h35, h40, h45, zz00, mm, mm5, mm10, mm15, mm20, mm25, mm30, mm35, mm40, mm45, mm00, hc] = main_animation(vis_start, vis_end, time, data)

ppp=figure('units','normalized','outerposition',[0 0 1 1]);

% Display Properties
ymargin = 1.05;
fontsize = 10;

% Title of whole figure
subplot(16,8,[1 8])
htitle = title('time');
set(htitle,'FontSize',fontsize);
axis off

% I_P plot
hip = subplot(16,8,[9 13 17 21 25 29]);
Ip = plot(time(vis_start:vis_end), data{1}(vis_start:vis_end),'LineWidth',2);
hold on
timecursor_ip = plot(NaN,NaN);
xlabel('Time [s]')
ylabel('I [A]')
hold off
pos = get(gca, 'Position');
pos(1) = pos(1);
pos(2) = ymargin*pos(2);
set(gca, 'Position', pos);
axis([-inf inf 0 inf])

% Normalized m1 to visualize quiescent period, from p5 to p45

m1_mean = data{106}(vis_start:vis_end);
m1_sig = data{107}(vis_start:vis_end);

subplot(16,8,[33 37 41 45 49 53])
plot(time(vis_start:vis_end)*1e6,m1_mean,'k','LineWidth',2);
hold on
plot(time(vis_start:vis_end)*1e6, m1_mean+m1_sig,'r');
plot(time(vis_start:vis_end)*1e6, 0.2*ones(length(time(vis_start:vis_end)),1),'k')
timecursor_m1m0 = plot(NaN,NaN);
xlabel('Time [\mus]')
ylabel('B_1/B_0')
pos = get(gca, 'Position');
pos(1) = pos(1);
pos(2) = ymargin*pos(2);
set(gca, 'Position', pos)
axis([-inf inf -inf 0.5])
ax = gca;
set(ax,'YTick',[0,0.2,0.4]);
hold off

% 3d plot of magnetic probe data
subplot(16,8,[14 16 22 24 30 32])% Top right bottom:[38 40 46 48 54 56]
h = plot3(NaN,NaN,NaN);
hold on
h5 = plot3(NaN,NaN,NaN);
h10 = plot3(NaN,NaN,NaN);
h15 = plot3(NaN,NaN,NaN);
h20 = plot3(NaN,NaN,NaN);
h25 = plot3(NaN,NaN,NaN);
h30 = plot3(NaN,NaN,NaN);
h35 = plot3(NaN,NaN,NaN);
h40 = plot3(NaN,NaN,NaN);
h45 = plot3(NaN,NaN,NaN);
zz00 = plot3(NaN,NaN,NaN);
axis([-max(data{56}) max(data{56}) 0 45 -max(data{56}) max(data{56})]);
xlabel('B [T]')
ylabel('z [cm]')
zlabel('B [T]')
daspect([1 20 1])
ax = gca;
set(ax,'YTick',[0,5,10,15,20,25,30,35,40,45]);
az = 130; el =5; % Point of view of the frame
view(az, el);
pos = get(gca, 'Position');
pos(1) = 1.*pos(1);
pos(2) = 1.08*pos(2);
set(gca, 'Position', pos)
grid on
hold off

% 3d plot of m1 and m0
subplot(16,8,[38 40 46 48 54 56])% Top right top:[14 16 22 24 30 32]
mm = plot3(NaN,NaN,NaN);
hold on
mm5 = plot3(NaN,NaN,NaN);
mm10 = plot3(NaN,NaN,NaN);
mm15 = plot3(NaN,NaN,NaN);
mm20 = plot3(NaN,NaN,NaN);
mm25 = plot3(NaN,NaN,NaN);
mm30 = plot3(NaN,NaN,NaN);
mm35 = plot3(NaN,NaN,NaN);
mm40 = plot3(NaN,NaN,NaN);
mm45 = plot3(NaN,NaN,NaN);
mm00 = plot3(NaN,NaN,NaN);
axis([-2.13/2*2.54 2.13/2*2.54 0 45 -2.13/2*2.54 2.13/2*2.54])
daspect([1 2 1])
xlabel('d [cm]')
ylabel('z [cm]')
zlabel('d [cm]')
% az = 50; el = 5; % 141202036
% az = 130; el = 15; % 141201018
view(az, el);
ax = gca;
set(ax,'YTick',[0,5,10,15,20,25,30,35,40,45]);
set(ax,'ZTick',[-2.7,-2,-1,0,1,2,2.7]);
set(ax,'XTick',[-2,0,2]);
set(1,'color','w');
pos = get(gca, 'Position');
pos(1) = pos(1);
pos(2) = 1.07*pos(2);
set(gca, 'Position', pos)
grid on
hold off

% Final figure parameters
ax=subplot(16,8,[57 128]);
[cont, hc] = contourf([], [], []);
axis([-4 49 -40 360])
xlabel('z [cm]')
ylabel(['\theta [' 176 ']'])
caxis([-1 1])
range = -1:.1:1;
ax.YTick = [0 45 90 135 180 225 270 315];
colorbar
colormap('jet')
grid on
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',fontsize)
end