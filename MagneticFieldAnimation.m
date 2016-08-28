close all;
clc; clear all;

% =========================================================================
% Part of the program you want to run
SaveLocal = 0; % Saves the workspace to a file (For use on a computer not connected to MDSplus)
LoadLocal = 1; % Load from existing file

% Starting and ending time parameters
tstart = 62.239 ; tend = 0; timesteps = 8;
% =========================================================================
shotnum = 151027024;
name = [num2str(shotnum) '_workspace'];
if LoadLocal == 0
    [time,data] = acquire(shotnum);
    if SaveLocal == 1
        save(name); % Save workspace to a local file
    end
else
    load(name)
end

Contour = 0;
M1y = 0;
m1m0 = 0; % Generate animation for m1m0
Mag = 1; % Generate animation for complete magnetic field
savevideo =0 ; % Save a video of the magnetic field lines
OnlyContour = 0; % Only plot contour on a separate figure
SingleFig =0; % Only plot and saves individual figures of the contour
% =========================================================================
% save video
% =========================================================================
if savevideo == 1
    video = VideoWriter(['Magnetic_Contour_' num2str(shotnum) '.avi']);
    open(video)
end
% ======================================================================
% Definition of the data
% ======================================================================

% Angle and axial position at which the probe are placed
theta = [0,pi/4,pi/2,3*pi/4,pi,5/4*pi,3/2*pi,7/4*pi,0];
zp0 = [0, 0, 0, 0, 0, 0, 0, 0 ,0]; % Position of the probe on the z axis
z0 = [0 5 10 15 20 25 30 35 40 45];
ppp=figure('units','normalized','outerposition',[0 0 1 1]); %,'visible','off'

% Display Properties
ymargin = 1.05;
fontsize = 10;
% Title (time stamp)
subplot(16,8,[1 8])
htitle = title('time');
set(htitle,'FontSize',fontsize);
axis off

% Time range of display
vis_start = find(data{1}>(max(data{1})/5),1,'first'); % Where to current first reach 10% of max current
vis_end = find(data{1}>(max(data{1})/5),1,'last'); % Where the current last reach 10% of max current

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

if OnlyContour ==1
    ax = gca;
    set(ax,'XTickLabel',[]);
    xlabel('')
    title('a)','FontSize',fontsize)
end

% Normalized m1 to visualize quiescent period, from p5 to p45
subplot(16,8,[33 37 41 45 49 53])
m1_mean = data{106}(vis_start:vis_end);
m1_sig = data{107}(vis_start:vis_end);
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

if OnlyContour ==1
    ax = gca;
    title('b)','FontSize',fontsize)
    pos(2) = .95*pos(2);
    set(gca, 'Position', pos)
end

if Mag == 1
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
    
    if OnlyContour ==1
        set(ax,'YTickLabel',[]);
        set(ax,'ZTickLabel',{-.2 0 .2});
        set(ax,'XTickLabel',[]);
        title('c)','FontSize',fontsize)
        pos(2) =.98*pos(2);
        set(gca, 'Position', pos)
    end
    hold off
end

if m1m0 == 1
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
    %     az = 130; el = 15; % 141201018
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
    if OnlyContour ==1
        set(ax,'YTickLabel',[]);
        set(ax,'ZTickLabel',{'',-2,'',0,'',2,''});
        set(ax,'XTickLabel',[]);
        %         xlabel('')
        %         ylabel('')
        %     zlabel('')
        pos(2) = .93*pos(2);
        set(gca, 'Position', pos)
        title('d)','FontSize',fontsize)
    end
end

if OnlyContour ==1
    onlyc= figure('units','normalized','outerposition',[0 0 1 1]);
    set(2,'color','w');
end

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

k = 0;
[t_pos,t_pos_start,t_pos_end] = startendtime(time,tstart,tend,timesteps);
tic

for t = t_pos
    k = k+1;
    fprintf('%2.2f %%\n',(time(t)-time(t_pos_start))/(time(t_pos_end)-time(t_pos_start))*100)
    % 8 array positions
    i = 2; % First position of the data in "node_string"
    b_p0 = [data{i}(t),data{i+1}(t),data{i+2}(t),data{i+3}(t),data{i+4}(t),data{i+5}(t),data{i+6}(t),data{i+7}(t),data{i}(t)]; % Value of the magnetic field at each angle with the last one repeating to have a complete circle
    [p0x, p0y] = pol2cart(theta, b_p0); % B value projected on the x and y axis
    i = 10; % Interpolation at data{i+5} because of non functionning b_p15_225_t probe
    b_p15 = [data{i}(t),data{i+1}(t),data{i+2}(t),data{i+3}(t),data{i+4}(t),0.5*data{i+4}(t)+0.5*data{i+6}(t),data{i+6}(t),data{i+7}(t),data{i}(t)];
    zp15 = 15*ones(1,9);
    [p15x, p15y] = pol2cart(theta, b_p15);
    i = 17;
    b_p30 = [data{i}(t),data{i+1}(t),data{i+2}(t),data{i+3}(t),data{i+4}(t),data{i+5}(t),data{i+6}(t),data{i+7}(t),data{i}(t)];
    zp30 = 30*ones(1,9);
    [p30x, p30y] = pol2cart(theta, b_p30);
    i = 25;
    b_p45 = [data{i}(t),data{i+1}(t),data{i+2}(t),data{i+3}(t),data{i+4}(t),data{i+5}(t),data{i+6}(t),data{i+7}(t),data{i}(t)];
    zp45 = 45*ones(1,9);
    [p45x, p45y] = pol2cart(theta, b_p45);
    
    % 4 array positions
    % The data point at the in-between angle is interpolated (average of the 2
    % adjacent data points)
    i = 33;
    b_p5 = [data{i}(t),   0.5*data{i}(t)+0.5*data{i+1}(t)   ,data{i+1}(t),   0.5*data{i+1}(t)+0.5*data{i+2}(t)   ,data{i+2}(t),   0.5*data{i+2}(t)+0.5*data{i+3}(t)   ,data{i+3}(t),   0.5*data{i+3}(t)+0.5*data{i}(t)    ,data{i}(t)];
    zp5 = 5*ones(1,9);
    [p5x, p5y] = pol2cart(theta, b_p5);
    i = 37; % p10_180 (i + 2) signal loss, so interpolating in the z direction.
    b_p10 = [data{i}(t),   0.5*data{i}(t)+0.5*data{i+1}(t)   ,data{i+1}(t),   0.5*data{i+1}(t)+0.5*(0.5*b_p5(5) + 0.5*b_p15(5)), 0.5*b_p5(5) + 0.5*b_p15(5),   0.5*(0.5*b_p5(5) + 0.5*b_p15(5))+0.5*data{i+2}(t)   ,data{i+2}(t),   0.5*data{i+2}(t)+0.5*data{i}(t)    ,data{i}(t)];
    zp10 = 10*ones(1,9);
    [p10x, p10y] = pol2cart(theta, b_p10);
    i = 40;
    b_p20 = [data{i}(t),   0.5*data{i}(t)+0.5*data{i+1}(t)   ,data{i+1}(t),   0.5*data{i+1}(t)+0.5*data{i+2}(t)   ,data{i+2}(t),   0.5*data{i+2}(t)+0.5*data{i+3}(t)   ,data{i+3}(t),   0.5*data{i+3}(t)+0.5*data{i}(t)    ,data{i}(t)];
    zp20 = 20*ones(1,9);
    [p20x, p20y] = pol2cart(theta, b_p20);
    i = 44;
    b_p25 = [data{i}(t),   0.5*data{i}(t)+0.5*data{i+1}(t)   ,data{i+1}(t),   0.5*data{i+1}(t)+0.5*data{i+2}(t)   ,data{i+2}(t),   0.5*data{i+2}(t)+0.5*data{i+3}(t)   ,data{i+3}(t),   0.5*data{i+3}(t)+0.5*data{i}(t)    ,data{i}(t)];
    zp25 = 25*ones(1,9);
    [p25x, p25y] = pol2cart(theta, b_p25);
    i = 48;
    b_p35 = [data{i}(t),   0.5*data{i}(t)+0.5*data{i+1}(t)   ,data{i+1}(t),   0.5*data{i+1}(t)+0.5*data{i+2}(t)   ,data{i+2}(t),   0.5*data{i+2}(t)+0.5*data{i+3}(t)   ,data{i+3}(t),   0.5*data{i+3}(t)+0.5*data{i}(t)    ,data{i}(t)];
    zp35 = 35*ones(1,9);
    [p35x, p35y] = pol2cart(theta, b_p35);
    i = 52;
    b_p40 = [data{i}(t),   0.5*data{i}(t)+0.5*data{i+1}(t)   ,data{i+1}(t),   0.5*data{i+1}(t)+0.5*data{i+2}(t)   ,data{i+2}(t),   0.5*data{i+2}(t)+0.5*data{i+3}(t)   ,data{i+3}(t),   0.5*data{i+3}(t)+0.5*data{i}(t)    ,data{i}(t)];
    zp40 = 40*ones(1,9);
    [p40x, p40y] = pol2cart(theta, b_p40);
    
    % Magnetic field at theta = 0
    b_th0 = [b_p5(1) b_p10(1) b_p15(1) b_p20(1) b_p25(1) b_p30(1) b_p35(1) b_p40(1) b_p45(1)]; % From 5 to 45
    % Magnetic field at theta = 180
    b_th180(:,k) = [b_p5(5) b_p10(5) b_p15(5) b_p20(5) b_p25(5) b_p30(5) b_p35(5) b_p40(5) b_p45(5)]; % From 5 to 45
    
    % m1 data from MDSplus
    i = 66;
    m1MDS = [data{i}(t),data{i+1}(t),data{i+2}(t),data{i+3}(t),data{i+4}(t),data{i+5}(t),data{i+6}(t),data{i+7}(t),data{i+8}(t), data{i+9}(t)];
    
    % Phase data for m1
    i = 76;
    phiMDS = [data{i}(t),data{i+1}(t),data{i+2}(t),data{i+3}(t),data{i+4}(t),data{i+5}(t),data{i+6}(t),data{i+7}(t),data{i+8}(t),data{i+9}(t)];
    
    % m0 data from MDSplus
    i = 56;
    m0MDS = [data{i}(t),data{i+1}(t),data{i+2}(t),data{i+3}(t),data{i+4}(t),data{i+5}(t),data{i+6}(t),data{i+7}(t),data{i+8}(t), data{i+9}(t)];
    
    % Defining the 'Z' values (In our case, the magnetic field at each theta
    % position) for the contour plot (matrix representing each
    % intersection of theta and axial position). The data doesn't use the last
    % data point of each magnetic field (b_pn) because the last data point is a
    % repetition of the first. The program was written this way to have
    % complete loop when ploted on a 3D plot.
    B = [b_p0(1:end-1)', b_p5(1:end-1)', b_p10(1:end-1)', b_p15(1:end-1)', b_p20(1:end-1)', b_p25(1:end-1)', b_p30(1:end-1)', b_p35(1:end-1)', b_p40(1:end-1)', b_p45(1:end-1)'];
    
    
    % ======================================================================
    % Time Stamp
    % ======================================================================
    % Time stamp in the title of the graphs
    fig_title = sprintf('%9.0f\n Time = %3.2f \\mus', shotnum, time(t)*1e6);
    set(htitle,'String',fig_title,'FontSize',fontsize);
    
    % Vertical timeline moving in the I_P and m1/m0 graphs
    tcx_ip = [time(t), time(t)];
    tcy_ip = [-10e4, 1.1*max(data{1})];
    set(timecursor_ip, 'XData',tcx_ip , 'YData',tcy_ip, 'Color', 'k', 'LineWidth', 1.5 )
    tcx_m1m0 = [time(t), time(t)]*1e6;
    tcy_m1m0 = [0, 0.5];
    set(timecursor_m1m0, 'XData',tcx_m1m0, 'YData',tcy_m1m0, 'Color', 'k', 'LineWidth', 1.5 )
    
    % ======================================================================
    % Plot - Magnetic data
    % ======================================================================
    if Mag == 1
        set(zz00, 'XData', zeros(1,10), 'YData', z0, 'ZData', zeros(1,10), 'Linewidth', 2, 'Color', 'k', 'LineStyle', '-.')
        %         set(h, 'XData', p0x, 'YData', zp0, 'ZData', p0y, 'Linewidth', 2.5, 'Color', 'b')
        set(h5, 'XData', p5x, 'YData', zp5, 'ZData', p5y, 'Linewidth', 2.5, 'Color', 'b')
        set(h10, 'XData', p10x, 'YData', zp10, 'ZData', p10y, 'Linewidth', 2.5, 'Color', 'b')
        set(h15, 'XData', p15x, 'YData', zp15, 'ZData', p15y, 'Linewidth', 2.5, 'Color', 'b')
        set(h20, 'XData', p20x, 'YData', zp20, 'ZData', p20y, 'Linewidth', 2.5, 'Color', 'b')
        set(h25, 'XData', p25x, 'YData', zp25, 'ZData', p25y, 'Linewidth', 2.5, 'Color', 'b')
        set(h30, 'XData', p30x, 'YData', zp30, 'ZData', p30y, 'Linewidth', 2.5, 'Color', 'b')
        set(h35, 'XData', p35x, 'YData', zp35, 'ZData', p35y, 'Linewidth', 2.5, 'Color', 'b')
        set(h40, 'XData', p40x, 'YData', zp40, 'ZData', p40y, 'Linewidth', 2.5, 'Color', 'b')
        set(h45, 'XData', p45x, 'YData', zp45, 'ZData', p45y, 'Linewidth', 2.5, 'Color', 'b')
    end
    
    if m1m0 == 1
        % ======================================================================
        % Plot - m0 and m1
        % ======================================================================
        
        ac = linspace(0,2*pi,100);
        zm(:,1) = 0*ones(1,length(ac));
        zm(:,2) = 5*ones(1,length(ac));
        zm(:,3) = 10*ones(1,length(ac));
        zm(:,4) = 15*ones(1,length(ac));
        zm(:,5) = 20*ones(1,length(ac));
        zm(:,6) = 25*ones(1,length(ac));
        zm(:,7) = 30*ones(1,length(ac));
        zm(:,8) = 35*ones(1,length(ac));
        zm(:,9) = 40*ones(1,length(ac));
        zm(:,10) = 45*ones(1,length(ac));
        
        m1y = m1MDS.*sin(phiMDS - pi/8);
        % Each column is a different Z position
        t_contour = linspace(0,2*pi,100);
        for jj = 1:length(m0MDS)
            r = 0.2*0.5*7.94*2.54/2; % Radius of circle representing m0
            
            % Displacement of column to be equal to real displacement.
            x0 = 0.5*(0.5*7.94*2.54)*min([m1MDS(jj)./max([m0MDS(jj), 0.01]),1])*cos(phiMDS(jj)-pi/8);
            y0 = 0.5*(0.5*7.94*2.54)*min([m1MDS(jj)./max([m0MDS(jj), 0.01]),1])*sin(phiMDS(jj)-pi/8);
            
            for ii = 1:length(t_contour)
                diff = min([abs(phiMDS(jj)-t_contour(ii)), abs(2*pi-abs(phiMDS(jj)-t_contour(ii)))]);
                if abs(diff) <=pi/2
                    Bcq(ii,jj) = min([m1MDS(jj)./max([m0MDS(jj), 0.01]),1]) * cos(diff);
                else
                    Bcq(ii,jj) = 0;
                end
            end
            
            % Determination of color depending on M1
            if m1MDS(jj) >= 0.2*m0MDS(jj)
                color(jj) = 'r';
            else
                color(jj) = 'k';
            end
            
            % Circle matrix where each column is a different axial (z) position.
            cx(:,jj) = r*cos(ac)+x0;
            cy(:,jj) = r*sin(ac)+y0;
            
            y(jj) = 0.5*(0.5*7.94*2.54)*min([m1y(jj)./max([m0MDS(jj), 0.01]),1]); %displacement = 0.5*(Inner diameter of tank)/2*(in2cm)*m1/m0
        end
        
        %         set(mm, 'XData', real(cx(:,1)), 'YData', zm(:,1), 'ZData', real(cy(:,1)), 'Linewidth', 1.5, 'Color', color(1), 'LineStyle', '-')
        set(mm5, 'XData', real(cx(:,2)), 'YData', zm(:,2), 'ZData', real(cy(:,2)), 'Linewidth', 1.5, 'Color', color(2), 'LineStyle', '-')
        set(mm10, 'XData', real(cx(:,3)), 'YData', zm(:,3), 'ZData', real(cy(:,3)), 'Linewidth', 1.5, 'Color', color(3), 'LineStyle', '-')
        set(mm15, 'XData', real(cx(:,4)), 'YData', zm(:,4), 'ZData', real(cy(:,4)), 'Linewidth', 1.5, 'Color', color(4), 'LineStyle', '-')
        set(mm20, 'XData', real(cx(:,5)), 'YData', zm(:,5), 'ZData', real(cy(:,5)), 'Linewidth', 1.5, 'Color', color(5), 'LineStyle', '-')
        set(mm25, 'XData', real(cx(:,6)), 'YData', zm(:,6), 'ZData', real(cy(:,6)), 'Linewidth', 1.5, 'Color', color(6), 'LineStyle', '-')
        set(mm30, 'XData', real(cx(:,7)), 'YData', zm(:,7), 'ZData', real(cy(:,7)), 'Linewidth', 1.5, 'Color', color(7), 'LineStyle', '-')
        set(mm35, 'XData', real(cx(:,8)), 'YData', zm(:,8), 'ZData', real(cy(:,8)), 'Linewidth', 1.5, 'Color', color(8), 'LineStyle', '-')
        set(mm40, 'XData', real(cx(:,9)), 'YData', zm(:,9), 'ZData', real(cy(:,9)), 'Linewidth', 1.5, 'Color', color(9), 'LineStyle', '-')
        set(mm45, 'XData', real(cx(:,10)), 'YData', zm(:,10), 'ZData', real(cy(:,10)), 'Linewidth', 1.5, 'Color', color(10), 'LineStyle', '-')
        set(mm00, 'XData', zeros(1,10), 'YData', z0, 'ZData', zeros(1,10), 'Linewidth', 2, 'Color', 'k', 'LineStyle', '-.')
        
        %         % GIF Generation
        %         frame=getframe(3);
        %         im=frame2im(frame);
        %         [imind,cmap]=rgb2ind(im,256);
        %         imwrite(imind,cmap,['gifs\',num2str(shotnum), '_m1m0.gif'],'gif','WriteMode','append','DelayTime',delay);
        %
    end
    % ======================================================================
    % Plot - Contour raw magnetic field
    % ======================================================================
    
    % Repeated angles so it is easier to find peaks.
    B = B(:,2:end);
    B_rep = B;
    B_rep(end+1:end+3,:) = B(1:3,:);
    theta_rep = theta;
    theta_rep(end:end+2) = 2*pi+theta(1:3);
    
    set(hc, 'XData', z0(2:end), 'YData', theta_rep/(2*pi)*360, 'ZData', B_rep/mean(B(:))-1,'LevelList', range, 'LineWidth', 2)
    drawnow
    if savevideo==1
        image_pert= frame2im(getframe(ppp));
        writeVideo(video,image_pert)
    end
    
    if SingleFig==1
        tlt = sprintf('%s_%2.0f_%4.0f_microsec.tif',num2str(shotnum),floor(time(t)*1e6),10000*(time(t)*1e6-floor(time(t)*1e6)));
        tlt_info = sprintf('Info_%s_%2.0f_%4.0f_microsec.tif',num2str(shotnum),floor(time(t)*1e6),10000*(time(t)*1e6-floor(time(t)*1e6)));
        if OnlyContour ==1
            imwrite(frame2im(getframe(onlyc)),tlt)
            imwrite(frame2im(getframe(ppp)),tlt_info)
        else
            imwrite(frame2im(getframe(ppp)),tlt)
        end
    end
    % =====================================================================
    % Fourier spectrum in the axial direction at each time
    % =====================================================================
    
    %     L = 0.4; % [m] From z = 5 cm to 45 cm
    %     n = length(b_th0(1:end-1));
    %     kf=(2*pi/L)*[0:(n/2-1) (-n/2):-1];
    %     kfs = fftshift(kf);
    %     fs = fft(b_th0(1:end-1)-mean(b_th0(1:end-1)));
    %
    %     figure(6)
    %     subplot(2,1,1)
    %     plot(z0(2:end), b_th0);
    %     title('B @ \theta = 0')
    %     subplot(2,1,2)
    %     plot(kfs, abs(fftshift(fs))')
    %     title('Frequency in space @ \theta = 0')
end
if savevideo==1
    close(video)
end

%%
% =========================================================================
% Principal Component Analysis with SVD only
% =========================================================================
[U,S,V] = svd(b_th180,'econ');

figure(8)
plot(diag(S)/sum(diag(S))*100,'o','LineWidth',3)
title('Principal Components [%]')

figure(9)
ii = 0;
for kk = [1 2 5 6]
    ii = ii +1;
    subplot(4,2,kk)
    plot(U(:,ii))
    ident = [num2str(ii) ' Spatial Mode'];
    title(ident)
    subplot(4,2,kk+2)
    plot(10^6*time(t_pos_start:t_pos_end),abs(V(:,ii)))
    ident = [num2str(ii) ' Temporal Mode'];
    title(ident)
end
% figure(10)
% for kk = [1 2 5 6]
%     ii = ii +1;
% subplot(4,2,kk)
% plot(U(:,ii))
% ident = [num2str(ii) ' Spatial Mode'];
% title(ident)
% subplot(4,2,kk+2)
% plot(10^6*time(t_pos_start:t_pos_end),abs(V(:,ii)))
% ident = [num2str(ii) ' Temporal Mode'];
% title(ident)
% end

% Representation of data in space time for one probe
[T,Z] = meshgrid(10^6*time(t_pos_start:t_pos_end),z0(2:end));
figure(7)
surf(T,Z,b_th180)
xlabel('t')
ylabel('x')
title('Probe data')

%    figure(7)
%    axis([5 45 0 0.40])
% for ii = 1:size(b_th180,2)
%  pause
%     plot(z0(2:end),b_th180(:,ii))
%     axis([5 45 0 0.40])
% end

% =========================================================================
% Dynamic Mode Decomposition Part (DMD)
% =========================================================================
clear u_modes u_dmd
dt = time(2) - time(1);

X1 = b_th180(:,1:end-1);
X2 = b_th180(:,2:end);

[U, Sigma, V] = svd(X1, 'econ');
U = U(:,1:3);
V = V(:,1:3);
Sigma = Sigma(1:3,1:3);
S = U'*X2*V*diag(1./diag(Sigma));
[eV,D] = eig(S);
mu = diag(D);
omega = log(mu)/dt;
Phi = U*eV;

u = b_th180(:,1);
y0 = Phi\u;

for iter=1:(size(b_th180,2))
    u_modes(:,iter) = (y0.*exp(omega*dt*(iter-1)));
    % X_inter(:,t) = diag(exp(omega*dt*(t-1)))*b;
end
u_dmd = Phi*u_modes;
% X_DMD = Phi*X_inter;

figure(11)
surf(T(:,1:end),Z(:,1:end),real(u_dmd))
xlabel('t')
ylabel('x')
title('DMD generated data')
figure(12)
surf(T(:,1:end),Z(:,1:end),real(u_dmd)-b_th180(:,1:end))
xlabel('t')
ylabel('x')
title('Error')

% =========================================================================
% Modes visualization
% =========================================================================
figure(13)
plot(real(omega),'ok','LineWidth',3)
title('Real part of omega')

figure(14)
for kk = 1:3
    clear Phim omegam y0 u_dmd u_modesm
    Phim = Phi(:,kk);
    y0 = Phim\u;
    omegam = omega(kk);
    
    for iter=1:(size(b_th180,2)-1)
        u_modesm(:,iter) = (y0.*exp(omegam*dt*(iter-1)));
        % X_inter(:,t) = diag(exp(omega*dt*(t-1)))*b;
    end
    u_dmd = Phim*u_modesm;
    mode{kk} = u_dmd;
    subplot(2,2,kk)
    surf(T(:,1:end-1),Z(:,1:end-1),real(u_dmd))
    tit = ['Mode ' num2str(kk)];
    title(tit)
end

% Frequencies of oscillations

f = imag(omega(2))/(2*pi); % [Hz]




