close all;
clc; clear all;

% =========================================================================
% Part of the program you want to run
SaveLocal = 0; % Saves the workspace to a file (For use on a computer not connected to MDSplus)
LoadLocal = 1; % Load from existing file
% =========================================================================
shotnum = 151027024;
name = [num2str(shotnum) '_workspace'];
if LoadLocal == 0
    [time,data,node_string] = acquire(shotnum);
else
    load(name)
end

if SaveLocal == 1
    save(name); % Save workspace to a local file
end
% Starting and ending time parameters  62.239
tstart = 62.239; tend = 0; timesteps = 8;
Contour = 0;
M1y = 0;
m1m0 = 0; % Generate animation for m1m0
Mag = 1; % Generate animation for complete magnetic field
savevideo =0 ; % Save a video of the magnetic field lines
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

% Time range of display
vis_start = find(data{1}>(max(data{1})/5),1,'first'); % Where to current first reach 10% of max current
vis_end = find(data{1}>(max(data{1})/5),1,'last'); % Where the current last reach 10% of max current

% Angle and axial position at which the probes are placed
theta = [0,pi/4,pi/2,3*pi/4,pi,5/4*pi,3/2*pi,7/4*pi,0];
zp0 = [0, 0, 0, 0, 0, 0, 0, 0 ,0]; % Position of the probes on the z axis
z0 = [0 5 10 15 20 25 30 35 40 45];

% Figure initialization parameters
[ppp, range, timecursor_ip, timecursor_m1m0, fontsize, htitle, h, h5, h10,...
    h15, h20, h25, h30, h35, h40, h45, zz00, mm, mm5, mm10, mm15, mm20, ...
    mm25, mm30, mm35, mm40, mm45, mm00, hc] = main_animation(vis_start, vis_end, time, data);

k = 0;
% Animation parameters. Start and end time.
[t_pos,t_pos_start,t_pos_end] = startendtime(time,tstart,tend,timesteps);
tic

for t = t_pos
    k = k+1;
    
    % Progress tracking
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
% [U,S,V] = svd(b_th180,'econ');

% % Representation of data in space time for one probe
[T,Z] = meshgrid(10^6*time(t_pos_start:t_pos_end),z0(2:end));
% figure(7)
% surf(T,Z,b_th180)
% xlabel('t')
% ylabel('x')
% title('Probe data')

% =========================================================================
% Dynamic Mode Decomposition Part (DMD)
% =========================================================================
% dt = time(2)-time(1);

rank = 8;
b_untrans = b_th180;
b_th180_T = b_th180';

% Three modes svd reconstruction
% [Ut, Sigmat, Vt] = svd(b_th180_T, 'econ');
% figure(7)
% subplot(1,2,1)
% surf(Ut(:,1:3)*Sigmat(1:3,1:3)*Vt(:,1:3)')
% title('USV reconstruction with 3 first modes')
% subplot(1,2,2)
% surf(Ut(:,1:3)*Sigmat(1:3,1:3)*Vt(:,1:3)'-b_th180_T)
% title('Error between data and 3 modes reconstructions')
% err_svd = mean(mean(abs(Ut(:,1:3)*Sigmat(1:3,1:3)*Vt(:,1:3)'-b_th180_T)));
% 
% figure(8)
% plot(diag(Sigmat)/sum(diag(Sigmat))*100,'.','LineWidth',3)
% title('Principal Components [%]')
% 
% figure(9)
% ii = 0;
% for kk = [1 2 5 6]
%     ii = ii +1;
%     subplot(4,2,kk)
%     plot(Ut(:,ii))
%     ident = [num2str(ii) ' Spatial Mode'];
%     title(ident)
%     subplot(4,2,kk+2)
%     plot(10^6*time(t_pos_start:t_pos_end),abs(Vt(:,ii)))
%     ident = [num2str(ii) ' Temporal Mode'];
%     title(ident)
% end

% X1 = b_th180_T(:,1:end-1);
% X2 = b_th180_T(:,2:end);
% u = b_th180_T(:,1);
% dt = 1;
% t = 0:9;
% 
% [U1, Sigma1, V1] = svd(X1, 'econ');
% U = U1(:,1:rank);
% V = V1(:,1:rank);
% Sigma = Sigma1(1:rank,1:rank);
% S = U'*X2*V*diag(1./diag(Sigma));
% [eV,D] = eig(S);
% mu = diag(D);
% omega = log(mu)/dt;
% Phi = U*eV;
% y0 = Phi\u;
% 
% for iter=1:(size(b_th180_T,2))-1
%     u_modes(:,iter) = (y0.*exp(omega*t(iter)));
% 
% end
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  data = b_th180_T;
  r = rank;
  t = linspace(0,9,9);
  fs = 0;
  
  x = data(:, 1:end-1);
  y = data(:, 2:end);
  first_frame = data(:, 1);

  % finalize our timebase
  dt = mean(diff(t));
  tfinal = t(end) + dt*fs;
  time = [t(1):dt:tfinal];

  % truncated SVD
  [U1, S1, V1] = svd(x, 'econ');

  if isempty(r)==1
    U = U1;
    V = V1;
    S = S1;
  else
    U = U1(:, 1:r);
    V = V1(:, 1:r);
    S = S1(1:r, 1:r);
  end

  % linear map s.t. y = Ax
  Atilde = U'*y*V*diag(1./diag(S));
  [eV, D] = eig(Atilde); mu = diag(D); omg = log(mu)/dt;
  mud = U*eV; y0 = mud\first_frame;

  % compute modes
  for iterate = 2:length(time)-1
    kip(:, iterate) = (y0.*exp(omg*time(iterate)));
  end

  sig = diag(S);
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% u_dmd = Phi*u_modes;
u_dmd = mud*kip;
% err_dmd = mean(mean(abs(real(u_dmd)-b_th180_T)));
% X_DMD = Phi*X_inter;
% rel_err_dmd = err_dmd/err_svd;

figure(11)
surf(real(u_dmd))
xlabel('t')
ylabel('x')
title('DMD generated data')
figure(12)
surf(real(u_dmd)-b_th180_T(:,1:end-1))
xlabel('t')
ylabel('x')
title('Error')

%%
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
    
    for iter=1:(size(b_th180_T,2)-1)
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

%%
% =========================================================================
% Information over the whole shot
% =========================================================================
clear b_th0

% b_th180 = [data{35} 0.5*data{35} + 0.5*data{14} data{14} data{42} data{46} data{21} data{50} data{54} data{29}];
node2 = cellstr(node_string); % Put the character array Node_string into a cell array

% Probe data at 180 degree location for all z
kk = 0;
b_th180 = zeros(length(time),8);
for zp = 5:5:45
    kk = kk+1;
    if zp == 10 % Broken probes at this location
        b_th180(:,kk) = zeros(length(time),1);
    else
        str = ['\b_p' num2str(zp) '_180_t']; % Name of probe for axial pos
        d = data{~cellfun('isempty', strfind(node2,str))};
        b_th180(:,kk) = data{~cellfun('isempty', strfind(node2,str))};
    end
end
b_th180(:,2) = 0.5*b_th180(:,1) + 0.5*b_th180(:,3); % Broken probe
b_th180 = b_th180' % So that rows are measurements and columns are time steps
% Moving window of data
k = 1;
for ii = 1:floor(length(time)-9)
    bb = b_th180(:,k:k+8);% window of data
    
    X1 = bb(:,1:end-1);
    X2 = bb(:,2:end);
    
    [U, Sigma, V] = svd(X1, 'econ');
    U = U(:,1:3);
    V = V(:,1:3);
    Sigma = Sigma(1:3,1:3);
    S = U'*X2*V*diag(1./diag(Sigma));
    [eV,D] = eig(S);
    mu = diag(D);
    omega = log(mu)/dt;
    Phi = U*eV;
    
    u = bb(:,1);
    y0 = Phi\u;
    
    for iter=1:(size(bb,2))
        u_modes(:,iter) = (y0.*exp(omega*dt*(iter-1)));
        % X_inter(:,t) = diag(exp(omega*dt*(t-1)))*b;
    end
    u_dmd = Phi*u_modes;
    
    % Error tracking between the 3 modes svd and dmd
    err_dmd(ii) = mean(mean(abs(real(u_dmd)-bb)));
    [Ut,St,Vt] = svd(bb,'econ');
    err_svd(ii) = mean(mean(abs(Ut(:,1:3)*St(1:3,1:3)*Vt(:,1:3)'-bb)));
    f(ii) = abs(imag(omega(2))/(2*pi)); % Frequency of oscillation
    I(ii) = real(omega(2)); % Intensity of oscillation?
    k = k+1;
end
figure(15)
subplot(2,1,1)
plot(time(vis_start:vis_end)*1e6,f(vis_start:vis_end))
title('FREQUENCY of oscillation')
axis([-inf inf -inf 5e6])
hold off
subplot(2,1,2)
% m1 data
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

figure(16)
plot(I)
title('Intensity of oscillation')

figure(17)
plot(I.*f)

figure(20)
plot(time(vis_start:vis_end)*1e6,err_dmd(vis_start:vis_end)./err_svd(vis_start:vis_end))
title('err dmd./err svd')
xlabel('time [us]')
axis([-inf inf 1 10])
% =========================================================================
% WAVELENGTH ANALYSIS
k = 1;
for ii = 1:floor(length(time)-9)
    bb = b_th180(:,k:k+8);% window of data
    bb = bb';
    
    X1 = bb(:,1:end-1);
    X2 = bb(:,2:end);
    
    [U, Sigma, V] = svd(X1, 'econ');
    U = U(:,1:3);
    V = V(:,1:3);
    Sigma = Sigma(1:3,1:3);
    S = U'*X2*V*diag(1./diag(Sigma));
    [eV,D] = eig(S);
    mu = diag(D);
    omega = log(mu)/dt;
    Phi = U*eV;
    
    u = bb(:,1);
    y0 = Phi\u;
    
    for iter=1:(size(bb,2))
        u_modes(:,iter) = (y0.*exp(omega*dt*(iter-1)));
        % X_inter(:,t) = diag(exp(omega*dt*(t-1)))*b;
    end
    u_dmd = Phi*u_modes;
    
    f(ii) = abs(imag(omega(2))/(2*pi)); % wavenumber of oscillation
    I(ii) = real(omega(2)); % Intensity of oscillation?
    k = k+1;
end
figure(18)
subplot(2,1,1)
plot(time(vis_start:vis_end)*1e6,f(vis_start:vis_end))
title('WAVENUMBER of oscillation')
axis([-inf inf -inf 5e6])
hold off
subplot(2,1,2)
% m1 data
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

figure(19)
plot(time(vis_start:vis_end)*1e6,I(vis_start:vis_end))
title('Intensity of oscillation')

figure(20)
plot(I.*f)





