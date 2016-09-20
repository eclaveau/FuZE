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


% ======================================================================
% Visualization parameters
% ======================================================================

% Starting and ending time parameters  62.239 for the animation
tstart = 50.239; tend = 0; timesteps = 8;
% Time range of display for I_P and m1/m0
vis_start = find(data{1}>(max(data{1})/1.9),1,'first'); % Where to current first reach 10% of max current
vis_end = find(data{1}>(max(data{1})/5),1,'last'); % Where the current last reach 10% of max current

% Animation of magnetic field
[b_th180, t_pos_start, t_pos_end, z0] = main_animation(vis_start, vis_end, shotnum, time, data, tstart, tend, timesteps);

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
dt = time(2)-time(1);

rank = 2;
b_untrans = b_th180;
b_th180_T = round(1e4*b_th180')/(1e4);

% Three modes svd reconstruction
[Ut, Sigmat, Vt] = svd(b_th180_T, 'econ');
figure(7)
subplot(1,2,1)
surf(Ut(:,1:rank)*Sigmat(1:rank,1:rank)*Vt(:,1:rank)')
title(['USV reconstruction with ' num2str(rank) ' first modes'])
subplot(1,2,2)
surf(Ut(:,1:rank)*Sigmat(1:rank,1:rank)*Vt(:,1:rank)'-b_th180_T)
title(['Error between data and ' num2str(rank) ' modes reconstructions'])
err_svd = mean(mean(abs(Ut(:,1:rank)*Sigmat(1:rank,1:rank)*Vt(:,1:rank)'-b_th180_T)));

figure(8)
plot(diag(Sigmat)/sum(diag(Sigmat))*100,'.','LineWidth',3)
title('Principal Components [%]')

figure(9)
ii = 0;
for kk = [1 2 5 6]
    ii = ii +1;
    subplot(4,2,kk)
    plot(Ut(:,ii))
    ident = [num2str(ii) ' Spatial Mode'];
    title(ident)
    subplot(4,2,kk+2)
    plot(10^6*time(t_pos_start:t_pos_end),abs(Vt(:,ii)))
    ident = [num2str(ii) ' Temporal Mode'];
    title(ident)
end

X1 = b_th180_T(:,1:end-1);
X2 = b_th180_T(:,2:end);
u = b_th180_T(:,1);
% dt = 1;
% t = 0:9;

[U1, Sigma1, V1] = svd(X1, 'econ');
U = U1(:,1:rank);
V = V1(:,1:rank);
Sigma = Sigma1(1:rank,1:rank);
S = U'*X2*V*diag(1./diag(Sigma));
[eV,D] = eig(S);
mu = diag(D);
omega = log(mu)/dt;
Phi = U*eV;
y0 = Phi\u;

for iter=1:(size(b_th180_T,2)-1)
    u_modes(:,iter) = (y0.*exp(omega*dt*(iter-1)));

end
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%   raw = b_th180_T;
%   r = rank;
%   t = linspace(0,9,9);
%   fs = 0;
%   
%   x = raw(:, 1:end-1);
%   y = raw(:, 2:end);
%   first_frame = raw(:, 1);
% 
%   % finalize our timebase
%   dt = mean(diff(t));
%   tfinal = t(end) + dt*fs;
%   time = [t(1):dt:tfinal];
% 
%   % truncated SVD
%   [U1, S1, V1] = svd(x, 'econ');
% 
%   if isempty(r)==1
%     U = U1;
%     V = V1;
%     S = S1;
%   else
%     U = U1(:, 1:r);
%     V = V1(:, 1:r);
%     S = S1(1:r, 1:r);
%   end
% 
%   % linear map s.t. y = Ax
%   Atilde = U'*y*V*diag(1./diag(S));
%   [eV, D] = eig(Atilde); mu = diag(D); omg = log(mu)/dt;
%   mud = U*eV; y0 = mud\first_frame;
% 
%   % compute modes
%   for iterate = 1:length(time)-1
%     kip(:, iterate) = (y0.*exp(omg*time(iterate)));
%   end
% 
%   sig = diag(S);
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
u_dmd = Phi*u_modes;
% u_dmd = mud*kip;
% err_dmd = mean(mean(abs(real(u_dmd)-b_th180_T)));

% rel_err_dmd = err_dmd/err_svd;

figure(11)
surf(real(u_dmd))
ylabel('t')
xlabel('x')
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
for kk = 1:rank
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
    subplot(3,3,kk)
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
clear b_th0 u_modes

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
b_th180 = b_th180'; % So that rows are measurements and columns are time steps
b_th180 = round(1e4*b_th180)/(1e4); % Remove the rounding error instability

% Moving window of data
k = vis_start;
for ii = 1:vis_end-vis_start
    bb = b_th180(:,k:k+8);% window of data
    
    X1 = bb(:,1:end-1);
    X2 = bb(:,2:end);
    
    [U, Sigma, V] = svd(X1, 'econ');
    U = U(:,1:rank);
    V = V(:,1:rank);
    Sigma = Sigma(1:rank,1:rank);
    S = U'*X2*V*diag(1./diag(Sigma));
    [eV,D] = eig(S);
    mu = diag(D);
    omega = log(mu)/dt;
    Phi = U*eV;
    
    u = bb(:,1);
    y0 = Phi\u;
    
    for iter=1:(size(bb,2))-1
        u_modes(:,iter) = (y0.*exp(omega*dt*(iter-1)));

    end
    u_dmd = Phi*u_modes;
    
    % Error tracking between the 3 modes svd and dmd
    err_dmd(ii) = mean(mean(abs(real(u_dmd)-bb(:,1:end-1))));
%     [Ut,St,Vt] = svd(bb,'econ');
%     err_svd(ii) = mean(mean(abs(Ut(:,1:3)*St(1:3,1:3)*Vt(:,1:3)'-bb)));
    f(ii) = abs(imag(omega(2))/(2*pi)); % Frequency of oscillation
    I(ii) = real(omega(2)); % Intensity of oscillation?
    k = k+1;
end
figure(15)
subplot(2,1,1)
plot(time(vis_start:vis_end-1)*1e6,f)
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
% pos = get(gca, 'Position');
% % pos(1) = pos(1);
% % pos(2) = ymargin*pos(2);
% set(gca, 'Position', pos)
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
plot(time(vis_start:vis_end-1)*1e6,err_dmd)
title('err dmd')
xlabel('time [us]')
axis([-inf inf 0 0.00001])
% =========================================================================
% WAVELENGTH ANALYSIS
k = vis_start;
for ii = 1:vis_end-vis_start
    bb = b_th180(:,k:k+8)';% window of data
    
    X1 = bb(:,1:end-1);
    X2 = bb(:,2:end);
    
    [U, Sigma, V] = svd(X1, 'econ');
    U = U(:,1:rank);
    V = V(:,1:rank);
    Sigma = Sigma(1:rank,1:rank);
    S = U'*X2*V*diag(1./diag(Sigma));
    [eV,D] = eig(S);
    mu = diag(D);
    omega = log(mu)/dt;
    Phi = U*eV;
    
    u = bb(:,1);
    y0 = Phi\u;
    
    for iter=1:(size(bb,2))-1
        u_modes(:,iter) = (y0.*exp(omega*dt*(iter-1)));

    end
    u_dmd = Phi*u_modes;
    
    % Error tracking between the 3 modes svd and dmd
    err_dmd(ii) = mean(mean(abs(real(u_dmd)-bb(:,1:end-1))));
%     [Ut,St,Vt] = svd(bb,'econ');
%     err_svd(ii) = mean(mean(abs(Ut(:,1:3)*St(1:3,1:3)*Vt(:,1:3)'-bb)));
    f(ii) = abs(imag(omega(2))/(2*pi)); % Frequency of oscillation
    I(ii) = real(omega(2)); % Intensity of oscillation?
    k = k+1;
end
figure(21)
subplot(2,1,1)
plot(time(vis_start:vis_end-1)*1e6,f)
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
% pos = get(gca, 'Position');
% % pos(1) = pos(1);
% % pos(2) = ymargin*pos(2);
% set(gca, 'Position', pos)
axis([-inf inf -inf 0.5])
ax = gca;
set(ax,'YTick',[0,0.2,0.4]);
hold off

figure(22)
plot(I)
title('Intensity of oscillation')

figure(23)
plot(I.*f)

figure(24)
plot(time(vis_start:vis_end-1)*1e6,err_dmd)
title('err dmd')
xlabel('time [us]')