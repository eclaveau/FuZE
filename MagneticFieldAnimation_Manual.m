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
tstart = 61.239; tend = 0; timesteps = 110;


% Time range of display for I_P and m1/m0
vis_start = find(data{1}>(max(data{1})/1.9),1,'first'); % Where to current first reach 10% of max current
vis_end = find(data{1}>(max(data{1})/5),1,'last'); % Where the current last reach 10% of max current

% Animation of magnetic field
[m1, phi, m0, b_th180, t_pos_start, t_pos_end, z0] = main_animation(vis_start, vis_end, shotnum, time, data, tstart, tend, timesteps);
window_time = time(t_pos_start:t_pos_end);
%%
% =========================================================================
% Peak finding
% =========================================================================
% for k = 1:size(b_th180,2)
%     figure(2)
%     plot(b_th180(:,k))
%     axis([-Inf Inf 0.05 0.35])
%     grid on
%     pause
% end

% figure(3)
% surf(b_th180)
% xlabel('t')
% ylabel('x')

m1m0 = m1./m0;
figure(4)
plot(m1m0(:,2),'r')
hold on
plot(m1m0(:,3),'b')
plot(m1m0(:,4),'g')
plot(m1m0(:,5),'k')
plot(m1m0(:,6),'m')
plot(m1m0(:,7),'y')
plot(m1m0(:,8),'.-k')
plot(m1m0(:,9),'.-r')
plot(m1m0(:,10),'.-b')
legend('p5', 'p10', 'p15', 'p20', 'p25', 'p30', 'p35', 'p40', 'p45')


phid = rad2deg(phi);
figure(5)
plot(phid(:,2),'r')
hold on
plot(phid(:,3),'b')
plot(phid(:,4),'g')
plot(phid(:,5),'k')
plot(phid(:,6),'m')
plot(phid(:,7),'y')
plot(phid(:,8),'.-k')
plot(phid(:,9),'.-r')
plot(phid(:,10),'.-b')
legend('p5', 'p10', 'p15', 'p20', 'p25', 'p30', 'p35', 'p40', 'p45')

% Find peaks and corresponding locations
for ii = 1:10
    [pks{ii}, t_ind{ii}] = findpeaks(m1m0(:,ii),'MinPeakProminence',.025);
end
% t_ind = indice of the time vector where the peaks are.

% All the peaks in a single graph
figure(6)
hold on
for kk = 1:10
    plot(1e6*window_time(t_ind{kk}),5*(kk-1)*ones(1,length(t_ind{kk})),'o')
end
axis([1e6*time(t_pos_start), 1e6*time(t_pos_end), 0, 50])

% STREAK FINDING
K = 1; % Number of "streaks", aka, travelling peaks.
K_s = 1; % Number of points in the streak.
L = 1; % Peak number at one array
ii = 1; % Z location
end_cond = 0; % Condition to end looking for peaks.
while end_cond == 0;
    if sum(t_ind{ii+1}>t_ind{ii}(L))>=1 % There is a future peak in the next axial array
        
        
        t_ind_next = t_ind{ii+1}(find(t_ind{ii+1}>t_ind{ii}(L),1)); % indice of time of the next peak at the next axial array
        
        % Parameter for the streak for 2 points. [times, pos]
        if K_s ==1 % First point in the streak
            streak{K} = [t_ind{ii}(L), ii];
        end
        
        streak{K}(K_s + 1, 1) = t_ind_next;  
        streak{K}(K_s + 1, 2) = ii + 1;
        %     phid_loc = phid(locs{ii}(1),ii); % phase of the local peak
        %     phid_next = phid(loc_next,ii+1); % phase of the next location
        %     if phid_loc-30<phid_next && phid_next<phid_loc+30 % Range of phase of the two peak
        %     end
        %     L = L + 1; % Look for the next peak in the same array.
        
        ii = ii + 1; % Look at the next array
        L = find(t_ind{ii+1}>t_ind{ii}(L),1); % Position at which the next point is
        K_s = K_s + 1; % Number of points in the streak.
    else
        ii = ii+1; % If no future peak, look from the next axial array.
        K = K + 1
    end
    
    if ii == 10 % Reached the end of the arrays.
    end_cond = 1;
    end
end
% Streaks, line connecting consecutive peaks
for kk = 1:size(streak{K},1)
streak_time{K}(kk) = 1e6*window_time(streak{K}(kk,1));
streak_pos{K}(kk) = 5*(streak{K}(kk,2)-1);
end
plot(streak_time{K}, streak_pos{K})
xlabel('time [us]')
ylabel('position [cm]')




















