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
    SaveLocal = 0;
    save(name); % Save workspace to a local file
end

Animate_and_save = 0;
% ======================================================================
% Visualization parameters
% ======================================================================
if Animate_and_save == 1
    % Starting and ending time parameters  62.239 for the animation
    tstart = 60.239; tend = 70; timesteps = 0;
    
    
    % Time range of display for I_P and m1/m0
    vis_start = find(data{1}>(max(data{1})/1.9),1,'first'); % Where to current first reach 10% of max current
    vis_end = find(data{1}>(max(data{1})/5),1,'last'); % Where the current last reach 10% of max current
    
    % Animation of magnetic field
    [m1, phi, m0, b_th180, t_pos_start, t_pos_end, z0] = main_animation(vis_start, vis_end, shotnum, time, data, tstart, tend, timesteps);
    save([name '_data_only'], 'm1', 'phi', 'm0', 'b_th180', 't_pos_start', 't_pos_end', 'z0')
else
    load([name '_data_only'])
end
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


phid = rad2deg(phi) + 180; % angle in degree from 0 to 360
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
s = 0; % Streak condition: If a streak is present, it is 1, if not, it is 0.
phi_range = 30; % Maximum phase difference between two consecutive peaks [degree]
while end_cond == 0;
    s = 0;
    
    if K_s == 1 % First point in the streak
        L_s = L; % Initial peak that you started
        ii_s = ii; % Initial array that you started from
        streak{K} = [t_ind{ii}(L), ii];
    end
   
    % Conditions for a sterak to exist    
    if  ii<10 && sum(t_ind{ii+1}>t_ind{ii}(L))>=1  
        % You are not starting at the last array
        % There is a future peak in the next axial array.
        
        phi_loc = phi(t_ind{ii}(L),ii); % phase of the local peak
        phi_next = phi(t_ind{ii+1}((find(t_ind{ii+1}>t_ind{ii}(L),1))),ii+1); % phase of the next location
        
        if abs(angdiff(phi_loc,phi_next))<=deg2rad(phi_range) % The future peak is in a range of phase angle
            s = 1;% There is a streak    
        end
    else
     % No future peaks (s = 0, no change)
    end
   
    switch s
        case 0 % There is no streak
            if L_s == length(t_ind{ii_s}) % If you are at the last peak of that array.
                ii = ii_s+1; % look from the next axial array.
                L = 1; % First point of the new array
            else % Else look for the next peak.
                L = L_s + 1; % Next peak
                ii = ii_s; % Initial position where you started
            end
            K = K + 1; % New streak
            K_s = 1; % number of points in the streak is restarted
            
        case 1 % There is a streak
            t_ind_next = t_ind{ii+1}(find(t_ind{ii+1}>t_ind{ii}(L),1)); % indice of time of the next peak at the next axial array
            
            % Parameter for the streak for 2 points. [times, pos]
            
            streak{K}(K_s + 1, 1) = t_ind_next;
            streak{K}(K_s + 1, 2) = ii + 1;

            %     if phid_loc-30<phid_next && phid_next<phid_loc+30 % Range of phase of the two peak
            %     end
            %     L = L + 1; % Look for the next peak in the same array.
            
            L = find(t_ind{ii+1}>t_ind{ii}(L),1); % Position at which the next point is
            ii = ii + 1; % Look at the next array
            K_s = K_s + 1; % Number of points in the streak.
    end
    
    ii_s
    if ii_s == 9 && L_s == length(t_ind{ii_s}) % Reached the end of the arrays.
        % Reached the last peak of the array before the last array
        end_cond = 1;
    end
    ii
end
% Streaks, line connecting consecutive peaks
nK = K; % Total number of streaks
for K = 1:nK
    if isempty(streak{K}) % Streak where no peaks were found.
        streak_time{K} = [];
        streak_pos{K} = [];
    end
    for kk = 1:size(streak{K},1)
        streak_time{K}(kk) = 1e6*window_time(streak{K}(kk,1));
        streak_pos{K}(kk) = 5*(streak{K}(kk,2)-1);
    end
    plot(streak_time{K}, streak_pos{K})
end

xlabel('time [us]')
ylabel('position [cm]')





















