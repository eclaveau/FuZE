function [t_pos_cam] = camera(ti,dt,tf)

% Input: camera(ti, dt, tf)
% ti - Initial time that the camera took a snapshot
% tf - last time of the snapshot
% dt - time between frames
% Output: Vector of the position of the snapshots

% t_camera = 54.99:0.1:80.5;
t_camera = ti:dt:tf;

t_up = [];
    t_low = [];
    t_p = [];
    for cc = 1:length(t_camera)
        cc
        tl = (t_camera(cc) - 0.020)*10^-6; % Time of photo - half a sampling period
        tu = (t_camera(cc) + 0.020)*10^-6; % Time of photo + half a sampling period
        t_up = time >= tl;
        t_low = time <= tu;
        t_p = t_up + t_low;
        t_pos_cam(cc) = find(t_p == 2);
        
        % %         w_cam_mon
        %     w_cam_mon = data{87};
        % %     Find the time at which the first frame was taken which is the point where
        % %     the w_cam_mon start to go to its maximum value
        %     posmax = find(w_cam_mon==max(w_cam_mon),1);
        %     posshot(1) = posmax - (11 - find(w_cam_mon(posmax-10:posmax)==0,1,'last'));
        %
        % %     Times at which picture have been taken 4 Mhz = every 0.25 micro sec =
        % %     every 5 time interval (0.05 micro sec)
        %     for s = 1:12
        %         posshot(s) = posshot(1) + 5*(s-1);
        %     end
        %
    end
end