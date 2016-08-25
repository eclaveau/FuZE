function [t_pos,t_pos_start,t_pos_end] = startendtime(time,tstart,tend,timesteps)

% startendtime(tstart,tend,timesteps)
% time - time vector
% tstart - Time [us] of the start of the data
% tend - Time [us] of the end or 0 if using timesteps
% timesteps - number of time steps to use from tstart. 0 if using tend
% Ouput: [t_pos, t_pos_start, t_pos,end]
% t_pos - Vector of the time values
% t_pos_start - First time
% t_pos_end - end time

    
    if tend == 0; % mode where it calculate for a finite number of timesteps
    tend = tstart+timesteps*0.0401; % Ending time microsec
    else
    end
    
    tl = (tstart - 0.020)*10^-6; % Time of photo - half a sampling period
    tu = (tstart + 0.020)*10^-6; % Time of photo + half a sampling period
    t_up = time >= tl;
    t_low = time <= tu;
    t_p = t_up + t_low;
    t_pos_start = find(t_p == 2);
    
    tl = (tend - 0.020)*10^-6; % Time of photo - half a sampling period
    tu = (tend + 0.020)*10^-6; % Time of photo + half a sampling period
    t_up = time >= tl;
    t_low = time <= tu;
    t_p = t_up + t_low;
    t_pos_end = find(t_p == 2);
    
    t_pos = t_pos_start:t_pos_end;
%     
%     if SingleFig ==1
%     tsnap = 40;
% %     tsnap = [62.239 62.6792 62.8793 63.4795 64.52 65.2403];
%     
%     for kk = 1:length(tsnap)
%         tl = (tsnap(kk) - 0.020)*10^-6; % Time of snapshot - half a sampling period
%         tu = (tsnap(kk) + 0.020)*10^-6; % Time of snapshot + half a sampling period
%         t_up = time >= tl;
%         t_low = time <= tu;
%         t_p = t_up + t_low;
%         t_pos(kk) = find(t_p == 2);
%     end
% else
%     t_pos = t_pos_start:t_pos_end;
% end