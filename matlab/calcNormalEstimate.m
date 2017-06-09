function [x_est,y_est,theta_est,t_est] = calcNormalEstimate(x,y,t,vx,vy,theta)
%
% [x_est,y_est,theta_est,t_est] = calcNormalEstimate(x,y,t,vx,vy,theta)
%
% [Input]:
% x:        x-coordinate of incoming event
% y:        y-coordinate of incoming event
% t:        time-stamp of incoming event
% vx:       estimated velocity in x-direction
% vy:       estimated velocity in y-direction
% t_disp:   estimated display-time of incoming event
% theta:    estimated local plane normal at position of incoming event
%
% [Output]:
% x_est:        predicted x-coordinates of next event in velocity direction
% y_est:        predicted x-coordinates of next event in velocity direction
% theta_est:    predicted plane normal at position of next event
% t_est:        predicted time-stamp of next event in velocity direction
%
% This function predicts the global position of the next (neighboring)
% pixel in direction of motion of the incoming event as well as the
% time-stamp of it's next incoming event.

% local position vector (from pixel to neighbor)
position = [-1 0;...
            -1 -1;...
            0 -1;...
            1 -1;...
            1 0;...
            1 1;...
            0 1;...
            -1 1;...
            -1 0];
v = sqrt(vx^2+vy^2);

% if incoming event considered to be noise -> no prediction
if(v == 0)
    nx = 0;
    ny = 0;
    theta_est = [0;0;1];
    t_est = 0;
else
    % local position of neighbor
    n = round(atan2d(vy,vx)/45)+5; % {1,2,3,...9} index in local position vector
    nx = position(n,1); %  = {-1,0,1}
    ny = position(n,2); %  = {-1,0,1}
    v = (vx*nx+vy*ny)/sqrt(nx^2+ny^2); % velocity component in normal direction
    
    theta_est = theta; % const. local velocity/plane normal in direction of velocity vector
    t_est = sqrt(nx^2+ny^2)/v+t; % distance between pixel centers divided by velocity
end

% global position of neighbor
x_est = x+(nx); 
y_est = y+(ny);

end
