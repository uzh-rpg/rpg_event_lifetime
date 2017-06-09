function [x_parent,y_parent,t_parent] = deleteParent(x,y,vx,vy,t_disp)
%
% [x_parent,y_parent,t_parent] = deleteParent(x,y,vx,vy,t_disp)
%
% [Input]:
% x:        x-coordinate of incoming event
% y:        y-coordinate of incoming event
% t:        time-stamp of incoming event
% vx:       estimated velocity in x-direction
% vy:       estimated velocity in y-direction
% t_disp:   estimated display-time of incoming event
%
% [Output]:
% x_parent: x-coordinate of last event in negative direction of motion
% y_parent: y-coordinate of last event in negative direction of motion
% t_parent: resetted display-time of last event
%
% This function calculates the global position of the last pixel in
% negative direction of motion of an incoming event and sets the
% display-time of this pixel to zero, if the event is not considered to be
% noise. 

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

% if incoming event considered to be noise -> no reset (t_disp stays)
if(v == 0)
    nx = 0;
    ny = 0;
    t_parent = t_disp;
else
    n = round(atan2d(vy,vx)/45)+5; % {1,2,3,...9} index in local position vector
    nx = position(n,1); %  = {-1,0,1}
    ny = position(n,2); %  = {-1,0,1}
    t_parent = 0; % reset display-time of last event
end

% global position of last event
x_parent = x-(nx); 
y_parent = y-(ny);

end
