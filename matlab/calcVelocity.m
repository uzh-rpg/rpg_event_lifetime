function [velEvents] = calcVelocity(allEvents,N,epsilon,mu,REG_ON,SHOW_PLOT)
%
% [velEvents] = calcVelocity(allEvents,N,epsilon,mu,dt,SHOW_PLOT)
%
% [Input]:
% allEvents:    [x,y,p,t] matrix containing address-events from DVS
% N:            Local Window Size NxN
% epsilon:      RANSAC estimated fraction of outliers [0,1]
% mu:           RANSAC distance threshold
% REG_ON:       Regularization ON (1) or OFF (0)
% SHOW_PLOT:    {0,1} boolean to show plot during algorithm
%
% [Output]:
% velEvents:    [x,y,p,t,vx,vy,t_disp,t_e]
% vx,vy:        Velocity-components of incoming event in [pixel/mus]
% t_disp:       Display time, based on velocity of incoming event in [mus]
% t_e:        	Time prediction error for incoming event
%
% This function calculates the velocity components and display time of
% events from the matrix allEvents. allEvents is generated using the
% file loadEvents.m and contains x,y,p,t of each event.
%
% These events are split up by polarity {-1,1} and processed separately. 
% Assuming local constant velocity, the planar velocity v = [vx vy] of each
% event is calculated in [pixel/mus] and therefore, the display time of a 
% pixel is 1/|v| [mus]. 


close all;
clc;

% N == odd number!
if(mod(N,2) == 0)
    disp('Error: N must be an odd number')
    return;
end


IMAGE_FRAME = [128,128]; % Imageframe

% preallocate space
surfPos = zeros(IMAGE_FRAME);
surfNeg = zeros(IMAGE_FRAME);
surface = zeros(IMAGE_FRAME);
theta_estimatedPos = zeros(IMAGE_FRAME(1),IMAGE_FRAME(2),3);
theta_estimatedNeg = zeros(IMAGE_FRAME(1),IMAGE_FRAME(2),3);
t_estimatedPos = zeros(IMAGE_FRAME);
t_estimatedNeg = zeros(IMAGE_FRAME);

theta = zeros(3,1);
theta_prior = theta;
data = zeros(N,N);
a = 0;
b = 0;
c = 0;
t_prior = 0;
vx = 0;
vy = 0;
t_disp = 0;


% set time to be relative -> starting at t=0
allEvents(:,4) = allEvents(:,4)-allEvents(1,4);
 
last = length(allEvents(:,4)); % index of last event
velEvents = zeros(last,8); % extended event data with velocity
velEvents(:,1:4) = allEvents;

for i = 1:last % loop trough input data
x = allEvents(i,1)+1;
y = allEvents(i,2)+1;
p = allEvents(i,3);
t = allEvents(i,4);
surface(y,x) = t;

    
%% Positive Events
if(p == -1)
    surfPos(y,x) = t; % assign incoming event to positive SAE
    % check if pixel in imageFrame, depending on window size N  
    if(((IMAGE_FRAME(1)-floor(N/2)>=x) && (x>=ceil(N/2))) && ((IMAGE_FRAME(2)-floor(N/2)>=y) && (y>=ceil(N/2))))
        % read out local NxN surface window
        data = surfPos((y-floor(N/2)):(y+floor(N/2)),...
                       (x-floor(N/2)):(x+floor(N/2)));
        % read predicted time and theta at incoming events' position (x,y)           
        theta_prior = [theta_estimatedPos(y,x,1);...
                       theta_estimatedPos(y,x,2);...
                       theta_estimatedPos(y,x,3)];
        t_est = t_estimatedPos(y,x);
        % estimate local plane normal theta
        [theta,flag] = fitPlane(data,theta_prior,t_est,epsilon,mu,REG_ON,SHOW_PLOT);
        
        if(flag) % if theta found
            % choose sign of normal s.t. it shows in positive z-direction
            theta = sign(theta(3))*theta;            
            a = theta(1);
            b = theta(2);
            c = theta(3);
            % calculate velocity components and display time
            vx = c*(-a)/(a^2+b^2);
            vy = c*(-b)/(a^2+b^2);
            t_disp = sqrt(a^2+b^2)/c;
            velEvents(i,5) = vx;
            velEvents(i,6) = vy;
            velEvents(i,7) = t_disp;            
            if(t_est ~= 0)  % if prediction available at this location
                velEvents(i,8) = t-t_est;   % time prediction error t_e
            else
                velEvents(i,8) = 0;
            end            
        else % if not enough inliers found (flag = 0): theta = [0;0;1]
           vx = 0;
           vy = 0;
           t_disp = 0;
           velEvents(i,5) = 0;  % v_x
           velEvents(i,6) = 0;  % v_y
           velEvents(i,7) = 0;  % t_disp
           velEvents(i,8) = 0;  % t_e
        end         
         % predict next pixel
         [x_estimate,y_estimate,theta_estimate,t_estimate] = calcNormalEstimate(x,y,t,vx,vy,theta);         
         theta_estimatedPos(y_estimate,x_estimate,1:3) = theta_estimate; 
         t_estimatedPos(y_estimate,x_estimate) = t_estimate;

    end
end

%% Negative Events
if(p == 1)
    surfNeg(y,x) = t; % assign incoming event to negative SAE
    % check if pixel in imageFrame, depending on window size N  
   if(((IMAGE_FRAME(1)-floor(N/2)>=x) && (x>=ceil(N/2))) && ((IMAGE_FRAME(2)-floor(N/2)>=y) && (y>=ceil(N/2))))
        % read out local NxN surface window
        data = surfNeg((y-floor(N/2)):(y+floor(N/2)),...
                       (x-floor(N/2)):(x+floor(N/2))); 
        % read predicted time and theta at incoming events' position (x,y) 
        theta_prior = [theta_estimatedNeg(y,x,1);...
                       theta_estimatedNeg(y,x,2);...
                       theta_estimatedNeg(y,x,3)];
        t_est = t_estimatedNeg(y,x);
        % estimate local plane normal theta
        [theta,flag] = fitPlane(data,theta_prior,t_est,epsilon,mu,REG_ON,SHOW_PLOT);        
        
        if(flag) % if theta found
            % choose sign of normal s.t. it shows in positive z-direction
            theta = sign(theta(3))*theta;            
            a = theta(1);
            b = theta(2);
            c = theta(3);
            % calculate velocity components and display time
            vx = c*(-a)/(a^2+b^2);
            vy = c*(-b)/(a^2+b^2);
            t_disp = sqrt(a^2+b^2)/c;
            velEvents(i,5) = vx;
            velEvents(i,6) = vy;
            velEvents(i,7) = t_disp;
            if(t_est ~= 0) % if prediction available at this location 
                velEvents(i,8) = t-t_est;  % time prediction error t_e
            else
                velEvents(i,8) = 0;
            end            
        else % if not enough inliers found (flag = 0): theta = [0;0;1]
           vx = 0;
           vy = 0;
           t_disp = 0;
           velEvents(i,5) = 0; % v_x
           velEvents(i,6) = 0; % v_y
           velEvents(i,7) = 0; % t_disp         
           velEvents(i,8) = 0; % t_e
        end
         % predict next pixel
         [x_estimate,y_estimate,theta_estimate,t_estimate] = calcNormalEstimate(x,y,t,vx,vy,theta);
         theta_estimatedNeg(y_estimate,x_estimate,1:3) = theta_estimate; 
         t_estimatedNeg(y_estimate,x_estimate) = t_estimate;
         
    end
end

% display progress in command window every 100 events
if mod(i, 100) == 0
    disp([num2str(round(100/last*i)) ' / 100%           (' num2str(i) ' / ' num2str(last) ')'])
end

end
