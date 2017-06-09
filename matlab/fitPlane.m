function [theta,flag] = fitPlane(data,theta_prior,t_prior,epsilon,mu,REG_ON,SHOW_PLOT)
%
% theta = fitPlane(data,theta_prior,t_prior,epsilon,mu,SHOW_PLOT)
%
% [Input]:
% data:             NxN matrix containing data to fit plane to
% theta_prior:      prior estimate for plane normal
% t_prior:          prior estimate for event time
% epsilon:          RANSAC estimated fraction of outliers [0,1]
% mu:               RANSAC distance threshold
% REG_ON:           {0,1} boolean: Regularization ON (1) or OFF (0)
% SHOW_PLOT:        {0,1} boolean to show plane for illustration
%
% [Output]:
% theta:            [a;b;c] plane normal for ax+by+cz = 1 plane model
% flag:             {0,1} boolean if inlier set found with RANSAC
%                   flag = 1 if found, flag = 0 if not found
%
% Fit plane to data given in NxN matrix data using RANSAC and
% regularization. X and Y coordinates are the pixel positions in the NxN
% image frame (origo is in upper left corner of image), the value at this
% position corresponds to the Z-value. 




%% format data: NxN matrix -> 3x(N^2) vector
N = size(data,1);
x = ceil(N/2); % y = x since it's the center point of the patch
t = data(x,x); % time-stamp of incoming event


% create 3x(N^2) matrix with data: A=[x;y;t] for linear model theta'*A=B
dataVector = reshape(data,1,N^2);
A = zeros(3,N^2);
% rows of data patch are y-entries, columns are x-entries
[A(2,:), A(1,:)] = ind2sub([N,N],1:N^2);
A(3,:) = dataVector;


%% Outlier rejection: only take upper median of data, set other data to 0 
med = median(A(3,:));
A(3,A(3,:)<med) = 0;


%% run RANSAC to estimate inliers and plane normal
[theta0,inliers,flag] = ransac(A,x,x,t,epsilon,mu);
% flag = 1 if consensus set found with RANSAC


%% run regularization over inlier set (if inlier set found)
if(flag)
    if(REG_ON)
        % time prediction error
        error = abs(t-t_prior);
        % exponential regularization weight
        lambda = 9+100*exp(-0.005*error);
    else
        % no regularization:
        lambda = 0;
    end
    % Nonlinear minimization using theta0 from RANSAC as an initial
    % starting point
    options = optimset('Display','off'); %,'Algorithm','levenberg-marquardt');
    [theta] = lsqnonlin(@(x) regularization(x,inliers,theta_prior,lambda),theta0,[],[],options);
else
    theta = theta0; % if no inliers found, theta0 = [0;0;1] -> incoming event is noise
end


% show Plot
if(SHOW_PLOT)
    figure(1)
    % plot points
    scatter3(A(1,:),A(2,:),A(3,:),50,'filled','MarkerFaceColor','b')
    hold on
    scatter3(inliers(1,:),inliers(2,:),inliers(3,:),100,'filled','MarkerFaceColor','r')
    %surf(data)
    [Y,X] = meshgrid(1:N, 1:N);
    a = theta(1);
    b = theta(2);
    c = theta(3);
    d = a*x+b*x+c*t;
    Z = -(a.*X+b.*Y-d)./c;
    surf(X,Y,Z,'Facecolor','g')
    hold off
    view(3)
    xlim([1,N])
    ylim([1,N])
    zlim([0,2*t+1])
end

end


function F = regularization(theta,A,theta_prior,lambda)
A = A'; % transpose inlier set A for linear model A*theta = b
F = zeros(size(A,1)+3,1);
% ordinary least squares term in minimization
F(4:end) = (A*theta-ones(length(A),1));
% regularization term in minimization
if(theta_prior(3) == 0 || isequal(theta_prior,[0;0;1]))
    % if invalid prior (0) or no prediction available ([0;0;1])
    F(1:3) = 0; % ignore regularization term
else
    F(1:3) = sqrt(lambda)*(theta-theta_prior);
end
end
