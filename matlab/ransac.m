function [theta,cleanData,flag] = ransac(A,x,y,t,epsilon,mu)
%
% [theta,cleanData,flag] = ransac(A,x,y,t,epsilon,mu)
%
% [Input]:
% A:            3x(M^2) data matrix [x;y;t] with x,y,t 1x(M^2) vectors
% x:            x-coordinate of patch center point
% y:            y-coordinate of patch center point (=x)
% t:            time-stamp of patch center point
% epsilon:      estimated fraction of outliers [0,1]
% mu:           distance threshold
%
% [Output]:
% theta:        plane normal theta (= [0;0;1] if no inlier set found)
% cleanData:    inlier set. cleanData = [x;y;0] if not enough inliers found
% flag:         {0,1}. 1 if inlier set found, else 0
%
% This function estimates the plane parameters fitted for a data matrix
% based on the RANSAC algorithm.

flag = 0;

%% Preprocessing
N = size(A,2); % number of data points in the whole patch

% Initialization
theta0 = [0;0;1];   % initialize theta0
theta = theta0;     % set theta to default theta0
cleanData = zeros(4,N); % set cleanData to default zeros
cleanData(1:2,:) = A(1:2,:);

% remove zeros from data
A(:,A(3,:) == 0) = [];
nPoints = size(A,2); % new data size
if(nPoints < 3) % if less than 3 points -> no estimation possible
    cleanData = cleanData(1:3,:);
    return
end


%% RANSAC 

tmp = 0; % RANSAC constant to detect consensus set

% preallocation
sample = zeros(3,3); % sample matrix for ransac algorithm
z = zeros(1,nPoints);
dist = zeros(1,nPoints); 
cleanDist = zeros(1,nPoints);

% calculate number of iterations m with epsilon outliers and P = 0.99
nInlierEstimates = (1-epsilon)*N/2;
P = 0.99;
m = ceil(log(1-P)/log(1-(1-epsilon)^3));

% sample columns of A: (x;y;t)
sampleMat = datasample(A,2*m,2); 
%2*m because 3 points are needed for plane estimation: 2 random sampled
%points and patch center point

for i=1:m % RANSAC loop
        sample(1,:) = [x,y,t]; % patch center point has to be in sample!
        sample(2:3,:) = sampleMat(:,i:i+1)';
        % check colinearity of sample points
        if(iscolinear(sample(1,1:2),sample(2,1:2),sample(3,1:2)))
            %disp('Sample is colinear')
            theta0 = [0;0;1];
            theta = theta0;
        else
            % ordinary least squares solution
            theta0 = sample\[1;1;1];
        end

       % determine distances to estimated plane
       % make normal vector point upwards (i.e. theta0(3) >= 0 !)
       theta0 = sign(theta0(3))*theta0;
       % z value on plane for distance calculation
       d = theta0(1)*x+theta0(2)*y+theta0(3)*t;
       z = ((d-theta0(1).*A(1,:)-theta0(2).*A(2,:))./theta0(3));
       % square distances from plane
       dist = ((A(3,:)-z.*ones(1,size(A,2))).*theta0(3)).^2;
       % delete points too far away from plane
       cleanDist = dist;
       cleanDist(cleanDist > mu) = [];
       nInliers = length(cleanDist);
       if(nInliers > nInlierEstimates) % number of Inliers threshold
           if(nInliers > tmp) % take sample with biggest number of inliers
               theta = theta0;
               tmp = nInliers;
               cleanData = [A;dist];
               flag = 1;
           end
       end    
end

% remove outliers
cleanData(:,cleanData(4,:) > mu) = [];
% remove additional row containing distances from points to plane
cleanData = cleanData(1:3,:);

end
