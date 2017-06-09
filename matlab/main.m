%% Event Lifetime Estimation
% E. Mueggler, C. Forster, N. Baumli, G. Gallego, D. Scaramuzza
% "Lifetime Estimation of Events from Dynamic Vision Sensors"
% In: IEEE International Conference on Robotics and Automation (ICRA), 
% Seattle, 2015.
% PDF: http://rpg.ifi.uzh.ch/docs/ICRA15_Mueggler.pdf

%% Load data
data_folder = 'data';

% Experiment 1: stripes dataset
%dataset = 'stripes.mat';

% Experiment 2: Garfield
%dataset = 'garfield.mat';

% Experiment 3: quadrotor flip
dataset = 'flip.mat';

% Experiment 4: building
%dataset = 'building.mat';

% load data
events = loadEvents([data_folder, '/', dataset]);

%% Parameters
% Window Size
N = 5;
% Estimated fraction of outliers for RANSAC algorithm
epsilon = 0.4;
% Euclidian distance threshold for RANSAC algorithm
mu = 0.0001;
% Regularization
reg = true;
% Visualization during computation
vis = false;
% Show velocity on visualization
show_vel = false;

%% Compute lifetime
events_with_lifetime = calcVelocity(events, N, epsilon, mu, reg, vis);

%% Create video
% visualization parameter for lifetime
if strcmp(dataset, 'flip.mat')
    cmax = 500;
else
    cmax = 14000;
end

% video using lifetime
dispOutput(events_with_lifetime, show_vel, -1, [dataset, '_lifetime'], cmax);
% video using fixed time interval of 30ms
dispOutput(events_with_lifetime, show_vel, 30, [dataset, '_dt30']);
