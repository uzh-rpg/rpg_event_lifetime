function dispOutput(velEvents, SHOW_VELOCITY, dt, filename, cmax)
%
% dispOutput(velEvents,SHOW_PLOT)
%
% [Input]:
% velEvents:        [x,y,p,t,vx,vy,t_disp,t_e]
%                   data matrix calculated with calcVelocity.m
% SHOW_VELOCITY:    {0,1} boolean to display velocity vectors in movie
% dt:               time interval [ms] for visualization, 
%                   negative value for lifetime
% filename:         filename without extension
% cmax:             visualization parameter
%
%
% This function  creates a movie based on the data in the matrix velEvents. 
% This data is generated using the function calcVelocity.m For the movies, 
% the incoming events are accumulated over a time interval of 1ms and
% display using a framerate of 100 fps, i.e. the movie is 10x slower than
% real-time. The Frames are displayed in a figure and captured for the
% movie.

if nargin < 5
    cmax = 10000;
end

close all;
clc;

iptsetpref('ImshowBorder','tight'); % remove borders around figure
s = 1e+05; % scale factor for quiver plot


% set time to be relative -> starting at t=0
velEvents(:,4) = velEvents(:,4)-velEvents(1,4);

endtime = round(velEvents(end,4)/1000); % endtime in ms
IMAGE_FRAME = [128,128]; % image frame
timemat = zeros(IMAGE_FRAME); % initialize display time matrix
colormat = zeros(IMAGE_FRAME); % initialize color matrix
velocitymat = zeros(IMAGE_FRAME(1),IMAGE_FRAME(2),3); % initialize velocity vector matrix

velEvents(:,4) = ceil(velEvents(:,4)./1000); % round data to ms
velEvents(:,7) = ceil(velEvents(:,7)./1000);

data = velEvents; 

% create linear colormap red to yellow
cmap = colormap(autumn);
cmap = zeros(length(cmap)+2,3);
cmap(2:end-1,:) = colormap(autumn);
caxis([0, cmax]);
cmap(end,:) = [1,1,1];
cmap(1,:) = [0,0,0];

% Set up the movie: try mp4, otherwise use avi
fmtList = VideoReader.getFileFormats();
if any(ismember({fmtList.Extension}, 'mp4'))
    writerObj = VideoWriter([filename, '.mp4'],'MPEG-4');
else
    writerObj = VideoWriter([filename, '.avi']);
end
writerObj.FrameRate = 100; % How many frames per second.
open(writerObj); 

%% event-based full frame visualization algorithm
if(dt < 0)
    
    for i = 0:endtime % run through data in milliseconds steps
        tmp = data;
        tmp = tmp(tmp(:,4)==i,:); % only take events at this particular time-stamp
        for l = 1:length(tmp(:,4)) % loop through these events and fill matrices
            x = tmp(l,1)+1;
            y = tmp(l,2)+1;
            vx = tmp(l,5);
            vy = tmp(l,6);
            t_disp = tmp(l,7);
            timemat(y,x) = t_disp; % assign integration time to event
            colormat(y,x) = t_disp; % assign color value to event
            velocitymat(y,x,1) = vx; % x velocity
            velocitymat(y,x,2) = vy; % y velocity
            % calculate global position of pixel in negative direction of motion
            [x_parent,y_parent,t_parent] = deleteParent(x,y,vx,vy,t_disp);
            % reset pixel in negative direction of motion
            timemat(y_parent,x_parent) = t_parent;
        end
        imagemat = logical(timemat).*colormat;
        % create white border around image
        imagemat(1,1:end) = cmax;
        imagemat(end,1:end) = cmax;
        imagemat(1:end,1) = cmax;
        imagemat(1:end,end) = cmax;
        % color image
        image = ind2rgb((imagemat),cmap);
        % black white image
        %image = ceil(mat2gray(timemat));
        
        % increase size of image by factor 3
        image = imresize(image,3,'method','nearest');
        
        if(SHOW_VELOCITY) % calculate velocity vectors
            xvelocitymat = logical((timemat)).*velocitymat(:,:,1);
            yvelocitymat = logical((timemat)).*velocitymat(:,:,2);
            xvel = reshape(xvelocitymat,1,IMAGE_FRAME(1)^2);
            yvel = reshape(yvelocitymat,1,IMAGE_FRAME(1)^2);
            timemat2 = timemat;
            % only take every 2nd vector
            timemat2(1:2:end,1:2:end) = 0;
            timemat2(2:2:end,2:2:end) = 0;
            flag = reshape(timemat2,1,IMAGE_FRAME(1)^2);
            A = zeros(5,IMAGE_FRAME(1)^2);
            [A(2,:), A(1,:)] = ind2sub([IMAGE_FRAME(1),IMAGE_FRAME(2)],1:IMAGE_FRAME(1)^2);
            A(3,:) = xvel;
            A(4,:) = yvel;
            A(5,:) = flag;
            A(:,A(5,:)==0) = [];
        end
        % write image to figure
        figure(1)
        imshow(image)
        hold on
        if(SHOW_VELOCITY)
            quiver(3*(A(1,:)-1),3*(A(2,:)-1),s*A(3,:),s*A(4,:),0,'w')
        end
        % display information in video
        text(10,10,[num2str(i) ' ms   |   N = 5   |   regularization ON'],'Color','g')
        hold off
        frame = getframe(gcf);
        writeVideo(writerObj, frame);
        
        % display progress in command window
        if mod(i, 100) == 0
            disp(['time = ' num2str(i) ' / ' num2str(endtime)])
        end
        
        timemat = timemat-1; % decrease integration time by 1ms
        timemat(timemat<0) = 0; % no negative times
    end

%% algorithm with fixed integration time dt
else
    for i = 0:endtime % run through data in milliseconds steps
        tmp = data;
        tmp = tmp(tmp(:,4)==i,:); % only take events at this particular time-stamp
        for l = 1:length(tmp(:,4))
            x = tmp(l,1)+1;
            y = tmp(l,2)+1;
            timemat(y,x) = dt; % assign integration time to event
        end
        % create white border around image
        timemat(1,1:end) = 0.5;
        timemat(end,1:end) = 0.5;
        timemat(1:end,1) = 0.5;
        timemat(1:end,end) = 0.5;
        image = ceil(mat2gray(timemat)); % gray value image
        % increase size of image by factor 3
        image = imresize(image,3,'method','nearest');
        
        % write image to figure
        figure(1)
        imshow(image)
        hold on
        % display information in video
        text(10,10,[num2str(i) ' ms   |   integration time = 30ms'],'Color','g')
        hold off
        frame = getframe(gcf);
        writeVideo(writerObj, frame);
        writeVideo(writerObj, image);
        
        % display progress in command window 
        if mod(i, 100) == 0
            disp(['time = ' num2str(i) ' / ' num2str(endtime)])
        end
        
        timemat = timemat-1; % decrease integration time by 1ms
        timemat(timemat<0) = 0; % no negative times
    end
end

% close movie object
close(writerObj);

end
