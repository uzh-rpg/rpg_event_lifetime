function histogram_figure(inputdata)
% 
% histogram_figure(inputdata)
%
% [Input]:
% inputData:    [x,y,p,t,vx,vy,t_disp,t_e]
%               data matrix calculated with calcVelocity.m
%
% This function creates histograms from the estimated velocity components
% and display time.


% remove initialization phase and end
inputdata = inputdata(end/4:end*3/4,:);
% remove zero entries
inputdata(inputdata(:,7) == 0,:) = [];


figure
% display time
data = inputdata;
%subplot(2,3,1), hist(round(data(:,7)),0:100:10e+05)
figure(1), hist(round(data(:,7)),0:100:10e+05)
xlim([0,2e+04])
ylim([0,3000])
set(gca,'XTick',0:1000:20000)
set(gca,'XTickLabel','0| |2| |4| |6| |8| |10| |12| |14| |16| |18| |20')
xlabel('Display Time [ms]')
ylabel('Number of counts')


% x velocity
data = inputdata;
data(data(:,5) == 0,:) = [];
%subplot(2,3,2), hist(data(:,5),-10e-03:10e-06:10e-03)
figure(2), hist(data(:,5),-10e-03:10e-06:10e-03)
xlim([-10e-04,10e-04])
xlabel('v_x [pix/ \mus]')
%ylabel('Number of counts')

% y velocity
data = inputdata;
data(data(:,6) == 0,:) = [];
%subplot(2,3,3), hist(data(:,6),-10e-03:10e-06:10e-03)
figure(3), hist(data(:,6),-10e-03:10e-06:10e-03)
xlim([-10e-04,10e-04])
xlabel('v_y [pix/\mus]')
ylabel('Number of counts')

end
