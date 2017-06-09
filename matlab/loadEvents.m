function [events] =loadEvents(fileName)
% loads events from log files (*.aedat) and simulation files (any other
% ending)
% returns matrix of events: each row is an event in the form:
%[x y polarity timeStamp]
%
% the timestamps are relative (timestamp of 1st event is at 0)
%
% author: basil.huber@gmail.com

% find out, if its simulation or real data

[~,~,suffix] = fileparts(fileName);

if(strcmp(suffix, '.aedat'))
    [allAddr, allTs] = loadaerdat(fileName);
    events = zeros(size(allAddr,1),4);
    for i = 1:size(events,1)
        [x,y,pol] = extractRetinaEventsFromAddr(allAddr(i));
        events(i,:) = [x 127-y pol allTs(i)];  % y axis must be flipped to point downwards
    end
    events(:,4) = events(:,4)-events(1,4);  % make timestamps relative
    events(:,3) = (events(:,3)-0.75)*4;  % make polarity {-1,1} instead of {0.5,1}
elseif(strcmp(suffix, '.mat'))
    load(fileName);
else
    events = importdata(fileName);
    if(size(events,1) < 1)
        events = [];
        return 
    end
    events = sortrows(events,4);          % sort the rows by the timestamp
    events(:,3) = (events(:,3)+1)/2;
    events(:,4) = events(:,4)-events(1,4);  % make timestamps relative
    events(:,3) = (events(:,3)-0.75)*4;  % make polarity {-1,1} instead of {0.5,1}
end
