function poisson
clear all;  
timeStepS = 0.009;                  
spikesPerS = 75;                   
durationS = 10.000;                  
times = 0:timeStepS:durationS;	
vt = rand(size(times));
spikes = makeSpikes(timeStepS, spikesPerS, durationS, 100);
rasterPlot(spikes, timeStepS);



durationS = 10000;                                       
spikes = makeSpikes(timeStepS, spikesPerS, durationS);  %  spike train
spikeTimes = find(spikes) * timeStepS * 1000;           %  times when spikes occurred (ms)
spikeIntervals = spikeTimes(2:length(spikeTimes)) - spikeTimes(1:length(spikeTimes) - 1);

% Histogram of the spike intervals, normalized to unit volume.

figure(2);                                              
binSize = 1;                                            % 1 ms bins
x = [1:binSize:100];
intervalDist = hist(spikeIntervals(spikeIntervals < 100), x);
intervalDist = intervalDist / sum(intervalDist) / binSize; % normalize by dividing by spike number
bar(x, intervalDist);


y = exppdf(x, 1 / (spikesPerS * timeStepS));            % exponential function, scaled to predicted max
axis([min(x) max(x) 0 max(y) * 1.1]);                   
xlabel('Interspike interval');
ylabel('Probability');
hold on;
plot(x, y, 'r');                                       
hold off;

function rasterPlot(spikes, timeStepS)

figure(1);

times = [0:timeStepS:timeStepS * (length(spikes) - 1)];
axes('position', [0.1, 0.1, 0.8, 0.8]);
axis([0, length(spikes) - 1, 0, 1]);
trains = size(spikes, 1); 
ticMargin = 0.01;                                      % gap between spike trains (full scale is 1)
ticHeight = (1.0 - (trains + 1) * ticMargin) / trains;

for train = 1:trains
    spikeTimes = find(spikes(train, :) == 1);
    yOffset = ticMargin + (train - 1) * (ticMargin + ticHeight);
    for i = 1:length(spikeTimes)
        line([spikeTimes(i), spikeTimes(i)], [yOffset, yOffset + ticHeight]);
    end
end

xlabel('Time (s)')
title('Raster plot of spikes');
end

function spikes = makeSpikes(timeStepS, spikesPerS, durationS, numTrains)

if (nargin < 4)
    numTrains = 1;
end
times = [0:timeStepS:durationS];
spikes = zeros(numTrains, length(times));
for train = 1:numTrains
    vt = rand(size(times));
    spikes(train, :) = (spikesPerS*timeStepS) > vt;
end
end
end