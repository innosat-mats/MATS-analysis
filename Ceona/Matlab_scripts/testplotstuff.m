%% Calculate mean values and max min
addpath("Monthdata\")
load("aprpeaksSH.mat")
load("aprpeaksNH.mat")

%set files for peaks
peaksNH = aprpeaksNH ;
peaksSH = aprpeaksSH ;

meanSH = mean(peaksSH.Mlat) ;
meanNH = mean(peaksNH.Mlat) ;
minNH = min(peaksNH.Mlat)  ;
maxNH = max(peaksNH.Mlat) ;
minSH = min(peaksSH.Mlat) ;
maxSH = max(peaksSH.Mlat) ;


%% Plots altitude vs time (uniform and non-uniform axis)
close all
addpath("Weekdata\Februarystrips\")
load("feb3Wpeaks.mat")

figure(1)
dateStrings = feb3Wpeaks.time; 
dateTimes = datetime(dateStrings, 'Format', 'dd/MM HH:mm:ss');
plot(dateTimes,feb3Wpeaks.alt,'.')
title('15 feb, limit, non-uniform time spacing of x-axis')
ylabel('Altitude (km)')
grid on
ylim([100,116]) ;
xtickangle(45);

% Select a subset of date strings to display as x-axis ticks, for
% non-uniform x-axis
numTicks = 12;  % Number of ticks to display
selectedIndices = round(linspace(1, length(dateTimes), numTicks));
xticks(dateTimes(selectedIndices));
xticklabels(datestr(dateTimes(selectedIndices), 'HH:MM:SS'));

figure(2)
plot(feb3Wpeaks.alt,'.')
title('15 feb, limit, uniform spacing of x-axis')
ylabel('Altitude (km)')
grid on
ylim([100,116]) ;
xtickangle(45);

