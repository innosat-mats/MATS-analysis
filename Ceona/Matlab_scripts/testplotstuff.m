%% Calculate mean values and max min for month data
addpath("Monthdata\Marchmonth\")
load("marpeaksSH.mat")
load("marpeaksNH.mat")

%set files for peaks
peaksNH = marpeaksNH ;
peaksSH = marpeaksSH ;

meanSH = mean(peaksSH.maxI) ;
meanNH = mean(peaksNH.maxI) ;
minNH = min(peaksNH.maxI);
maxNH = max(peaksNH.maxI) ;
minSH = min(peaksSH.maxI) ;
maxSH = max(peaksSH.maxI) ;

%% 
addpath("Weekdata\Maystrips\")
load("may1WpeaksSH.mat")
load("may1WpeaksNH.mat")

%set files for peaks
peaksNH = may1WpeaksNH ;
peaksSH = may1WpeaksSH ;

meanSH = mean(peaksSH.kp) ;
meanNH = mean(peaksNH.kp) ;
minNH = min(peaksNH.kp)  ;
maxNH = max(peaksNH.kp) ;
minSH = min(peaksSH.kp) ;
maxSH = max(peaksSH.kp) ;

%% Plots altitude vs time (uniform and non-uniform axis)
close all
addpath("Monthdata\Aprilmonth\")
load("aprpeaksSH.mat")
load("aprpeaks.mat")
peaks = aprpeaks ;

figure(1)
dateTimes = aprpeaksSH.time; 
plot(dateTimes,aprpeaksSH.alt,'.')
title('Non-uniform time spacing of x-axis')
ylabel('Altitude (km)')
grid on
ylim([90,116]) ;
xtickangle(45);

% Select a subset of date strings to display as x-axis ticks, for
% non-uniform x-axis
numTicks = 12;  % Number of ticks to display
selectedIndices = round(linspace(1, length(dateTimes), numTicks));
xticks(dateTimes(selectedIndices));
xticklabels(datestr(dateTimes(selectedIndices), 'dd/mm, HH:MM'));

% figure(2)
% plot(aprpeaksSH.alt,'.')
% title('Uniform spacing of the points')
% ylabel('Altitude (km)')
% grid on
% ylim([100,116]) ;
% xtickangle(45);

figure(2)
plot(peaks.time,peaks.alt,'.')
title('April uniform spacing of x-axis')
ylabel('Altitude (km)')
grid on
ylim([90,116]) ;
xtickangle(45);

startDate = datetime(2023, 4, 1,'Format','dd/MM HH:mm:ss');
endDate = datetime(2023, 4, 30, 'Format','dd/MM HH:mm:ss');
timeline = startDate:endDate;

numTicks = 12;  % Number of ticks to display
selectedIndices = round(linspace(1, length(timeline), numTicks));
xticklabels(datestr(timeline(selectedIndices), 'dd/mm'));
