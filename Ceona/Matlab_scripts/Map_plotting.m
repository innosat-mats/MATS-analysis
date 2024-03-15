%% Plot bubbles on square map, with longitude and latitude.
clf
clear all
clc
%Files to plot the peak points
addpath("Weekdata\Februarystrips\")
load("feb3Wpeaks.mat") ;
load('coastlines')

peaks = feb3Wpeaks ;

colormap default
Intensity = discretize(peaks.maxI,6,'categorical');
figure(1)
gb = geobubble(peaks.maxlat,peaks.maxlon,peaks.kp,ColorData=Intensity);
gb.SizeLegendTitle = 'Kp';
gb.Title = 'Aurora distribution 3 week February' ;
gb.ColorLegendTitle = 'Intensity 10E13 ph/m2/sr/nm/s'; %Sum of intensity values or specific
%color for certain day
geobasemap satellite

figure(2)
plot(feb3Wpeaks.kp,abs(feb3Wpeaks.Mlat), '.')
%% Create polar plots with magnetic local time and magnetic latitude coordinates.
clear workspace
clc
clf
%Files to plot the peak points
addpath("Weekdata\Maystrips\")
load("may1WpeaksNH.mat");
load("may1WpeaksSH.mat") ;

%Files with all strips to use for TP path
load("may1WallstripsNH.mat")
load("may1WallstripsSH.mat")

%set files for path
stripsNH = may1WallstripsNH ;
stripsSH = may1WallstripsSH;

%set files for peaks
peaksNH = may1WpeaksNH ;
peaksSH = may1WpeaksSH ;

type = 2 ;  %1 means only peaks, 2 both peaks and satellite path
MLTplotfunc(type,stripsNH,stripsSH,peaksNH,peaksSH)

%% Polar MLT plots for a whole month
clear workspace
clc
clf
%Files to plot the peak points
addpath("Monthdata\")
load("aprpeaksNH.mat");
load("aprpeaksSH.mat") ;

%set files for the peaks you want to plot
peaksNH = aprpeaksNH ;
peaksSH = aprpeaksSH ;

%set files for the path you want to plot
stripsNH = aprpeaksNH ;
stripsSH = aprpeaksSH;

type = 1 ;  %1 means only peaks are plotted, 2 both peaks and satellite path
MLTplotfunc(1,stripsNH,stripsSH,peaksNH,peaksSH)

