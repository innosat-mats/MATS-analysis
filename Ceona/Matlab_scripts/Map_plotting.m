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
%% Polar plots with Magnetic local time and magnetic latitude coordinates.
clear workspace
clc
clf
%Files to plot the peak points
addpath("Weekdata\Aprilstrips\")
load("apr3WpeaksNH.mat");
load("apr3WpeaksSH.mat") ;

%Files with all strips to use for TP path
load("apr3WallstripsNH.mat")
load("apr3WallstripsSH.mat")

%set files for path
stripsNH = apr3WallstripsNH ;
stripsSH = apr3WallstripsSH;

%set files for peaks
peaksNH = apr3WpeaksNH ;
peaksSH = apr3WpeaksSH ;

type = 1 ;  %1 means only peaks, 2 both peaks and satellite path
MLTplotfunc(1,stripsNH,stripsSH,peaksNH,peaksSH)

%% Polar MLT plots for a whole month
clear workspace
clc
clf
%Files to plot the peak points
addpath("Weekdata\Aprilstrips\")
load("aprpeaksNH.mat");
load("aprpeaksSH.mat") ;

%set files for peaks
peaksNH = aprpeaksNH ;
peaksSH = aprpeaksSH ;

%set files for path
stripsNH = aprpeaksNH ;
stripsSH = aprpeaksSH;

type = 1 ;  %1 means only peaks, 2 both peaks and satellite path
MLTplotfunc(1,stripsNH,stripsSH,peaksNH,peaksSH)

