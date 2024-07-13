@@ -7,9 +7,9 @@ days = 30 ;
% Reshape yourValues into a matrix
KPdata = reshape(value, 8, days);

figure(1)
bar(time,value)
sgtitle({'\bf Kp-index vs time';'\rm Period: May 1st week'},fontsize=15)
ylabel('Kp', FontWeight='bold', FontSize=12)
xlabel('Time', FontWeight='bold', FontSize=12)
ylim([0 9])

@@ -17,51 +17,56 @@ grid on

%% Plot KP coorelations
clear
close all
addpath("Monthdata\Aprilmonth")
load("aprpeaks.mat")
addpath("Monthdata\Marchmonth")
load("marpeaks.mat")
addpath("Monthdata\Februarymonth\")
load("febpeaks.mat")
strips = aprpeaks ;

figure(1)
plot(strips.kp,strips.alt , '.')
sgtitle({'\bf Altitude vs Kp-level';'\rm Period: April Northern Hemisphere'},fontsize=15)
ylabel('Altitude (km)', FontWeight='bold', FontSize=12)
xlabel('Kp-level', FontWeight='bold', FontSize=12)
grid on

figure(2)
plot(strips.kp,strips.Mlat , '.')
sgtitle({'\bf MLat vs Kp-level';'\rm Period: April Northern Hemisphere'},fontsize=15)
ylabel('Geomagnetic latitude (degrees)', FontWeight='bold', FontSize=12)
xlabel('Kp-level', FontWeight='bold', FontSize=12)
grid on

figure(3)
plot(strips.kp,strips.maxI , '.')
sgtitle({'\bf Intensity vs Kp-level';'\rm Period: April Northern Hemisphere'},fontsize=15)
ylabel('Intensity 10^{13}/ (nm \cdot m^2 \cdot str \cdot s)',FontWeight='bold', FontSize=12) ;
xlabel('Kp-level', FontWeight='bold', FontSize=12)
grid on


%%
close all
%addpath("Monthdata\Marchmonth\")
addpath("Weekdata\Maystrips\")
%load("marpeaks.mat")
load("may1Wpeaks")
peaks = may1Wpeaks ;

startDate = datetime(2023, 5, 1,'Format','dd/MM HH:mm:ss');
endDate = datetime(2023, 5, 7, 'Format','dd/MM HH:mm:ss');
timeline = startDate:endDate;

figure(1)
plot(peaks.time,peaks.kp,'.')
sgtitle({'\bf Plot of peak points Kp-index vs time';'\rm Period: March month '},fontsize=16)
ylabel('Kp-index')
grid on