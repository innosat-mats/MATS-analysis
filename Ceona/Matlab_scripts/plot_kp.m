%% plot kp-index for whole period
% From https://kp.gfz-potsdam.de/en/data
clear
[time, value, status] = getKpindex('2023-02-01', '2023-02-28', 'Kp') ;
days = 28 ;

% Reshape yourValues into a matrix
KPdata = reshape(value, 8, days);

figure(1)
bar(time,value)
sgtitle({'\bf Kp-index vs time';'\rm Period: February'},fontsize=15)
ylabel('Kp', FontWeight='bold', FontSize=12)
xlabel('Time', FontWeight='bold', FontSize=12)
grid on

%% Plot parameter coorelations
clear
close all
addpath("Monthdata\")
load("aprpeaks.mat")

strips = aprpeaks ;

figure(1)
plot(strips.kp,strips.maxI , '.')
sgtitle({'\bf Intensity vs Kp-level';'\rm Period: April Northern Hemisphere'},fontsize=15)
ylabel('Intensity 10^{13}/ (nm \cdot m^2 \cdot str \cdot s)',FontWeight='bold', FontSize=12) ;
xlabel('Kp-level', FontWeight='bold', FontSize=12)
grid on

figure(2)
plot(strips.kp,strips.alt , '.')
sgtitle({'\bf Altitude vs Kp-level';'\rm Period: April Northern Hemisphere'},fontsize=15)
ylabel('Altitude (km)', FontWeight='bold', FontSize=12)
xlabel('Kp-level', FontWeight='bold', FontSize=12)
grid on

figure(3)
plot(strips.kp,strips.Mlat , '.')
sgtitle({'\bf MLat vs Kp-level';'\rm Period: April Northern Hemisphere'},fontsize=15)
ylabel('Geomagnetic latitude (degrees)', FontWeight='bold', FontSize=12)
xlabel('Kp-level', FontWeight='bold', FontSize=12)
grid on

figure(4)
plot(strips.alt,strips.maxI , '.')
sgtitle({'\bf Altitude vs Intensity';'\rm Period: April Northern Hemisphere'},fontsize=15)
ylabel('Intensity 10^{13}/ (nm \cdot m^2 \cdot str \cdot s)', FontWeight='bold', FontSize=12)
xlabel('Altitude (km)', FontWeight='bold', FontSize=12)
grid on
