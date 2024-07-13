%% Plot Vertical profiles
clear
close all
addpath("Monthdata\Aprilmonth")
addpath("Monthdata\Marchmonth")
addpath("Monthdata\Februarymonth\")
addpath("Weekdata\Maystrips\")

load("may1Wpeaks.mat")
load("aprpeaks.mat")
load("marpeaks.mat")
load("febpeaks.mat")

figure(1)
title('\bf Vertical profiles with integrated intensity', fontsize=16)
subplot(2,2,1)
plot(febpeaks.maxI,febpeaks.alt, '.')
title({'\bf February 8th to 28th'},fontsize=13)
xlabel('Intensity 10^{13}/ (nm \cdot m^2 \cdot str \cdot s)', FontSize=12)
ylabel('Altitude (km)', FontSize=12)
ylim([90 115]);
grid on

subplot(2,2,2)
plot(marpeaks.maxI,marpeaks.alt, '.')
title({'\bf March'},fontsize=13)
xlabel('Intensity 10^{13}/ (nm \cdot m^2 \cdot str \cdot s)', FontSize=12)
ylabel('Altitude (km)', FontSize=12)
ylim([90 115]);
grid on

subplot(2,2,3)
plot(aprpeaks.maxI,aprpeaks.alt, '.')
title({'\bf April'},fontsize=13)
xlabel('Intensity 10^{13}/ (nm \cdot m^2 \cdot str \cdot s)', FontSize=12)
ylabel('Altitude (km)', FontSize=12)
ylim([90 115]);
grid on

subplot(2,2,4)
plot(may1Wpeaks.maxI,may1Wpeaks.alt, '.')
title({'\bf May 1st to 8th'},fontsize=13)
xlabel('Intensity 10^{13}/ (nm \cdot m^2 \cdot str \cdot s)',FontSize=12)
ylabel('Altitude (km)', FontSize=12)
ylim([90 115]);
grid on

%% PLot intensity vs MLT
clear
close all
addpath("Monthdata\Aprilmonth")
addpath("Monthdata\Marchmonth")
addpath("Monthdata\Februarymonth\")
addpath("Weekdata\Maystrips\")

load("may1Wpeaks.mat")
load("aprpeaks.mat")
load("marpeaks.mat")
load("febpeaks.mat")

figure(1)
title('\bf Integrated intensity vs Magnetic Local Time', fontsize=16)
subplot(2,2,1)
plot(febpeaks.MLT,febpeaks.maxI, '.')
title({'\bf February 8th to 28th'},fontsize=13)
ylabel('Intensity 10^{13}/ (nm \cdot m^2 \cdot str \cdot s)',FontSize=12)
xlabel('MLT', FontSize=12)
xlim([0 24])
xticks([0,4,8,12,16,20,24])
grid on

subplot(2,2,2)
plot(marpeaks.MLT,marpeaks.maxI, '.')
title({'\bf March'},fontsize=13)
ylabel('Intensity 10^{13}/ (nm \cdot m^2 \cdot str \cdot s)',FontSize=12)
xlabel('MLT', FontSize=12)
xlim([0 24])
xticks([0,4,8,12,16,20,24])
grid on

subplot(2,2,3)
plot(aprpeaks.MLT,aprpeaks.maxI, '.')
title({'\bf April'},fontsize=13)
ylabel('Intensity 10^{13}/ (nm \cdot m^2 \cdot str \cdot s)',FontSize=12)
xlabel('MLT', FontSize=12)
xlim([0 24])
xticks([0,4,8,12,16,20,24])
grid on

subplot(2,2,4)
plot(may1Wpeaks.MLT,may1Wpeaks.maxI, '.')
title({'\bf May 1st to 8th'},fontsize=13)
ylabel('Intensity 10^{13}/ (nm \cdot m^2 \cdot str \cdot s)',FontSize=12)
xlabel('MLT', FontSize=12)
xlim([0 24])
xticks([0,4,8,12,16,20,24])
grid on


%% All peaks
clear
close all
addpath("Monthdata\Aprilmonth")
addpath("Monthdata\Marchmonth")
addpath("Monthdata\Februarymonth\")
addpath("Weekdata\Maystrips\")

load("may1Wpeaks.mat")
load("aprpeaks.mat")
load("marpeaks.mat")
load("febpeaks.mat")
allAlt=[aprpeaks.alt marpeaks.alt febpeaks.alt may1Wpeaks.alt] ;
allI = [aprpeaks.maxI marpeaks.maxI febpeaks.maxI may1Wpeaks.maxI] ;
allMLT = [aprpeaks.MLT marpeaks.MLT febpeaks.MLT may1Wpeaks.MLT] ;
allMlat =[aprpeaks.Mlat marpeaks.Mlat febpeaks.Mlat may1Wpeaks.Mlat] ;
figure(1)
plot(allI,allAlt, '.')
title({'\bf All peaks'},fontsize=13)
xlabel('Intensity 10^{13}/ (nm \cdot m^2 \cdot str \cdot s)',FontSize=12)
ylabel('Altitude (km)', FontSize=12)
ylim([90 115]);
grid on

figure(2)
plot(allMLT,allI, '.')
title({'\bf All peaks'},fontsize=13)
ylabel('Intensity 10^{13}/ (nm \cdot m^2 \cdot str \cdot s)',FontSize=12)
xlabel('MLT', FontSize=12)
xlim([0 24])
xticks([0,4,8,12,16,20,24])
grid on
