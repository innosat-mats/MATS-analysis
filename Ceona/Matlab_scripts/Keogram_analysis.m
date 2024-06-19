
%% Plot single image
close
clear all
load("mar31_00_58IR1.mat")
clims = [min(min(ccdimage)) max(max(ccdimage))];
figure(1)
imagesc(ccdimage,clims)
ax1 = gca;
ax1.YDir = 'normal' ;
title('March 31st Image 00:58:38 ')
hold on

peakrow = 160 ;
crop = ccdimage(peakrow-4:peakrow+3,15:30) ;

%add central column line
xline([15 30], Color=[1 0 0], LineWidth=2)
hold on
yline([peakrow-4 peakrow+3],Color=[1 0 0], LineWidth=2);

S = sum(crop, "all") ;
numpixels = 128 ;
av = S/numpixels
%% Plot a keogram and a corresponding position plot below (ex altitude)
keogramgrad = keogram15feborb11 ;
clims = [min(min(keogramgrad)) max(max(keogramgrad))];
figure(1)
subplot(2,1,1) ;
imagesc(keogramgrad,clims)
ax1 = gca;
ax1.YDir = 'normal' ;
title('Keogram 3 Mars orbit 10')

dateStrings = time3mars; 
dateTimes = datetime(dateStrings, 'Format','dd/MM HH:mm:ss');  %dd/MM

subplot(2,1,2)
plot(alt3mars, '.') ;
title('Altitude 15 feb')
ylabel('Altitude (km)')
grid on
ylim([80,115])
xtickangle(45);
linkaxes([ax1 gca], 'x');  % Link x-axes of the two subplots

%%  Plotting the pixel values of different row of a keogram.
close all 
addpath('OneOrbit\')
load('mar2orb10Row170.mat')
load('mar2orb10Row150.mat')
load('mar2orb10Row130.mat')
load('mar2orb10regRow170.mat')
load('mar2orb10regRow150.mat')
load('mar2orb10regRow130.mat')
load('mar2orb10keogram.mat')
load('mar2orb10keograd.mat')
image = keogram;

[xData, yData] = prepareCurveData( [], Row130 );

% Fit model to data.
f = fit( xData, yData, 'poly1' );

% Plot fit with data.
figure(1);
plot(f);
hold on
plot(image(130,:),'DisplayName','Row 130',Marker='.',LineStyle='none')
hold on
plot(image(140,:),'DisplayName','Row 140',Marker='.',LineStyle='none')
hold on
plot(image(150,:),'DisplayName','Row 150',Marker='.',LineStyle='none')
hold on
plot(image(155,:), 'DisplayName','Row 155',Marker='.',LineStyle='none')
hold on
plot(image(160,:),'DisplayName','Row 160',Marker='.',LineStyle='none')
xlim([0,230])
ylabel('Pixel value',FontSize=13)
xlabel('Column','FontSize',13)
legend 
grid on
title('Multiple rows plotted and a linear regression of row 130')

figure(2)
plot(Row150-regRow150,'DisplayName','New row 150',Marker='.',LineStyle='none')
hold on
plot(Row150,'DisplayName','Row 150',Marker='.',LineStyle='none')
hold on
plot(regRow150, 'DisplayName','Regression')
xlim([0,230])
ylabel('Pixel value', FontSize=13)
xlabel('Column', FontSize=13)
grid on
legend 
title('Polynomial regression of row 150 and normalized curve')

%% Plot keogram image and some of the polynomial regression rows
close all
addpath('OneOrbit\')
load('mar2orb10Row170.mat')
load('mar2orb10Row150.mat')
load('mar2orb10Row130.mat')
load('mar2orb10regRow170.mat')
load('mar2orb10regRow150.mat')
load('mar2orb10regRow130.mat')
load('mar2orb10keogram.mat')
load('mar2orb10keograd.mat')

figure(1)
subplot(2,1,1)
clims = [min(min(keograd)) max(max(keograd))];
imagesc(keograd,clims)
ax1 = gca;
ax1.YDir = 'normal' ;
title('Keogram 2nd March orbit 10')

subplot(2,1,2)
%plot(Row150-regRow150,'DisplayName','Row 150-diff')
hold on
plot(Row150,'DisplayName','Row 150',Marker='.',LineStyle='none')
hold on
plot(regRow150, 'DisplayName','Row 150-reg')
hold on
plot(Row130,'DisplayName','Row 130',Marker='.',LineStyle='none')
hold on
plot(regRow130,'DisplayName','Row 130reg')
hold on
%plot(Row130-regRow130, 'DisplayName','Row 130-diff')
hold on
%plot(Row170, 'DisplayName','170')
hold on
%plot(Row170-regRow170,'DisplayName','170-diff')

title('Plot showing Keogram rows')
ylabel('Pixel value')
xlabel('Column')
ylim([0,350])
xlim([0,226])
linkaxes([ax1 gca], 'x');  % Link x-axes of the two subplots
grid on
legend

figure(2)
clims = [min(min(keogram)) max(max(keogram))];
imagesc(keogram,clims)
ax1 = gca;
ax1.YDir = 'normal' ;
title('Keogram 2nd March orbit 10')
%% Check means and maximum values of specific keogram columns
 
col = keogram(151:end,106) ;
coltest = keogram(161:end,109);

sum(col)/length(col) 
max(col) 
sum(coltest)/length(coltest) 
max(coltest) 