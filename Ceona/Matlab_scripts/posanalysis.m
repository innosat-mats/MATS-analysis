%% Plot positions of strips, latitude, longitude and altitude vs time.
close all
clc
clear
addpath("Weekdata\Maystrips\")
load("may1WpeaksSH.mat")
figure(1)
strips = may1WpeaksSH ;
dateTimes = strips.time; % Assuming strips.time is in 'dd/mm, HH:MM' format

% Select a subset of date strings to display as x-axis ticks
numTicks = 12;  % Number of ticks to display
selectedIndices = round(linspace(1, length(dateTimes), numTicks));

% Altitude plot
ax1 = subplot(3,1,1) ;
p1 = plot(strips.alt,'.') ;
p1.DataTipTemplate.DataTipRows(1) = dataTipTextRow("Altitude",strips.alt);
p1.DataTipTemplate.DataTipRows(2) = dataTipTextRow("Time",strips.time);

title({'\bf Peak points in geodetic coordinates'; '\rm Period: May 1st to 7th'}, fontsize=16)
grid on
grid minor
ylabel('\bf Altitude (km)','FontSize',11)
ylim([95,110]);
set(gca,'xticklabel',[],'xtick',selectedIndices)

% Latitude plot
ax2 = subplot(3,1,2) ;
p2 = plot(strips.lat, '.') ;
p2.DataTipTemplate.DataTipRows(1) = dataTipTextRow("Latitude",strips.maxlat);
p2.DataTipTemplate.DataTipRows(2) = dataTipTextRow("Time",strips.time);
grid on
grid minor

ylabel('\bf Latitude (degrees)','FontSize',11)
ylim([-90,-40])
% Set the x-axis ticks to be the selected date strings
set(gca,'xticklabel',[],'xtick',selectedIndices)


% Longitude plot
ax3 = subplot(3,1,3) ;
plot(strips.maxlon, '.')
grid on
grid minor

ylabel('\bf Longitude (degrees)','FontSize',11)
ylim([-200,200])
% Set the x-axis ticks to be the selected date strings
ax3.XTick = selectedIndices;
ax3.XTickLabel= datestr(dateTimes(selectedIndices), 'dd/mm, HH:MM');
xtickangle(ax3, 45);


% %Kp plot
% ax3 = subplot(3,1,3) ;
% p3 = plot(strips.kp, '.') ;
% grid on
% ylabel('\bf Kp index','FontSize',11)
% % Set the x-axis ticks to be the selected date strings
% ax3.XTick = selectedIndices;
% ax3.XTickLabel= datestr(dateTimes(selectedIndices), 'dd/mm, HH:MM');
% xtickangle(ax3, 45);
% p3.DataTipTemplate.DataTipRows(1) = dataTipTextRow("Kp",strips.kp);
% p3.DataTipTemplate.DataTipRows(2) = dataTipTextRow("Time",strips.time);

linkaxes([ax1 ax2 ax3],'x');  % Link x-axes of the three subplots