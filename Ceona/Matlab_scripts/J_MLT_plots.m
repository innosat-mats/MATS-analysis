%% Polar MLT plots for the whole period
what = 'maxI'; % Change between 'kp' 'alt' 'maxI'

clear workspace
clc

%Files to plot the peak points
addpath("Monthdata\Februarymonth\")
load("febpeaksNH.mat");
load("febpeaksSH.mat") ;
load("feballstripsNH.mat")
load("feballstripsSH.mat")
addpath("Monthdata\Marchmonth\")
load("marpeaksNH.mat");
load("marpeaksSH.mat") ;
load("marallstripsNH.mat")
load("marallstripsSH.mat")
addpath("Monthdata\Aprilmonth\")
load("aprpeaksNH.mat");
load("aprpeaksSH.mat");
load("aprallstripsNH.mat")
load("aprallstripsSH.mat")

peaksNH_MLT =  [febpeaksNH.MLT,marpeaksNH.MLT,aprpeaksNH.MLT];
peaksSH_MLT =  [febpeaksSH.MLT,marpeaksSH.MLT,aprpeaksSH.MLT];
peaksNH_Mlat =  [febpeaksNH.Mlat,marpeaksNH.Mlat,aprpeaksNH.Mlat];
peaksSH_Mlat =  [febpeaksSH.Mlat,marpeaksSH.Mlat,aprpeaksSH.Mlat];
peaksNH_Mlon =  [febpeaksNH.Mlon,marpeaksNH.Mlon,aprpeaksNH.Mlon];
peaksSH_Mlon =  [febpeaksSH.Mlon,marpeaksSH.Mlon,aprpeaksSH.Mlon];
peaksNH_kp =  [febpeaksNH.kp,marpeaksNH.kp,aprpeaksNH.kp];
peaksSH_kp =  [febpeaksSH.kp,marpeaksSH.kp,aprpeaksSH.kp];
peaksNH_time =  [febpeaksNH.time,marpeaksNH.time,aprpeaksNH.time];
peaksSH_time =  [febpeaksSH.time,marpeaksSH.time,aprpeaksSH.time];
peaksNH_maxI =  [febpeaksNH.maxI,marpeaksNH.maxI,aprpeaksNH.maxI];
peaksSH_maxI =  [febpeaksSH.maxI,marpeaksSH.maxI,aprpeaksSH.maxI];
peaksNH_alt =  [febpeaksNH.alt,marpeaksNH.alt,aprpeaksNH.alt];
peaksSH_alt =  [febpeaksSH.alt,marpeaksSH.alt,aprpeaksSH.alt];

stripsNH_MLT = [feballstripsNH.MLT,marallstripsNH.MLT,aprallstripsNH.MLT];
stripsSH_MLT = [feballstripsSH.MLT,marallstripsSH.MLT,aprallstripsSH.MLT];
stripsNH_Mlat = [feballstripsNH.Mlat,marallstripsNH.Mlat,aprallstripsNH.Mlat];
stripsSH_Mlat = [feballstripsSH.Mlat,marallstripsSH.Mlat,aprallstripsSH.Mlat];


fig = figure(position=[20 20 1200 600]);
switch what
    case 'kp'
        sgtitle({'\bf MLT plot of peak points correlated with kp';'\rm Period: February 8th to April 30th'},fontsize=16)
        savename ='C:\Users\judit\OneDrive - KTH\MATS\Plots\kp_all.png';
    case 'maxI'
        sgtitle({'\bf MLT plot of peak points correlated with max Intensity';'\rm Period: February 8th to April 30th'},fontsize=16) 
        savename ='C:\Users\judit\OneDrive - KTH\MATS\Plots\maxI_all.png';
    case 'alt'
        sgtitle({'\bf MLT plot of peak points correlated with altitude';'\rm Period: February 8th to April 30th'},fontsize=16)
        savename ='C:\Users\judit\OneDrive - KTH\MATS\Plots\alt_all.png';
end
%%%% North Hemisphere polar plot
subplot(1,2,1);
%plot settings
s1 = polarscatter(stripsNH_MLT*(pi/12),stripsNH_Mlat,0.8,'filled','black','MarkerFaceAlpha',.1);
hold on

switch what
    case 'kp'
        p1 = polarscatter(peaksNH_MLT*(pi/12),peaksNH_Mlat,[],peaksNH_kp,'filled');
        set(gca,"CLim",[0 9])  %Altitude limit = 90 115,  kplimit = 0-9 intensity limit = 1*10^3 3*10^4
        cb = colorbar ;
        colormap jet
        cb.Label.String = 'Kp-index' ;
    case 'maxI'
        p1 = polarscatter(peaksNH_MLT*(pi/12),peaksNH_Mlat,[],peaksNH_maxI,'filled');
        set(gca,"CLim",[1*10^3 3*10^4]) 
        cb = colorbar ;
        colormap jet
        cb.Label.String = '10^{13} photons / (nm \cdot m^2 \cdot sr \cdot s)' ; 
    case 'alt'
        p1 = polarscatter(peaksNH_MLT*(pi/12),peaksNH_Mlat,[],peaksNH_alt,'filled');
        set(gca,"CLim",[90 115]) 
        cb = colorbar ;
        colormap jet
        cb.Label.String = 'Altitude (km)';
end
        
ax1 = gca;
% MLT (theta settings)
ax1.ThetaDir = 'counterclockwise' ;
ax1.ThetaZeroLocation = 'bottom' ;
ax1.ThetaTick = [0,45,90,135,180,225,270,315,360];
ax1.ThetaTickLabel = {'00';'03';'06';'09';'12';'15';'18';'21'} ;
% Latitude settings
ax1.RLim = [40, 90] ;
ax1.RDir = "reverse" ;
rtickformat(ax1,"degrees") 
rticks(40:10:90);
        
%Settings for the peak points data tips
p1.DataTipTemplate.DataTipRows(1) = dataTipTextRow("MLT",peaksNH_MLT);
p1.DataTipTemplate.DataTipRows(2) = dataTipTextRow("MLat",peaksNH_Mlat);
p1.DataTipTemplate.DataTipRows(3) = dataTipTextRow("Altitude",peaksNH_alt);
p1.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Kp",peaksNH_kp);
p1.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Time",peaksNH_time);
p1.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Intensity",peaksSH_maxI);
%Settings for the path data tips
s1.DataTipTemplate.DataTipRows(1) = dataTipTextRow("MLT",stripsNH_MLT);
s1.DataTipTemplate.DataTipRows(2) = dataTipTextRow("MLat",stripsNH_Mlat);
        
subtitle('Northern Hemisphere','FontWeight','bold',FontSize=12);

%%%% South Hemisphere polar plot
subplot(1,2,2) ;
s2 = polarscatter(stripsSH_MLT*(pi/12),stripsSH_Mlat,0.8,'filled','black','MarkerFaceAlpha',.1);
hold on
switch what
    case 'kp'
        p2 = polarscatter(peaksSH_MLT*(pi/12),peaksSH_Mlat,[],peaksSH_kp,'filled'); %change the last attribute as see fit
        set(gca,"CLim",[0 9])  %[0 9] , [2*10^3 4*10^4] , [90 115]
        cb = colorbar ;
        cb.Label.String = 'Kp-index' ;  %'10^{13} photons / (nm \cdot m^2 \cdot sr \cdot s)'; 'Altitude (km)';'Kp-index' ; 
    case 'maxI'
        p2 = polarscatter(peaksSH_MLT*(pi/12),peaksSH_Mlat,[],peaksSH_maxI,'filled'); %change the last attribute as see fit
        set(gca,"CLim",[1*10^3 3*10^4]) 
        cb = colorbar ;
        cb.Label.String = '10^{13} photons / (nm \cdot m^2 \cdot sr \cdot s)' ;
    case 'alt'
        p2 = polarscatter(peaksSH_MLT*(pi/12),peaksSH_Mlat,[],peaksSH_alt,'filled'); %change the last attribute as see fit
        set(gca,"CLim",[90 115])
        cb = colorbar ;
        cb.Label.String = 'Altitude (km)' ; 
end

ax2 = gca;
ax2.ThetaDir = 'counterclockwise' ;
ax2.ThetaZeroLocation = 'bottom' ;
ax2.ThetaTick = [0,45,90,135,180,225,270,315,360];
ax2.ThetaTickLabel = {'00';'03';'06';'09';'12';'15';'18';'21'} ;
ax2.RLim = [-90, -40] ;
ax2.RDir = "normal" ;
rtickformat(ax2,"degrees")
rticks(ax2,(-90:10:-40));
        
%Settings for the peak points data tips
p2.DataTipTemplate.DataTipRows(1) = dataTipTextRow("MLT",peaksSH_MLT);
p2.DataTipTemplate.DataTipRows(2) = dataTipTextRow("MLat",peaksSH_Mlat);
p2.DataTipTemplate.DataTipRows(3) = dataTipTextRow("Altitude",peaksSH_alt);
p2.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Kp",peaksSH_kp);
p2.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Time",peaksSH_time);
p2.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Intensity",peaksSH_maxI);
%Settings for the path data tips
s2.DataTipTemplate.DataTipRows(1) = dataTipTextRow("MLT",stripsSH_MLT);
s2.DataTipTemplate.DataTipRows(2) = dataTipTextRow("MLat",stripsSH_Mlat);

subtitle('Southern Hemisphere','FontWeight','bold',FontSize=12);

saveas(fig,savename)