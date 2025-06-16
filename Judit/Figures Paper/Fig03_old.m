%% 

clear workspace
clc

savename ='\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub datafiles\MATS\alt_MLT.png';
savefig = '\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub datafiles\MATS\Fig03.png';
savefigfeb = '\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub datafiles\MATS\Fig03feb.png';
savefigap = '\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub datafiles\MATS\Fig03apr.png';
savefigmar = '\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub datafiles\MATS\Fig03mar.png';
seph = 1;

c1= [0.0 0.3 0.9];
c2= [0.0 0.8 0.4];
c3= [0.9 0.2 0.5];


%Files to plot the peak points
addpath("\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Ceona\Matlab_scripts\Monthdata\Februarymonth\")
load("febpeaksNH.mat");
load("febpeaksSH.mat") ;
load("feballstripsNH.mat")
load("feballstripsSH.mat")
addpath("\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Ceona\Matlab_scripts\Monthdata\Marchmonth\")
load("marpeaksNH.mat");
load("marpeaksSH.mat") ;
load("marallstripsNH.mat")
load("marallstripsSH.mat")
addpath("\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Ceona\Matlab_scripts\Monthdata\Aprilmonth\")
load("aprpeaksNH.mat");
load("aprpeaksSH.mat");
load("aprallstripsNH.mat")
load("aprallstripsSH.mat")

peaksNH_MLT =  [febpeaksNH.MLT,marpeaksNH.MLT,aprpeaksNH.MLT];
peaksSH_MLT =  [febpeaksSH.MLT,marpeaksSH.MLT,aprpeaksSH.MLT];
peaksNH_kp =  [febpeaksNH.kp,marpeaksNH.kp,aprpeaksNH.kp];
peaksSH_kp =  [febpeaksSH.kp,marpeaksSH.kp,aprpeaksSH.kp];
peaksNH_maxI =  [febpeaksNH.maxI,marpeaksNH.maxI,aprpeaksNH.maxI];
peaksSH_maxI =  [febpeaksSH.maxI,marpeaksSH.maxI,aprpeaksSH.maxI];
peaksNH_alt =  [febpeaksNH.alt,marpeaksNH.alt,aprpeaksNH.alt];
peaksSH_alt =  [febpeaksSH.alt,marpeaksSH.alt,aprpeaksSH.alt];

% Separate kp index in 3 groups for each hemisphere
in03NH = find(peaksNH_kp<=3);
in36NH = find(peaksNH_kp>3 & peaksNH_kp<=6);
in69NH = find(peaksNH_kp>6 & peaksNH_kp<=9);
in03SH = find(peaksSH_kp<=3);
in36SH = find(peaksSH_kp>3 & peaksSH_kp<=6);
in69SH = find(peaksSH_kp>6 & peaksSH_kp<=9);

%Intentisites for each group
maxI_NH03 = peaksNH_maxI(in03NH); alt_NH03 = peaksNH_alt(in03NH);
maxI_NH36 = peaksNH_maxI(in36NH); alt_NH36 = peaksNH_alt(in36NH);
maxI_NH69 = peaksNH_maxI(in69NH); alt_NH69 = peaksNH_alt(in69NH);
maxI_SH03 = peaksSH_maxI(in03SH); alt_SH03 = peaksSH_alt(in03SH);
maxI_SH36 = peaksSH_maxI(in36SH); alt_SH36 = peaksSH_alt(in36SH);
maxI_SH69 = peaksSH_maxI(in69SH); alt_SH69 = peaksSH_alt(in69SH);

%Find average and standard deviation for each group
avg_altSH_h_03 = zeros([24/seph 1]);
std_altSH_h_03 = zeros([24/seph 1]);
avg_altSH_h_36 = zeros([24/seph 1]);
std_altSH_h_36 = zeros([24/seph 1]);
avg_altSH_h_69 = zeros([24/seph 1]);
std_altSH_h_69 = zeros([24/seph 1]);
avg_altNH_h_03 = zeros([24/seph 1]);
std_altNH_h_03 = zeros([24/seph 1]);
avg_altNH_h_36 = zeros([24/seph 1]);
std_altNH_h_36 = zeros([24/seph 1]);
avg_altNH_h_69 = zeros([24/seph 1]);
std_altNH_h_69 = zeros([24/seph 1]);
hNH_03 = zeros([24/seph 1]);
hSH_03 = zeros([24/seph 1]);
hNH_36 = zeros([24/seph 1]);
hSH_36 = zeros([24/seph 1]);
hNH_69 = zeros([24/seph 1]);
hSH_69 = zeros([24/seph 1]);

avgalt_03 = 0.5*(sum(peaksSH_alt(in03SH))/length(in03SH)+sum(peaksNH_alt(in03NH))/length(in03NH));
avgalt_36 = 0.5*(sum(peaksSH_alt(in36SH))/length(in36SH)+sum(peaksNH_alt(in36NH))/length(in36NH));
avgalt_69 = 0.5*(sum(peaksSH_alt(in69SH))/length(in69SH)+sum(peaksNH_alt(in69NH))/length(in69NH));

for i = 1:24/seph
    %SH
    peaks_alt_kp = peaksSH_alt(in03SH);
    peaks_MLT_kp = peaksSH_MLT(in03SH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >2
        avg_altSH_h_03(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altSH_h_03(i)= std(peaks_alt_kp(ind_h))/sqrt(length(ind_h));
        hSH_03(i) = i-seph/2;
    end
    peaks_alt_kp = peaksSH_alt(in36SH);
    peaks_MLT_kp = peaksSH_MLT(in36SH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >2
        avg_altSH_h_36(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altSH_h_36(i)= std(peaks_alt_kp(ind_h))/sqrt(length(ind_h));
        hSH_36(i) = i-seph/2;
    end
    peaks_alt_kp = peaksSH_alt(in69SH);
    peaks_MLT_kp = peaksSH_MLT(in69SH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >2
        avg_altSH_h_69(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altSH_h_69(i)= std(peaks_alt_kp(ind_h))/sqrt(length(ind_h));
        hSH_69(i) = i-seph/2;
    end
    % NH
    peaks_alt_kp = peaksNH_alt(in03NH);
    peaks_MLT_kp = peaksNH_MLT(in03NH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >2
        avg_altNH_h_03(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altNH_h_03(i)= std(peaks_alt_kp(ind_h))/sqrt(length(ind_h));
        hNH_03(i) = i-seph/2;
    end
    peaks_alt_kp = peaksNH_alt(in36NH);
    peaks_MLT_kp = peaksNH_MLT(in36NH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >2
        avg_altNH_h_36(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altNH_h_36(i)= std(peaks_alt_kp(ind_h))/sqrt(length(ind_h));
        hNH_36(i) = i-seph/2;
    end
    peaks_alt_kp = peaksNH_alt(in69NH);
    peaks_MLT_kp = peaksNH_MLT(in69NH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >2
        avg_altNH_h_69(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altNH_h_69(i)= std(peaks_alt_kp(ind_h))/sqrt(length(ind_h));
        hNH_69(i) = i-seph/2;
    end
    length(ind_h)
end


hNH_03(hNH_03>14)=hNH_03(hNH_03>14)-4;
hNH_36(hNH_36>14)=hNH_36(hNH_36>14)-4;
hNH_69(hNH_69>14)=hNH_69(hNH_69>14)-4;
hSH_03(hSH_03>14)=hSH_03(hSH_03>14)-4;
hSH_36(hSH_36>14)=hSH_36(hSH_36>14)-4;
hSH_69(hSH_69>14)=hSH_69(hSH_69>14)-4;

fig = figure(position=[20 20 1600 550]); 
tiledlayout(1,5)
ax1 = nexttile([1 3]);hold on; grid; legend('NumColumns', 3, Location='southeast');
errorbar(nonzeros(hNH_03)-0.25 ,nonzeros(avg_altNH_h_03),nonzeros(std_altNH_h_03),'*', color = c1, DisplayName='NH kp=0-3')
errorbar(nonzeros(hNH_36)      ,nonzeros(avg_altNH_h_36),nonzeros(std_altNH_h_36),'*', color = c2, DisplayName='NH kp=3-6')
errorbar(nonzeros(hNH_69)+0.25 ,nonzeros(avg_altNH_h_69),nonzeros(std_altNH_h_69),'*', color = c3, DisplayName='NH kp=6-9')
errorbar(nonzeros(hSH_03)-0.25 ,nonzeros(avg_altSH_h_03),nonzeros(std_altSH_h_03),'^', color = c1, DisplayName='SH kp=0-3')
errorbar(nonzeros(hSH_36)      ,nonzeros(avg_altSH_h_36),nonzeros(std_altSH_h_36),'^', color = c2, DisplayName='SH kp=3-6')
errorbar(nonzeros(hSH_69)+0.25 ,nonzeros(avg_altSH_h_69),nonzeros(std_altSH_h_69),'^', color = c3, DisplayName='SH kp=6-9')
yline(avgalt_03, '-.',color = c1,DisplayName = append(num2str(avgalt_03,'%.1f'),' km'))
yline(avgalt_36, '-.',color = c2,DisplayName = append(num2str(avgalt_36,'%.1f'),' km'))
yline(avgalt_69, '-.',color = c3,DisplayName = append(num2str(avgalt_69,'%.1f'),' km'))
text(9.4,90,'//',fontsize=15)
ylabel('Altitude (km)')
xlabel('MLT')
xlim([0 20])
ylim([90 115])
xticks([0:1:9, 10:1:20]);
xticklabels([0:1:9, 14:1:24]);

ax2 = nexttile([1 2]);hold on; grid;
plot(maxI_NH03,alt_NH03,'*', color = c1,DisplayName='NH kp=0-3')
plot(maxI_NH36,alt_NH36,'*', color = c2,DisplayName='NH kp=3-6')
plot(maxI_NH69,alt_NH69,'*', color = c3,DisplayName='NH kp=6-9')
plot(maxI_SH03,alt_SH03,'^', color = c1,DisplayName='SH kp=0-3')
plot(maxI_SH36,alt_SH36,'^', color = c2,DisplayName='SH kp=3-6')
plot(maxI_SH69,alt_SH69,'^', color = c3,DisplayName='SH kp=6-9')
xlabel('Intensity (10^{13} photons \cdot nm^{-1} \cdot m^{-2} \cdot str^{-1} \cdot s^{-1})') ;
ylabel('Altitude (km)')

text(ax1,'Units', 'Normalized','Position', [0.96, 0.96],'string','a)')
text(ax2,'Units', 'Normalized','Position', [0.95, 0.96],'string','b)')

set(findall(gcf,'-property','FontSize'),'FontSize',14)
saveas(fig,savefig)



%% Feb/April separately

peaksNH_MLT =   [febpeaksNH.MLT];
peaksSH_MLT =   [febpeaksSH.MLT];
peaksNH_kp =    [febpeaksNH.kp ];
peaksSH_kp =    [febpeaksSH.kp ];
peaksNH_maxI =  [febpeaksNH.maxI];
peaksSH_maxI =  [febpeaksSH.maxI];
peaksNH_alt =   [febpeaksNH.alt];
peaksSH_alt =   [febpeaksSH.alt];

% Separate kp index in 3 groups for each hemisphere
in03NH = find(peaksNH_kp<=3);
in36NH = find(peaksNH_kp>3 & peaksNH_kp<=6);
in69NH = find(peaksNH_kp>6 & peaksNH_kp<=9);
in03SH = find(peaksSH_kp<=3);
in36SH = find(peaksSH_kp>3 & peaksSH_kp<=6);
in69SH = find(peaksSH_kp>6 & peaksSH_kp<=9);

%Intentisites for each group
maxI_NH03 = peaksNH_maxI(in03NH); alt_NH03 = peaksNH_alt(in03NH);
maxI_NH36 = peaksNH_maxI(in36NH); alt_NH36 = peaksNH_alt(in36NH);
maxI_NH69 = peaksNH_maxI(in69NH); alt_NH69 = peaksNH_alt(in69NH);
maxI_SH03 = peaksSH_maxI(in03SH); alt_SH03 = peaksSH_alt(in03SH);
maxI_SH36 = peaksSH_maxI(in36SH); alt_SH36 = peaksSH_alt(in36SH);
maxI_SH69 = peaksSH_maxI(in69SH); alt_SH69 = peaksSH_alt(in69SH);

%Find average and standard deviation for each group
avg_altSH_h_03 = zeros([24/seph 1]);
std_altSH_h_03 = zeros([24/seph 1]);
avg_altSH_h_36 = zeros([24/seph 1]);
std_altSH_h_36 = zeros([24/seph 1]);
avg_altSH_h_69 = zeros([24/seph 1]);
std_altSH_h_69 = zeros([24/seph 1]);
avg_altNH_h_03 = zeros([24/seph 1]);
std_altNH_h_03 = zeros([24/seph 1]);
avg_altNH_h_36 = zeros([24/seph 1]);
std_altNH_h_36 = zeros([24/seph 1]);
avg_altNH_h_69 = zeros([24/seph 1]);
std_altNH_h_69 = zeros([24/seph 1]);
hNH_03 = zeros([24/seph 1]);
hSH_03 = zeros([24/seph 1]);
hNH_36 = zeros([24/seph 1]);
hSH_36 = zeros([24/seph 1]);
hNH_69 = zeros([24/seph 1]);
hSH_69 = zeros([24/seph 1]);

avgalt_03 = 0.5*(sum(peaksSH_alt(in03SH))/length(in03SH)+sum(peaksNH_alt(in03NH))/length(in03NH));
avgalt_36 = 0.5*(sum(peaksSH_alt(in36SH))/length(in36SH)+sum(peaksNH_alt(in36NH))/length(in36NH));
avgalt_69 = 0.5*(sum(peaksSH_alt(in69SH))/length(in69SH)+sum(peaksNH_alt(in69NH))/length(in69NH));

for i = 1:24/seph
    %SH
    peaks_alt_kp = peaksSH_alt(in03SH);
    peaks_MLT_kp = peaksSH_MLT(in03SH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >2
        avg_altSH_h_03(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altSH_h_03(i)= std(peaks_alt_kp(ind_h))/sqrt(length(ind_h));
        hSH_03(i) = i-seph/2;
    end
    peaks_alt_kp = peaksSH_alt(in36SH);
    peaks_MLT_kp = peaksSH_MLT(in36SH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >2
        avg_altSH_h_36(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altSH_h_36(i)= std(peaks_alt_kp(ind_h))/sqrt(length(ind_h));
        hSH_36(i) = i-seph/2;
    end
    peaks_alt_kp = peaksSH_alt(in69SH);
    peaks_MLT_kp = peaksSH_MLT(in69SH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >2
        avg_altSH_h_69(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altSH_h_69(i)= std(peaks_alt_kp(ind_h))/sqrt(length(ind_h));
        hSH_69(i) = i-seph/2;
    end
    % NH
    peaks_alt_kp = peaksNH_alt(in03NH);
    peaks_MLT_kp = peaksNH_MLT(in03NH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >2
        avg_altNH_h_03(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altNH_h_03(i)= std(peaks_alt_kp(ind_h))/sqrt(length(ind_h));
        hNH_03(i) = i-seph/2;
    end
    peaks_alt_kp = peaksNH_alt(in36NH);
    peaks_MLT_kp = peaksNH_MLT(in36NH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >2
        avg_altNH_h_36(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altNH_h_36(i)= std(peaks_alt_kp(ind_h))/sqrt(length(ind_h));
        hNH_36(i) = i-seph/2;
    end
    peaks_alt_kp = peaksNH_alt(in69NH);
    peaks_MLT_kp = peaksNH_MLT(in69NH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >2
        avg_altNH_h_69(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altNH_h_69(i)= std(peaks_alt_kp(ind_h))/sqrt(length(ind_h));
        hNH_69(i) = i-seph/2;
    end
    length(ind_h)
end


hNH_03(hNH_03>14)=hNH_03(hNH_03>14)-5;
hNH_36(hNH_36>14)=hNH_36(hNH_36>14)-5;
hNH_69(hNH_69>14)=hNH_69(hNH_69>14)-5;
hSH_03(hSH_03>14)=hSH_03(hSH_03>14)-5;
hSH_36(hSH_36>14)=hSH_36(hSH_36>14)-5;
hSH_69(hSH_69>14)=hSH_69(hSH_69>14)-5;

fig = figure(position=[20 20 1600 550]); 
tiledlayout(1,5)
ax1 = nexttile([1 3]);hold on; grid; legend('NumColumns', 3, Location='southeast');
errorbar(nonzeros(hNH_03)-0.25 ,nonzeros(avg_altNH_h_03),nonzeros(std_altNH_h_03),'*', color = c1, DisplayName='NH kp=0-3')
errorbar(nonzeros(hNH_36)      ,nonzeros(avg_altNH_h_36),nonzeros(std_altNH_h_36),'*', color = c2, DisplayName='NH kp=3-6')
errorbar(nonzeros(hNH_69)+0.25 ,nonzeros(avg_altNH_h_69),nonzeros(std_altNH_h_69),'*', color = c3, DisplayName='NH kp=6-9')
errorbar(nonzeros(hSH_03)-0.25 ,nonzeros(avg_altSH_h_03),nonzeros(std_altSH_h_03),'^', color = c1, DisplayName='SH kp=0-3')
errorbar(nonzeros(hSH_36)      ,nonzeros(avg_altSH_h_36),nonzeros(std_altSH_h_36),'^', color = c2, DisplayName='SH kp=3-6')
errorbar(nonzeros(hSH_69)+0.25 ,nonzeros(avg_altSH_h_69),nonzeros(std_altSH_h_69),'^', color = c3, DisplayName='SH kp=6-9')
yline(avgalt_03, '-.',color = c1,DisplayName = append(num2str(avgalt_03,'%.1f'),' km'))
yline(avgalt_36, '-.',color = c2,DisplayName = append(num2str(avgalt_36,'%.1f'),' km'))
yline(avgalt_69, '-.',color = c3,DisplayName = append(num2str(avgalt_69,'%.1f'),' km'))
text(9.4,90,'//',fontsize=15)
ylabel('Altitude (km)')
xlabel('MLT')
xlim([0 20])
ylim([90 115])
xticks([0:1:9, 10:1:19]);
xticklabels([0:1:9, 14:1:24]);
title('February')

ax2 = nexttile([1 2]);hold on; grid;
plot(maxI_NH03,alt_NH03,'*', color = c1,DisplayName='NH kp=0-3')
plot(maxI_NH36,alt_NH36,'*', color = c2,DisplayName='NH kp=3-6')
plot(maxI_NH69,alt_NH69,'*', color = c3,DisplayName='NH kp=6-9')
plot(maxI_SH03,alt_SH03,'^', color = c1,DisplayName='SH kp=0-3')
plot(maxI_SH36,alt_SH36,'^', color = c2,DisplayName='SH kp=3-6')
plot(maxI_SH69,alt_SH69,'^', color = c3,DisplayName='SH kp=6-9')
xlabel('Intensity (10^{13} photons \cdot nm^{-1} \cdot m^{-2} \cdot str^{-1} \cdot s^{-1})') ;
ylabel('Altitude (km)')
ylim([90 115])
text(ax1,'Units', 'Normalized','Position', [0.96, 0.96],'string','a)')
text(ax2,'Units', 'Normalized','Position', [0.95, 0.96],'string','b)')

set(findall(gcf,'-property','FontSize'),'FontSize',14)
saveas(fig,savefigfeb)


%% April
peaksNH_MLT =   [aprpeaksNH.MLT ];
peaksSH_MLT =   [aprpeaksSH.MLT ];
peaksNH_kp =    [aprpeaksNH.kp  ];
peaksSH_kp =    [aprpeaksSH.kp  ];
peaksNH_maxI =  [aprpeaksNH.maxI];
peaksSH_maxI =  [aprpeaksSH.maxI];
peaksNH_alt =   [aprpeaksNH.alt ];
peaksSH_alt =   [aprpeaksSH.alt ];
in03NH = find(peaksNH_kp<=3);
in36NH = find(peaksNH_kp>3 & peaksNH_kp<=6);
in69NH = find(peaksNH_kp>6 & peaksNH_kp<=9);
in03SH = find(peaksSH_kp<=3);
in36SH = find(peaksSH_kp>3 & peaksSH_kp<=6);
in69SH = find(peaksSH_kp>6 & peaksSH_kp<=9);
maxI_NH03 = peaksNH_maxI(in03NH); alt_NH03 = peaksNH_alt(in03NH);
maxI_NH36 = peaksNH_maxI(in36NH); alt_NH36 = peaksNH_alt(in36NH);
maxI_NH69 = peaksNH_maxI(in69NH); alt_NH69 = peaksNH_alt(in69NH);
maxI_SH03 = peaksSH_maxI(in03SH); alt_SH03 = peaksSH_alt(in03SH);
maxI_SH36 = peaksSH_maxI(in36SH); alt_SH36 = peaksSH_alt(in36SH);
maxI_SH69 = peaksSH_maxI(in69SH); alt_SH69 = peaksSH_alt(in69SH);
avg_altSH_h_03 = zeros([24/seph 1]);
std_altSH_h_03 = zeros([24/seph 1]);
avg_altSH_h_36 = zeros([24/seph 1]);
std_altSH_h_36 = zeros([24/seph 1]);
avg_altSH_h_69 = zeros([24/seph 1]);
std_altSH_h_69 = zeros([24/seph 1]);
avg_altNH_h_03 = zeros([24/seph 1]);
std_altNH_h_03 = zeros([24/seph 1]);
avg_altNH_h_36 = zeros([24/seph 1]);
std_altNH_h_36 = zeros([24/seph 1]);
avg_altNH_h_69 = zeros([24/seph 1]);
std_altNH_h_69 = zeros([24/seph 1]);
hNH_03 = zeros([24/seph 1]);
hSH_03 = zeros([24/seph 1]);
hNH_36 = zeros([24/seph 1]);
hSH_36 = zeros([24/seph 1]);
hNH_69 = zeros([24/seph 1]);
hSH_69 = zeros([24/seph 1]);

avgalt_03 = 0.5*(sum(peaksSH_alt(in03SH))/length(in03SH)+sum(peaksNH_alt(in03NH))/length(in03NH));
avgalt_36 = 0.5*(sum(peaksSH_alt(in36SH))/length(in36SH)+sum(peaksNH_alt(in36NH))/length(in36NH));
avgalt_69 = 0.5*(sum(peaksSH_alt(in69SH))/length(in69SH)+sum(peaksNH_alt(in69NH))/length(in69NH));

for i = 1:24/seph
    %SH
    peaks_alt_kp = peaksSH_alt(in03SH);
    peaks_MLT_kp = peaksSH_MLT(in03SH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >2
        avg_altSH_h_03(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altSH_h_03(i)= std(peaks_alt_kp(ind_h))/sqrt(length(ind_h));
        hSH_03(i) = i-seph/2;
    end
    peaks_alt_kp = peaksSH_alt(in36SH);
    peaks_MLT_kp = peaksSH_MLT(in36SH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >2
        avg_altSH_h_36(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altSH_h_36(i)= std(peaks_alt_kp(ind_h))/sqrt(length(ind_h));
        hSH_36(i) = i-seph/2;
    end
    peaks_alt_kp = peaksSH_alt(in69SH);
    peaks_MLT_kp = peaksSH_MLT(in69SH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >2
        avg_altSH_h_69(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altSH_h_69(i)= std(peaks_alt_kp(ind_h))/sqrt(length(ind_h));
        hSH_69(i) = i-seph/2;
    end
    % NH
    peaks_alt_kp = peaksNH_alt(in03NH);
    peaks_MLT_kp = peaksNH_MLT(in03NH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >2
        avg_altNH_h_03(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altNH_h_03(i)= std(peaks_alt_kp(ind_h))/sqrt(length(ind_h));
        hNH_03(i) = i-seph/2;
    end
    peaks_alt_kp = peaksNH_alt(in36NH);
    peaks_MLT_kp = peaksNH_MLT(in36NH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >2
        avg_altNH_h_36(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altNH_h_36(i)= std(peaks_alt_kp(ind_h))/sqrt(length(ind_h));
        hNH_36(i) = i-seph/2;
    end
    peaks_alt_kp = peaksNH_alt(in69NH);
    peaks_MLT_kp = peaksNH_MLT(in69NH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >2
        avg_altNH_h_69(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altNH_h_69(i)= std(peaks_alt_kp(ind_h))/sqrt(length(ind_h));
        hNH_69(i) = i-seph/2;
    end
    length(ind_h)
end
hNH_03(hNH_03>14)=hNH_03(hNH_03>14)-5;
hNH_36(hNH_36>14)=hNH_36(hNH_36>14)-5;
hNH_69(hNH_69>14)=hNH_69(hNH_69>14)-5;
hSH_03(hSH_03>14)=hSH_03(hSH_03>14)-5;
hSH_36(hSH_36>14)=hSH_36(hSH_36>14)-5;
hSH_69(hSH_69>14)=hSH_69(hSH_69>14)-5;

fig = figure(position=[20 20 1600 550]); 
tiledlayout(1,5)
ax1 = nexttile([1 3]);hold on; grid; legend('NumColumns', 3, Location='southeast');
errorbar(nonzeros(hNH_03)-0.25 ,nonzeros(avg_altNH_h_03),nonzeros(std_altNH_h_03),'*', color = c1, DisplayName='NH kp=0-3')
errorbar(nonzeros(hNH_36)      ,nonzeros(avg_altNH_h_36),nonzeros(std_altNH_h_36),'*', color = c2, DisplayName='NH kp=3-6')
errorbar(nonzeros(hNH_69)+0.25 ,nonzeros(avg_altNH_h_69),nonzeros(std_altNH_h_69),'*', color = c3, DisplayName='NH kp=6-9')
errorbar(nonzeros(hSH_03)-0.25 ,nonzeros(avg_altSH_h_03),nonzeros(std_altSH_h_03),'^', color = c1, DisplayName='SH kp=0-3')
errorbar(nonzeros(hSH_36)      ,nonzeros(avg_altSH_h_36),nonzeros(std_altSH_h_36),'^', color = c2, DisplayName='SH kp=3-6')
errorbar(nonzeros(hSH_69)+0.25 ,nonzeros(avg_altSH_h_69),nonzeros(std_altSH_h_69),'^', color = c3, DisplayName='SH kp=6-9')
yline(avgalt_03, '-.',color = c1,DisplayName = append(num2str(avgalt_03,'%.1f'),' km'))
yline(avgalt_36, '-.',color = c2,DisplayName = append(num2str(avgalt_36,'%.1f'),' km'))
yline(avgalt_69, '-.',color = c3,DisplayName = append(num2str(avgalt_69,'%.1f'),' km'))
text(9.4,90,'//',fontsize=15)
ylabel('Altitude (km)')
xlabel('MLT')
xlim([0 20])
ylim([90 115])
xticks([0:1:9, 10:1:19]);
xticklabels([0:1:9, 14:1:24]);
title('April')

ax2 = nexttile([1 2]);hold on; grid;
plot(maxI_NH03,alt_NH03,'*', color = c1,DisplayName='NH kp=0-3')
plot(maxI_NH36,alt_NH36,'*', color = c2,DisplayName='NH kp=3-6')
plot(maxI_NH69,alt_NH69,'*', color = c3,DisplayName='NH kp=6-9')
plot(maxI_SH03,alt_SH03,'^', color = c1,DisplayName='SH kp=0-3')
plot(maxI_SH36,alt_SH36,'^', color = c2,DisplayName='SH kp=3-6')
plot(maxI_SH69,alt_SH69,'^', color = c3,DisplayName='SH kp=6-9')
xlabel('Intensity (10^{13} photons \cdot nm^{-1} \cdot m^{-2} \cdot str^{-1} \cdot s^{-1})') ;
ylabel('Altitude (km)')
ylim([90 115])

text(ax1,'Units', 'Normalized','Position', [0.96, 0.96],'string','a)')
text(ax2,'Units', 'Normalized','Position', [0.95, 0.96],'string','b)')

set(findall(gcf,'-property','FontSize'),'FontSize',14)
saveas(fig,savefigap)

%% March
peaksNH_MLT =   [marpeaksNH.MLT ];
peaksSH_MLT =   [marpeaksSH.MLT ];
peaksNH_kp =    [marpeaksNH.kp  ];
peaksSH_kp =    [marpeaksSH.kp  ];
peaksNH_maxI =  [marpeaksNH.maxI];
peaksSH_maxI =  [marpeaksSH.maxI];
peaksNH_alt =   [marpeaksNH.alt ];
peaksSH_alt =   [marpeaksSH.alt ];
in03NH = find(peaksNH_kp<=3);
in36NH = find(peaksNH_kp>3 & peaksNH_kp<=6);
in69NH = find(peaksNH_kp>6 & peaksNH_kp<=9);
in03SH = find(peaksSH_kp<=3);
in36SH = find(peaksSH_kp>3 & peaksSH_kp<=6);
in69SH = find(peaksSH_kp>6 & peaksSH_kp<=9);
maxI_NH03 = peaksNH_maxI(in03NH); alt_NH03 = peaksNH_alt(in03NH);
maxI_NH36 = peaksNH_maxI(in36NH); alt_NH36 = peaksNH_alt(in36NH);
maxI_NH69 = peaksNH_maxI(in69NH); alt_NH69 = peaksNH_alt(in69NH);
maxI_SH03 = peaksSH_maxI(in03SH); alt_SH03 = peaksSH_alt(in03SH);
maxI_SH36 = peaksSH_maxI(in36SH); alt_SH36 = peaksSH_alt(in36SH);
maxI_SH69 = peaksSH_maxI(in69SH); alt_SH69 = peaksSH_alt(in69SH);
avg_altSH_h_03 = zeros([24/seph 1]);
std_altSH_h_03 = zeros([24/seph 1]);
avg_altSH_h_36 = zeros([24/seph 1]);
std_altSH_h_36 = zeros([24/seph 1]);
avg_altSH_h_69 = zeros([24/seph 1]);
std_altSH_h_69 = zeros([24/seph 1]);
avg_altNH_h_03 = zeros([24/seph 1]);
std_altNH_h_03 = zeros([24/seph 1]);
avg_altNH_h_36 = zeros([24/seph 1]);
std_altNH_h_36 = zeros([24/seph 1]);
avg_altNH_h_69 = zeros([24/seph 1]);
std_altNH_h_69 = zeros([24/seph 1]);
hNH_03 = zeros([24/seph 1]);
hSH_03 = zeros([24/seph 1]);
hNH_36 = zeros([24/seph 1]);
hSH_36 = zeros([24/seph 1]);
hNH_69 = zeros([24/seph 1]);
hSH_69 = zeros([24/seph 1]);

avgalt_03 = 0.5*(sum(peaksSH_alt(in03SH))/length(in03SH)+sum(peaksNH_alt(in03NH))/length(in03NH));
avgalt_36 = 0.5*(sum(peaksSH_alt(in36SH))/length(in36SH)+sum(peaksNH_alt(in36NH))/length(in36NH));
avgalt_69 = 0.5*(sum(peaksSH_alt(in69SH))/length(in69SH)+sum(peaksNH_alt(in69NH))/length(in69NH));

for i = 1:24/seph
    %SH
    peaks_alt_kp = peaksSH_alt(in03SH);
    peaks_MLT_kp = peaksSH_MLT(in03SH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >2
        avg_altSH_h_03(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altSH_h_03(i)= std(peaks_alt_kp(ind_h))/sqrt(length(ind_h));
        hSH_03(i) = i-seph/2;
    end
    peaks_alt_kp = peaksSH_alt(in36SH);
    peaks_MLT_kp = peaksSH_MLT(in36SH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >2
        avg_altSH_h_36(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altSH_h_36(i)= std(peaks_alt_kp(ind_h))/sqrt(length(ind_h));
        hSH_36(i) = i-seph/2;
    end
    peaks_alt_kp = peaksSH_alt(in69SH);
    peaks_MLT_kp = peaksSH_MLT(in69SH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >2
        avg_altSH_h_69(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altSH_h_69(i)= std(peaks_alt_kp(ind_h))/sqrt(length(ind_h));
        hSH_69(i) = i-seph/2;
    end
    % NH
    peaks_alt_kp = peaksNH_alt(in03NH);
    peaks_MLT_kp = peaksNH_MLT(in03NH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >2
        avg_altNH_h_03(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altNH_h_03(i)= std(peaks_alt_kp(ind_h))/sqrt(length(ind_h));
        hNH_03(i) = i-seph/2;
    end
    peaks_alt_kp = peaksNH_alt(in36NH);
    peaks_MLT_kp = peaksNH_MLT(in36NH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >2
        avg_altNH_h_36(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altNH_h_36(i)= std(peaks_alt_kp(ind_h))/sqrt(length(ind_h));
        hNH_36(i) = i-seph/2;
    end
    peaks_alt_kp = peaksNH_alt(in69NH);
    peaks_MLT_kp = peaksNH_MLT(in69NH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >2
        avg_altNH_h_69(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altNH_h_69(i)= std(peaks_alt_kp(ind_h))/sqrt(length(ind_h));
        hNH_69(i) = i-seph/2;
    end
    length(ind_h)
end
hNH_03(hNH_03>14)=hNH_03(hNH_03>14)-5;
hNH_36(hNH_36>14)=hNH_36(hNH_36>14)-5;
hNH_69(hNH_69>14)=hNH_69(hNH_69>14)-5;
hSH_03(hSH_03>14)=hSH_03(hSH_03>14)-5;
hSH_36(hSH_36>14)=hSH_36(hSH_36>14)-5;
hSH_69(hSH_69>14)=hSH_69(hSH_69>14)-5;

fig = figure(position=[20 20 1600 550]); 
tiledlayout(1,5)
ax1 = nexttile([1 3]);hold on; grid; legend('NumColumns', 3, Location='southeast');
errorbar(nonzeros(hNH_03)-0.25 ,nonzeros(avg_altNH_h_03),nonzeros(std_altNH_h_03),'*', color = c1, DisplayName='NH kp=0-3')
errorbar(nonzeros(hNH_36)      ,nonzeros(avg_altNH_h_36),nonzeros(std_altNH_h_36),'*', color = c2, DisplayName='NH kp=3-6')
errorbar(nonzeros(hNH_69)+0.25 ,nonzeros(avg_altNH_h_69),nonzeros(std_altNH_h_69),'*', color = c3, DisplayName='NH kp=6-9')
errorbar(nonzeros(hSH_03)-0.25 ,nonzeros(avg_altSH_h_03),nonzeros(std_altSH_h_03),'^', color = c1, DisplayName='SH kp=0-3')
errorbar(nonzeros(hSH_36)      ,nonzeros(avg_altSH_h_36),nonzeros(std_altSH_h_36),'^', color = c2, DisplayName='SH kp=3-6')
errorbar(nonzeros(hSH_69)+0.25 ,nonzeros(avg_altSH_h_69),nonzeros(std_altSH_h_69),'^', color = c3, DisplayName='SH kp=6-9')
yline(avgalt_03, '-.',color = c1,DisplayName = append(num2str(avgalt_03,'%.1f'),' km'))
yline(avgalt_36, '-.',color = c2,DisplayName = append(num2str(avgalt_36,'%.1f'),' km'))
yline(avgalt_69, '-.',color = c3,DisplayName = append(num2str(avgalt_69,'%.1f'),' km'))
text(9.4,90,'//',fontsize=15)
ylabel('Altitude (km)')
xlabel('MLT')
xlim([0 20])
ylim([90 115])
xticks([0:1:9, 10:1:19]);
xticklabels([0:1:9, 14:1:24]);
title('March')

ax2 = nexttile([1 2]);hold on; grid;
plot(maxI_NH03,alt_NH03,'*', color = c1,DisplayName='NH kp=0-3')
plot(maxI_NH36,alt_NH36,'*', color = c2,DisplayName='NH kp=3-6')
plot(maxI_NH69,alt_NH69,'*', color = c3,DisplayName='NH kp=6-9')
plot(maxI_SH03,alt_SH03,'^', color = c1,DisplayName='SH kp=0-3')
plot(maxI_SH36,alt_SH36,'^', color = c2,DisplayName='SH kp=3-6')
plot(maxI_SH69,alt_SH69,'^', color = c3,DisplayName='SH kp=6-9')
xlabel('Intensity (10^{13} photons \cdot nm^{-1} \cdot m^{-2} \cdot str^{-1} \cdot s^{-1})') ;
ylabel('Altitude (km)')
ylim([90 115])

text(ax1,'Units', 'Normalized','Position', [0.96, 0.96],'string','a)')
text(ax2,'Units', 'Normalized','Position', [0.95, 0.96],'string','b)')

set(findall(gcf,'-property','FontSize'),'FontSize',14)
saveas(fig,savefigmar)



