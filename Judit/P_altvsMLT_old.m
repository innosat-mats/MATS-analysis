%% 

clear workspace
clc

savename ='\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub datafiles\MATS\alt_MLT.png';
savefig = '\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub datafiles\MATS\Fig03.png';
seph = 1;

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


%% Figure 1 - scatter color plot
fig = figure(position=[20 20 1200 600]); hold on; grid; legend();
indbefmidnight = find(peaksNH_MLT<24 & peaksNH_MLT>12);
indaftmidnight = find(peaksNH_MLT>0 & peaksNH_MLT<12);
scatter(peaksNH_MLT(indbefmidnight),peaksNH_alt(indbefmidnight),[],peaksNH_kp(indbefmidnight),'*',DisplayName = 'Northern Hemisphere')
scatter(peaksNH_MLT(indaftmidnight)+24,peaksNH_alt(indaftmidnight),[],peaksNH_kp(indaftmidnight),'*',HandleVisibility='off')
indbefmidnight = find(peaksSH_MLT<24 & peaksSH_MLT>12);
indaftmidnight = find(peaksSH_MLT>0  & peaksSH_MLT<12);
scatter(peaksSH_MLT(indbefmidnight),peaksSH_alt(indbefmidnight),[],peaksSH_kp(indbefmidnight),'o',DisplayName = 'Southern Hemisphere')
scatter(peaksSH_MLT(indaftmidnight)+24,peaksSH_alt(indaftmidnight),[],peaksSH_kp(indaftmidnight),'o',HandleVisibility='off')
xticks([12:2:22, 24:2:36]);
xticklabels([12:2:22, 0:2:12]);
xlim([12 36]);
set(gca,"CLim",[0 9])
cb = colorbar ;
colormap jet
cb.Label.String = 'Kp-index' ;
ylabel('Altitude (km)')
xlabel('MLT')
saveas(fig,savename)


%% Figure 2 - error bars all together

avg_alt_h = zeros([24/seph 1]);
std_alt_h = zeros([24/seph 1]);
h = zeros([24/seph 1]);
peaks_MLT = [peaksSH_MLT peaksNH_MLT];
peaks_alt = [peaksSH_alt peaksNH_alt];

for i = 1:24/seph
    ind_h = find(peaks_MLT >= i-seph & peaks_MLT < i);
    avg_alt_h(i) = sum(peaks_alt(ind_h))/length(ind_h);
    std_alt_h(i) = std(peaks_alt(ind_h));
    h(i) = i-seph/2;
    length(ind_h)
end

fig = figure(position=[20 20 1200 600]); hold on; grid;
indbefmidnight = find(h<24 & h>12);
indaftmidnight = find(h>00 & h<12);
errorbar(h(indbefmidnight),avg_alt_h(indbefmidnight),std_alt_h(indbefmidnight),'b*')
errorbar(h(indaftmidnight)+24,avg_alt_h(indaftmidnight),std_alt_h(indaftmidnight),'b*')
xticks([12:2:22, 24:2:36]);
xticklabels([12:2:22, 0:2:12]);
xlim([12 36]);
ylabel('Altitude (km)')
xlabel('MLT')
%saveas(fig,savename)

%% Figure 3 - With error bars NH / SH

avg_altSH_h = zeros([24/seph 1]);
std_altSH_h = zeros([24/seph 1]);
avg_altNH_h = zeros([24/seph 1]);
std_altNH_h = zeros([24/seph 1]);
hNH = zeros([24/seph 1]);
hSH = zeros([24/seph 1]);

for i = 1:24/seph
    ind_h = find(peaksSH_MLT >= i-seph & peaksSH_MLT < i);
    if length(ind_h) >4
        avg_altSH_h(i) = sum(peaksSH_alt(ind_h))/length(ind_h);
        std_altSH_h(i) = std(peaksSH_alt(ind_h));
        hSH(i) = i-seph/2;
    end
    ind_h = find(peaksNH_MLT >= i-seph & peaksNH_MLT < i);
    if length(ind_h) >4
        avg_altNH_h(i) = sum(peaksNH_alt(ind_h))/length(ind_h);
        std_altNH_h(i) = std(peaksNH_alt(ind_h));
        hNH(i) = i-seph/2;
    end
    length(ind_h)
end

fig = figure(position=[20 20 1200 600]); hold on; grid; legend;
indbefmidnight = find(hNH<24 & hNH>12);
indaftmidnight = find(hNH>00 & hNH<12);

h = nonzeros(hNH(indbefmidnight));
avg_alt = nonzeros(avg_altNH_h(indbefmidnight));
std_alt = nonzeros(std_altNH_h(indbefmidnight));
errorbar(h   ,avg_alt,std_alt,'b*', DisplayName='Northern Hemisphere')

h = nonzeros(hNH(indaftmidnight));
avg_alt = nonzeros(avg_altNH_h(indaftmidnight));
std_alt = nonzeros(std_altNH_h(indaftmidnight));

indbefmidnight = find(hSH<24 & hSH>12);
indaftmidnight = find(hSH>00 & hSH<12);
errorbar(h+24,avg_alt,std_alt,'b*',HandleVisibility='off')
h = nonzeros(hSH(indbefmidnight));
avg_alt = nonzeros(avg_altSH_h(indbefmidnight));
std_alt = nonzeros(std_altSH_h(indbefmidnight));
errorbar(h   ,avg_alt,std_alt,'r*', DisplayName='Southern Hemisphere')

h = nonzeros(hSH(indaftmidnight));
avg_alt = nonzeros(avg_altSH_h(indaftmidnight));
std_alt = nonzeros(std_altSH_h(indaftmidnight));
errorbar(h+24,avg_alt,std_alt,'r*',HandleVisibility='off')

xticks([12:2:22, 24:2:36]);
xticklabels([12:2:22, 0:2:12]);
xlim([12 36]);
ylabel('Altitude (km)')
xlabel('MLT')
%saveas(fig,savename)

%% Figure 4 - With error bars separating by kp index

% Separate kp index in 3 groups for each hemisphere
in03NH = find(peaksNH_kp<=3);
in36NH = find(peaksNH_kp>3 & peaksNH_kp<=6);
in69NH = find(peaksNH_kp>6 & peaksNH_kp<=9);
in03SH = find(peaksSH_kp<=3);
in36SH = find(peaksSH_kp>3 & peaksSH_kp<=6);
in69SH = find(peaksSH_kp>6 & peaksSH_kp<=9);

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

for i = 1:24/seph
    %SH
    peaks_alt_kp = peaksSH_alt(in03SH);
    peaks_MLT_kp = peaksSH_MLT(in03SH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >3
        avg_altSH_h_03(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altSH_h_03(i)= std(peaks_alt_kp(ind_h));
        hSH_03(i) = i-seph/2;
    end
    peaks_alt_kp = peaksSH_alt(in36SH);
    peaks_MLT_kp = peaksSH_MLT(in36SH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >3
        avg_altSH_h_36(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altSH_h_36(i)= std(peaks_alt_kp(ind_h));
        hSH_36(i) = i-seph/2;
    end
    peaks_alt_kp = peaksSH_alt(in69SH);
    peaks_MLT_kp = peaksSH_MLT(in69SH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >3
        avg_altSH_h_69(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altSH_h_69(i)= std(peaks_alt_kp(ind_h));
        hSH_69(i) = i-seph/2;
    end
    % NH
    peaks_alt_kp = peaksNH_alt(in03NH);
    peaks_MLT_kp = peaksNH_MLT(in03NH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >3
        avg_altNH_h_03(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altNH_h_03(i)= std(peaks_alt_kp(ind_h));
        hNH_03(i) = i-seph/2;
    end
    peaks_alt_kp = peaksNH_alt(in36NH);
    peaks_MLT_kp = peaksNH_MLT(in36NH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >3
        avg_altNH_h_36(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altNH_h_36(i)= std(peaks_alt_kp(ind_h));
        hNH_36(i) = i-seph/2;
    end
    peaks_alt_kp = peaksNH_alt(in69NH);
    peaks_MLT_kp = peaksNH_MLT(in69NH);
    ind_h = find(peaks_MLT_kp >= i-seph & peaks_MLT_kp < i);
    if length(ind_h) >3
        avg_altNH_h_69(i)= sum(peaks_alt_kp(ind_h))/length(ind_h);
        std_altNH_h_69(i)= std(peaks_alt_kp(ind_h));
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

fig = figure(position=[20 20 1500 500]); hold on; grid; legend(Location='southeast');
errorbar(nonzeros(hNH_03)-0.25 ,nonzeros(avg_altNH_h_03),nonzeros(std_altNH_h_03),'b*', DisplayName='NH kp=0-3')
errorbar(nonzeros(hNH_36) ,nonzeros(avg_altNH_h_36),nonzeros(std_altNH_h_36),'g*', DisplayName='NH kp=3-6')
errorbar(nonzeros(hNH_69)+0.25 ,nonzeros(avg_altNH_h_69),nonzeros(std_altNH_h_69),'r*', DisplayName='NH kp=6-9')
errorbar(nonzeros(hSH_03)-0.25 ,nonzeros(avg_altSH_h_03),nonzeros(std_altSH_h_03),'b^', DisplayName='SH kp=0-3')
errorbar(nonzeros(hSH_36) ,nonzeros(avg_altSH_h_36),nonzeros(std_altSH_h_36),'g^', DisplayName='SH kp=3-6')
errorbar(nonzeros(hSH_69)+0.25 ,nonzeros(avg_altSH_h_69),nonzeros(std_altSH_h_69),'r^', DisplayName='SH kp=6-9')
text(9.4,94,'//',fontsize=15)
ylabel('Altitude (km)')
xlabel('MLT')
xlim([0 20])
xticks([0:1:9, 10:1:19]);
xticklabels([0:1:9, 14:1:24]);

set(findall(gcf,'-property','FontSize'),'FontSize',16)
saveas(fig,savefig)




