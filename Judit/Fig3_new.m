%% 

clear workspace
clc

savename ='\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub datafiles\MATS\alt_MLT.png';
savefig = '\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub datafiles\MATS\Fig03.png';
seph = 1;

c1= [0.0 0.3 0.9];
c2= [0.0 0.8 0.4];
c3= [0.9 0.2 0.5];

newintensities = 1;

sigma_I = 0.1;
sigma_h = 1;

lambda = sigma_I/sigma_h;

%Files to plot the peak points
if newintensities
    addpath("\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structuresnew\IR2\")
    load("febpeaksNH.mat");
    load("febpeaksSH.mat") ;
    load("marpeaksNH.mat");
    load("marpeaksSH.mat") ;
    load("aprpeaksNH.mat");
    load("aprpeaksSH.mat");
    peaksNH_MLT =  [febNH.MLT  ,marNH.MLT  ,aprNH.MLT];
    peaksSH_MLT =  [febSH.MLT  ,marSH.MLT  ,aprSH.MLT];
    peaksNH_kp =   [febNH.kp   ,marNH.kp   ,aprNH.kp];
    peaksSH_kp =   [febSH.kp   ,marSH.kp   ,aprSH.kp];
    peaksNH_maxI = [febNH.newI ,marNH.newI ,aprNH.newI]*3; %The 3 is fixing an indexing problem for the algorithm
    peaksSH_maxI = [febSH.newI ,marSH.newI ,aprSH.newI]*3; %The 3 is fixing an indexing problem for the algorithm
    peaksNH_alt =  [febNH.altnew,marNH.altnew  ,aprNH.altnew]/1000;
    peaksSH_alt =  [febSH.altnew,marSH.altnew  ,aprSH.altnew]/1000;
    peaksNH_tim =  [febNH.time,marNH.time  ,aprNH.time];
    peaksSH_tim =  [febSH.time,marSH.time  ,aprSH.time];
else
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
end

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

fig = figure(position=[20 20 1650 500]); 
tiledlayout(1,3,'TileSpacing','tight')
ax1 = nexttile([1 2]);hold on; grid; legend('NumColumns', 3, Location='southeast');
errorbar(nonzeros(hNH_03)-0.25 ,nonzeros(avg_altNH_h_03),nonzeros(std_altNH_h_03),'*', color = c1, DisplayName='NH kp=0-3')
errorbar(nonzeros(hNH_36)      ,nonzeros(avg_altNH_h_36),nonzeros(std_altNH_h_36),'*', color = c2, DisplayName='NH kp=3-6')
errorbar(nonzeros(hNH_69)+0.25 ,nonzeros(avg_altNH_h_69),nonzeros(std_altNH_h_69),'*', color = c3, DisplayName='NH kp=6-9')
errorbar(nonzeros(hSH_03)-0.25 ,nonzeros(avg_altSH_h_03),nonzeros(std_altSH_h_03),'^', color = c1, DisplayName='SH kp=0-3')
errorbar(nonzeros(hSH_36)      ,nonzeros(avg_altSH_h_36),nonzeros(std_altSH_h_36),'^', color = c2, DisplayName='SH kp=3-6')
errorbar(nonzeros(hSH_69)+0.25 ,nonzeros(avg_altSH_h_69),nonzeros(std_altSH_h_69),'^', color = c3, DisplayName='SH kp=6-9')
yline(avgalt_03, '-.',color = c1,DisplayName = append(num2str(avgalt_03,'%.1f'),' km'))
yline(avgalt_36, '-.',color = c2,DisplayName = append(num2str(avgalt_36,'%.1f'),' km'))
yline(avgalt_69, '-.',color = c3,DisplayName = append(num2str(avgalt_69,'%.1f'),' km'))
text(9.4,95,'//',fontsize=15)
ylabel('Altitude (km)')
xlabel('MLT')
xlim([0 20])
ylim([95 110])
xticks([0:1:9, 10:1:20]);
xticklabels([0:1:9, 14:1:24]);




fit_INH = [maxI_NH03,maxI_NH36,maxI_NH69,maxI_SH03,maxI_SH36,maxI_SH69];
fit_hNH = [alt_NH03 ,alt_NH36 ,alt_NH69,alt_SH03 ,alt_SH36 ,alt_SH69 ];
[b, sigma2_x, hfitNH, IfitNH, stats]  = deming(fit_hNH',fit_INH',lambda);
display(sigma2_x); display(stats)
%hfitNH = linspace(min(fit_hNH), max(fit_hNH), 100);
%IfitNH = b(1)+ b(2) * hfitNH;
mNH = b(2);


ax2 = nexttile; hold on; grid;
plot(alt_NH03(maxI_NH03>0),nonzeros(maxI_NH03)','*', color = c1,HandleVisibility = 'off')
plot(alt_NH36(maxI_NH36>0),nonzeros(maxI_NH36)','*', color = c2,HandleVisibility = 'off')
plot(alt_NH69(maxI_NH69>0),nonzeros(maxI_NH69)','*', color = c3,HandleVisibility = 'off')
plot(alt_SH03(maxI_SH03>0),nonzeros(maxI_SH03)','^', color = c1,HandleVisibility = 'off')
plot(alt_SH36(maxI_SH36>0),nonzeros(maxI_SH36)','^', color = c2,HandleVisibility = 'off')
plot(alt_SH69(maxI_SH69>0),nonzeros(maxI_SH69)','^', color = c3,HandleVisibility = 'off')
plot(hfitNH,IfitNH,'k-',DisplayName=append('m = ',num2str(mNH,'%.2f')))

ylabel('ph \cdot nm^{-1} \cdot m^{-2} \cdot sr^{-1} \cdot s^{-1}') ;
xlabel('Altitude (km)')
xlim([90 115])
ylim([0 5.25e14])

text(ax1,'Units', 'Normalized','Position', [0.97, 0.96],'string','a)')
text(ax2,'Units', 'Normalized','Position', [0.95, 0.96],'string','b)')

set(findall(gcf,'-property','FontSize'),'FontSize',14)
saveas(fig,savefig)





function [xfit, yfit, m_sym] = fitsimetric(x, y)
    % Ajust "clàssic": y = m1*x + b1
    p_xy = polyfit(x, y, 1);
    m1 = p_xy(1);
    b1 = p_xy(2);
    
    % Ajust invers: x = m2*y + b2 → invertim després per obtenir y = ...
    p_yx = polyfit(y, x, 1);
    m2 = p_yx(1);
    b2 = p_yx(2);
    
    % Pendent simètrica (mitjana geomètrica dels pendents directes i inversos)
    m_sym = sqrt(m1 / m2);  % Com que m2 = dx/dy, 1/m2 ≈ dy/dx
    % Opcionalment pots agafar signe de m1:
    if m1 < 0
        m_sym = -abs(m_sym);
    else
        m_sym = abs(m_sym);
    end
    
    % Per trobar la intersecció amb m_sym, fem servir el centre de masses
    x0 = mean(x);
    y0 = mean(y);
    
    % Recta: y = y0 + m_sym * (x - x0)
    xfit = linspace(min(x), max(x), 100);
    yfit = y0 + m_sym * (xfit - x0);
end
  