

altitude = 90:1:120;
latitude = -70;
longitude = 10;

savefig = append('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub datafiles\MATS\Figmodel',num2str(latitude),'.png');

savefig2 = append('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub datafiles\MATS\Figmodel',num2str(latitude),'r.png');
t_00 = datetime(2023,01,01,00,00,00);
t_f1 = datetime(2023,02,15,06,00,00);
t_f2 = datetime(2023,02,15,18,00,00);
t_m1 = datetime(2023,03,15,06,00,00);
t_m2 = datetime(2023,03,15,18,00,00);
t_a1 = datetime(2023,04,15,06,00,00);
t_a2 = datetime(2023,04,15,18,00,00);

t = [t_f1, t_f2, t_m1, t_m2, t_a1, t_a2];

year = 2023;


fig=figure(Position=[10 10 1500 800]); hold on; grid; 
aa = tiledlayout(1,4, 'tilespacing','tight');
ax1= nexttile; hold on; grid;subtitle('T'); legend(Location='southeast')
ax2= nexttile; hold on; grid;subtitle('O density')  ;legend();yticks([]);%xscale('log')
ax3= nexttile; hold on; grid;subtitle('N_2 density');legend();yticks([]);%xscale('log')
ax4= nexttile; hold on; grid;subtitle('O_2 density');legend();yticks([]);%xscale('log')
title(aa,append('Latitude: ',num2str(latitude),' ; Longitude: ', num2str(longitude)))

for i = 1:length(t)
    dayOfYear = floor(days(t(i)-t_00));
    UTseconds = (days(t(i)-t_00)-dayOfYear)*24*60*60;
    [T, rho] = atmosnrlmsise00(altitude*1000,latitude,longitude,year,dayOfYear,UTseconds);
    plot(ax1,T(:,2),altitude,'DisplayName',datestr(t(i)))
    plot(ax2,rho(:,2),altitude,'DisplayName',datestr(t(i)))
    plot(ax3,rho(:,3),altitude,'DisplayName',datestr(t(i)))
    plot(ax4,rho(:,4),altitude,'DisplayName',datestr(t(i)))
end

saveas(fig,savefig)


%%

fig=figure(Position=[10 10 1200 800]); hold on; grid; 
aa = tiledlayout(1,3, 'tilespacing','tight');
ax1= nexttile; hold on; grid;subtitle('T'); legend(Location='southeast')
ax2= nexttile; hold on; grid;subtitle('O density')  ;legend();yticks([]);%xscale('log')
ax3= nexttile; hold on; grid;subtitle('O_2/N_2 ');legend();yticks([]);xscale('log')
title(aa,append('Latitude: ',num2str(latitude),' ; Longitude: ', num2str(longitude)))

for i = 1:length(t)
    dayOfYear = floor(days(t(i)-t_00));
    UTseconds = (days(t(i)-t_00)-dayOfYear)*24*60*60;
    [T, rho] = atmosnrlmsise00(altitude*1000,latitude,longitude,year,dayOfYear,UTseconds);
    plot(ax1,T(:,2),altitude,'DisplayName',datestr(t(i)))
    plot(ax2,rho(:,2),altitude,'DisplayName',datestr(t(i)))
    plot(ax3,rho(:,4)./rho(:,3),altitude,'DisplayName',datestr(t(i)))
end

saveas(fig,savefig2)