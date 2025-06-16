%This code plots the IR2, IR3 and IR4 keograms with and without background
%Inclulding or not the downloading of the data

download = 1;

ti = datetime(2023,02,27,21,10,00);
te = datetime(2023,02,27,21,20,00);

filter = 'IR2';

%% Download nc file with images
if download
    getzarr = 'python \\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\get_zarr.py';
    changedirpath = 'C:\Nobackup\juditpcj\MATS\Images_J\All\';
    cd(changedirpath)
    str_ti = append('2023 ',num2str(month(ti)), ' ',num2str(day(ti)), ' ', num2str(hour(ti)), ' ', num2str(minute(ti)), ' ', num2str(second(ti)));
    str_te = append('2023 ',num2str(month(te)), ' ',num2str(day(te)), ' ', num2str(hour(te)), ' ', num2str(minute(te)), ' ', num2str(second(te)));
    call_python = append(getzarr,' -c ',filter,' -b ',str_ti,' -e ', str_te)
    system(call_python);
    cd \\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\


    %% Read nc file
    str_read_ti = append('2023_',num2str(month(ti)), '_',num2str(day(ti)), '_', num2str(hour(ti)), '_', num2str(minute(ti)), '_', num2str(second(ti)));
    str_read_te = append('2023_',num2str(month(te)), '_',num2str(day(te)), '_', num2str(hour(te)), '_', num2str(minute(te)), '_', num2str(second(te)));
    source = append('C:\Nobackup\juditpcj\MATS\Images_J\All\L1b_',filter,'-',str_read_ti,'-',str_read_te,'.nc');
    
    imcal = 'ImageCalibrated';
    imtime = 'time';
    imtplat = 'TPlat';
    imtplon = 'TPlon';
    imtpalt = 'TPheight';
    vardata = ncread(source,imcal);
    vartime = ncread(source,imtime);
    varlat = ncread(source,imtplat);
    varlon = ncread(source,imtplon);
    times = vartime*10^(-9);
    timepl = datetime(2000,01,01,00,00,00) + seconds(times);
end

for j = 1:length(vardata(1,1,:))
      mid_col = vardata(22,:,j)';
      keogram2(:,j) = mid_col';
end

keogram2_wob = zeros(size(keogram2));
for rown = 1:length(keogram2(:,1))
    row_raw = keogram2(rown,:);
    row = smoothdata(row_raw,'gaussian',4);
    tsecs = seconds(timepl-timepl(1));
    
    [fit1,~,mu] = polyfit(tsecs,row',3);
    line1 = polyval(fit1,tsecs,[],mu);
    row_wob = row-line1';
    keogram2_wob(rown,:) = row_wob;
end


pxi = 100;pyi = 100;pxf = 1400;pyf = 600;  
fig = figure(Position= [pxi pyi pxf pyf]); 
tiledlayout(2,1)
nexttile;
imagesc(timepl,1:length(vardata(1,:,1)),keogram2);
set(gca,'YDir','normal');
ylim([1 length(mid_col)])
title('IR2')
colorbar;
nexttile;
imagesc(timepl,1:length(vardata(1,:,1)),keogram2_wob,[0 max(keogram2_wob,[],'all')]);
set(gca,'YDir','normal');
ylim([1 length(mid_col)])
title('IR2 wob')
colorbar;

savefig = '\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub datafiles\MATS\IR2.png';
saveas(fig,savefig)



filter = 'IR3';

%% Download nc file with images
if download
    getzarr = 'python \\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\get_zarr.py';
    changedirpath = 'C:\Nobackup\juditpcj\MATS\Images_J\All\';
    cd(changedirpath)
    str_ti = append('2023 ',num2str(month(ti)), ' ',num2str(day(ti)), ' ', num2str(hour(ti)), ' ', num2str(minute(ti)), ' ', num2str(second(ti)));
    str_te = append('2023 ',num2str(month(te)), ' ',num2str(day(te)), ' ', num2str(hour(te)), ' ', num2str(minute(te)), ' ', num2str(second(te)));
    call_python = append(getzarr,' -c ',filter,' -b ',str_ti,' -e ', str_te)
    system(call_python);
    cd \\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\


    %% Read nc file
    str_read_ti = append('2023_',num2str(month(ti)), '_',num2str(day(ti)), '_', num2str(hour(ti)), '_', num2str(minute(ti)), '_', num2str(second(ti)));
    str_read_te = append('2023_',num2str(month(te)), '_',num2str(day(te)), '_', num2str(hour(te)), '_', num2str(minute(te)), '_', num2str(second(te)));
    source = append('C:\Nobackup\juditpcj\MATS\Images_J\All\L1b_',filter,'-',str_read_ti,'-',str_read_te,'.nc');
    
    imcal = 'ImageCalibrated';
    imtime = 'time';
    imtplat = 'TPlat';
    imtplon = 'TPlon';
    imtpalt = 'TPheight';
    vardata3 = ncread(source,imcal);
    vartime3 = ncread(source,imtime);
    varlat3 = ncread(source,imtplat);
    varlon3 = ncread(source,imtplon);
    times3 = vartime3*10^(-9);
    timepl3 = datetime(2000,01,01,00,00,00) + seconds(times3);
end

for j = 1:length(vardata3(1,1,:))
      mid_col = vardata3(5,:,j)';
      keogram3(:,j) = mid_col';
end

keogram3_wob = zeros(size(keogram3));
for rown = 1:length(keogram3(:,1))
    row_raw = keogram3(rown,:);
    row = smoothdata(row_raw,'gaussian',4);
    tsecs = seconds(timepl3-timepl3(1));
    
    [fit1,~,mu] = polyfit(tsecs,row',3);
    line1 = polyval(fit1,tsecs,[],mu);
    row_wob = row-line1';
    keogram3_wob(rown,:) = row_wob;
end

%%
pxi = 100;pyi = 100;pxf = 1400;pyf = 600;  
fig = figure(Position= [pxi pyi pxf pyf]); 
tiledlayout(2,1)
nexttile;
imagesc(timepl,1:length(vardata3(1,:,1)),keogram3,[0 2e14]);
set(gca,'YDir','normal');
ylim([1 length(mid_col)])
title('IR3')
colorbar;
nexttile;
imagesc(timepl,1:length(vardata3(1,:,1)),keogram3_wob,[0 0.1e14]);
set(gca,'YDir','normal');
ylim([1 length(mid_col)])
title('IR3 wob')
colorbar;
savefig = '\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub datafiles\MATS\IR3.png';
saveas(fig, savefig)

%%

filter = 'IR4';

%% Download nc file with images
if download
    getzarr = 'python \\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\get_zarr.py';
    changedirpath = 'C:\Nobackup\juditpcj\MATS\Images_J\All\';
    cd(changedirpath)
    str_ti = append('2023 ',num2str(month(ti)), ' ',num2str(day(ti)), ' ', num2str(hour(ti)), ' ', num2str(minute(ti)), ' ', num2str(second(ti)));
    str_te = append('2023 ',num2str(month(te)), ' ',num2str(day(te)), ' ', num2str(hour(te)), ' ', num2str(minute(te)), ' ', num2str(second(te)));
    call_python = append(getzarr,' -c ',filter,' -b ',str_ti,' -e ', str_te)
    system(call_python);
    cd \\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\


    %% Read nc file
    str_read_ti = append('2023_',num2str(month(ti)), '_',num2str(day(ti)), '_', num2str(hour(ti)), '_', num2str(minute(ti)), '_', num2str(second(ti)));
    str_read_te = append('2023_',num2str(month(te)), '_',num2str(day(te)), '_', num2str(hour(te)), '_', num2str(minute(te)), '_', num2str(second(te)));
    source = append('C:\Nobackup\juditpcj\MATS\Images_J\All\L1b_',filter,'-',str_read_ti,'-',str_read_te,'.nc');
    
    vardata4 = ncread(source,imcal);
    vartime4 = ncread(source,imtime);
    varlat4 = ncread(source,imtplat);
    varlon4 = ncread(source,imtplon);
    times4 = vartime4*10^(-9);
    timepl4 = datetime(2000,01,01,00,00,00) + seconds(times4);
end

%% Make keogram
for j = 1:length(vardata4(1,1,:))
      mid_col = vardata4(5,:,j)';
      keogram4(:,j) = mid_col';
end

keogram4_wob = zeros(size(keogram4));
for rown = 1:length(keogram4(:,1))
    row_raw = keogram4(rown,:);
    row = smoothdata(row_raw,'gaussian',4);
    tsecs = seconds(timepl4-timepl4(1));
    
    [fit1,~,mu] = polyfit(tsecs,row',3);
    line1 = polyval(fit1,tsecs,[],mu);
    row_wob = row-line1';
    keogram4_wob(rown,:) = row_wob;
end


pxi = 100;pyi = 100;pxf = 1400;pyf = 600;  
fig = figure(Position= [pxi pyi pxf pyf]); 
tiledlayout(2,1)
nexttile;
imagesc(timepl,1:length(vardata4(1,:,1)),keogram4);
set(gca,'YDir','normal');
ylim([1 length(mid_col)])
title('IR4')
colorbar;
nexttile;
imagesc(timepl,1:length(vardata4(1,:,1)),keogram4_wob,[0 max(keogram4_wob,[],'all')]);
set(gca,'YDir','normal');
ylim([1 length(mid_col)])
title('IR4 wob')
colorbar;

savefig = '\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub datafiles\MATS\IR4.png';
saveas(fig, savefig)

