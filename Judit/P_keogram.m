% Download data and Create a keogram from date and filter

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
end

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

%% Make keogram
for j = 1:length(vardata(1,1,:))
      mid_col = vardata(22,:,j)';
      keogram(:,j) = mid_col';
end

pxi = 100;pyi = 100;pxf = 1700;pyf = 600;  
figure(Position= [pxi pyi pxf pyf]); 
imagesc(timepl,1:length(vardata(1,:,1)),keogram,[0 6e14]);
set(gca,'YDir','normal');
ylim([1 length(mid_col)])
title(filter)
colorbar;
%xlim([timepl(1) timepl(end)])
%xticks(timepl(1):minutes(10):timepl(end))

ticklabelstim = timepl(1):minutes(5):timepl(end);
ticklabelslat = varlat(ismember(timepl,ticklabelstim));
labelArray = [compose('%5s',datestr(ticklabelstim,'HH:MM'))';compose('%.fº',ticklabelslat')];
labelArray = strjust(pad(labelArray),'center');
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
ax = gca(); 
ax.XTick = timepl(1):minutes(5):timepl(end); 
ax.XLim = [timepl(1) timepl(end)];
ax.XTickLabel = tickLabels; 

