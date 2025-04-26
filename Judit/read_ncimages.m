


source = '\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\L1b-2023_4_23_11_28_0-2023_4_23_11_52_0.nc';
%source = '\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\L1b_IR1-2023_4_23_0_0_0-2023_4_24_0_0_0.nc';
%source = '\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\February\L1b_IR1-2023_2_10_2_14_30-2023_2_10_2_22_30.nc';

%source = '\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\February\L1b_IR1-2023_2_9_23_0_44-2023_2_9_23_8_44.nc';
imcal = 'ImageCalibrated';

%vardata = ncinfo(source);
vardata = ncread(source,imcal);
figure;
imagesc(vardata(:,:,40)')

%%
%source = '\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\L1b-2023_3_2_19_16_50-2023_3_2_19_17_0.nc';
%source = '\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\L1b_IR1-2023_3_2_19_16_50-2023_3_2_19_17_0.nc';
%source = '\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\L1b-2023_4_23_11_28_0-2023_4_23_11_52_0.nc';
source = '\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\L1b_IR1-2023_4_23_0_0_0-2023_4_24_0_0_0.nc';
%source = '\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\L1b_IR1-2023_4_23_11_28_0-2023_4_23_11_52_0.nc';


imcal = 'ImageCalibrated';
imjpg = 'JPEGQ';
imtime = 'time';
i = 20
%varinfo = ncinfo(source);
vardata = ncread(source,imcal);
vartime = ncread(source,imtime);
    times = vartime*10^(-9);
    time = datetime(2000,01,01,00,00,00) + seconds(times);
image = vardata(:,end:-1:1,i);
%keogram_i=keogram_i(end:-1:1,:);
figure;
imagesc(image')
title(datestr(time(i)))
ncdisp(source)

%%
source = '\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\L1b_IR1-2023_4_23_11_28_0-2023_4_23_11_52_0.nc';


imcal = 'ImageCalibrated';
imtime = 'time';

%varinfo = ncinfo(source);
vardata = ncread(source,imcal);
vartime = ncread(source,imtime);
times = vartime*10^(-9);
time = datetime(2000,01,01,00,00,00) + seconds(times);
image = vardata(:,end:-1:1,1);
%keogram_i=keogram_i(end:-1:1,:);
figure;
imagesc(image')
title(datestr(time(1)))
ncdisp(source)