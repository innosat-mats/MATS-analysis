
addpath("\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structuresnew\IR2\")
load("febpeaksNH.mat");
load("febpeaksSH.mat") ;
load("marpeaksNH.mat");
load("marpeaksSH.mat") ;
load("aprpeaksNH.mat");
load("aprpeaksSH.mat");

times_list = [febNH.time,febSH.time,marNH.time,marSH.time,aprNH.time,aprSH.time];
rows_list =  [febNH.row ,febSH.row ,marNH.row ,marSH.row ,aprNH.row ,aprSH.row ];
newI_list =  [febNH.newI,febSH.newI,marNH.newI,marSH.newI,aprNH.newI,aprSH.newI];
MLT_list =   [febNH.MLT ,febSH.MLT ,marNH.MLT ,marSH.MLT ,aprNH.MLT ,aprSH.MLT ];


alt_list = [];
for i = 1:length(times_list)
    getzarr = 'python \\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\get_zarr.py';
    changedirpath = 'C:\Nobackup\juditpcj\MATS\Images_alt\';
    cd(changedirpath)
    disp(times_list(i));
    ti = times_list(i)-seconds(4);
    te = times_list(i)+seconds(4);
    str_ti = append('2023 ',num2str(month(ti)), ' ',num2str(day(ti)), ' ', num2str(hour(ti)), ' ', num2str(minute(ti)), ' ', num2str(second(ti)));
    str_te = append('2023 ',num2str(month(te)), ' ',num2str(day(te)), ' ', num2str(hour(te)), ' ', num2str(minute(te)), ' ', num2str(second(te)));
    
    call_python = append(getzarr,' -c IR2 -b ',str_ti,' -e ', str_te, '  --tp_coords')
    system(call_python);
    cd \\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\
    
    source = append('C:\Nobackup\juditpcj\MATS\Images_alt\L1b_IR2-2023_',num2str(month(ti)),'_',num2str(day(ti)), '_', num2str(hour(ti)), '_', num2str(minute(ti)), '_', num2str(second(ti)),...
        '-2023_',num2str(month(te)),'_',num2str(day(te)), '_', num2str(hour(te)), '_', num2str(minute(te)), '_', num2str(second(te)),'.nc');
    altpix = ncread(source,'TPheightPixel');
    alt = altpix(22,rows_list(i));

    alt_list(i) = alt;
end

lfn =length(febNH.time);
lfs = length(febSH.time);
lmn = length(marNH.time);
lms = length(marSH.time);
lan = length(aprNH.time);
las = length(aprSH.time);

febNH.altnew = alt_list(1:lfn)
febSH.altnew = alt_list(lfn+1:lfn+lfs)
marNH.altnew = alt_list(lfn+lfs+1:lfn+lfs+lmn)
marSH.altnew = alt_list(lfn+lfs+lmn+1:lfn+lfs+lmn+lms)
aprNH.altnew = alt_list(lfn+lfs+lmn+lms+1:lfn+lfs+lmn+lms+lan)
aprSH.altnew = alt_list(lfn+lfs+lmn+lms+lan+1:end)
%
save('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structuresnew\IR2\febpeaksNH.mat','febNH')
save('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structuresnew\IR2\febpeaksSH.mat','febSH')
save('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structuresnew\IR2\marpeaksNH.mat','marNH')
save('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structuresnew\IR2\marpeaksSH.mat','marSH')
save('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structuresnew\IR2\aprpeaksNH.mat','aprNH')
save('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structuresnew\IR2\aprpeaksSH.mat','aprSH')