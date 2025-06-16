% This code plots the intensity gotten from IR1 vs the one gotten from IR2
close all; clear all;
retrieve = 0;
just_update = 1;
already_saved=0;

month_chosen = 'April';
hemisphere = 'S';
filter = 'IR2';


imcal = 'ImageCalibrated';
imtime = 'time';
colors= [0.9 0.2 0.8; 
         0.5 0.9 0.2; 
         0.2 0.2 0.2;
         0.9 0.5 0.8; 
         0.5 0.9 0.5; 
         0.2 0.5 0.5];


%Get times of all auroral events
load("\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structuresnew\IR2\marpeaksNH.mat");
load("\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structuresnew\IR2\marpeaksSH.mat");
load("\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structuresnew\IR2\febpeaksNH.mat");
load("\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structuresnew\IR2\febpeaksSH.mat");
load("\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structuresnew\IR2\aprpeaksNH.mat");
load("\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structuresnew\IR2\aprpeaksSH.mat");

full_tlist = [febNH.time,febSH.time,marNH.time,marSH.time,aprNH.time,aprSH.time];
full_rowlist = [febNH.row,febSH.row,marNH.row,marSH.row,aprNH.row,aprSH.row];
full_latlist = [febNH.lat,febSH.lat,marNH.lat,marSH.lat,aprNH.lat,aprSH.lat];

intensity_array_IR1 = [];
intensity_array_IR2 = [];

for filtern = 1:2
    if filtern == 1
        filter = 'IR1';
    else
        filter = 'IR2';
    end
for i = 1:length(full_tlist)
    if full_latlist(i)>0
        hemisphere = 'N';
    else
        hemisphere = 'S';
    end
    tmonth=month(full_tlist(i));
    switch tmonth
        case 2
            month_l = 'February';
        case 3 
            month_l = 'March';
        case 4
            month_l = 'April';
    end
    
    t = full_tlist(i)-minutes(5);
    while t < full_tlist(i)+minutes(4)
        filename = append('L1b_',filter,'-',num2str(year(t)),'_',num2str(month(t)),'_',num2str(day(t)),'_',num2str(hour(t)),'_',num2str(minute(t)),'_',num2str(second(t)),'-*');
        dirlist = dir(append('C:\Nobackup\juditpcj\MATS\Images_J\',filter,'\', month_l, '\',hemisphere,'H\',filename));
        if ~isempty(dirlist)
            break
        end
        t = t+seconds(1);
    end
    if isempty(dirlist)
        disp('The file has not been found');
        continue
    else
        disp('Found')
    end
    source = append('C:\Nobackup\juditpcj\MATS\Images_J\',filter,'\', month_l, '\',hemisphere,'H\',dirlist.name);
    
    vardata = ncread(source,imcal);
    vartime = ncread(source,imtime);
    times = vartime*10^(-9);
    time = datetime(2000,01,01,00,00,00) + seconds(times);
    [val,ind_col] = min(abs(time-full_tlist(i)));

    ir = 1;
    ic = 1;
    ind_row = full_rowlist(i);
    intensity_matrix = zeros([5,3]);
    for coln = 21:23
        nexttile; hold on; grid;
        % We need to create 3 different keograms
        keogramb = zeros([length(vardata(1,:,1)),length(vardata(1,1,:))]);
        for j = 1:length(vardata(1,1,:))
            mid_col = vardata(coln,:,j)';
            keogramb(:,j) = mid_col';
        end
        if ind_row>185
            ind_row = 185;
        end
        for rown = ind_row-2:ind_row+2 
            
            row = keogramb(rown,:);
            tsecs = seconds(time-time(1));
            [fit1,~,mu] = polyfit(tsecs,row',1);
            line1 = polyval(fit1,tsecs,[],mu);
            %Find max deviation of the background
            clear1 = row-line1';
            imax = find(clear1 == max(clear1(5:end-5)));
    
            %Find where aurora starts and ends to separate the background
            is = find(row(1:imax)'<line1(1:imax),1,'last')-1;
            ie = find(row(imax:end)'<line1(imax:end),1,'first')+imax+1;
            if isempty(ie)
                ie = length(row)-5;
            end
            if isempty(is) || is == 0
                is = 5;
            end
            y = [row(1:is) row(ie:end)];
            x = [tsecs(1:is)' tsecs(ie:end)'];
    
            %Fit again without the aurora to find background trend
            [fit2,~,mu] = polyfit(x',y',1);
            avg_woaurora = polyval(fit2,tsecs,[],mu);

            %Row witout background
            row_wob = row-avg_woaurora';

            %Fill intensity matrix
            intensity_matrix(ir,ic) = row_wob(ind_col);

            if ir == 5
                ir = 0;
            end
            ir = ir+1;
        end
        ic = ic+1;              
    end
    intensity_new = mean(intensity_matrix,'all');
    %Save in array
    if filter == 'IR1'
        intensity_array_IR1(i)= intensity_new;
    else
        intensity_array_IR2(i)= intensity_new;
    end
end
end


figure;hold on;
scatter(intensity_array_IR1,intensity_array_IR2)
xlabel('IR1')
ylabel('IR2')

figure;hold on;
scatter(full_tlist,intensity_array_IR2)
xlabel('IR1')
ylabel('IR2')

