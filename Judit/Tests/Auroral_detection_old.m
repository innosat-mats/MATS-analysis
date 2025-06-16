% Intensities Re-run
close all; %clear all;

%Write 1 to the things you want to do and 0 if not:
retrieve = 1;
remove_lat_daylight = 0;
find_aurora=0;
remove_background = 0;
monthdo = 'mar';
plotyesno = 0;



colors= [0.9 0.2 0.8; 
         0.5 0.9 0.2; 
         0.2 0.2 0.2;
         0.9 0.5 0.8; 
         0.5 0.9 0.5; 
         0.2 0.5 0.5];

%% Download images for all three months
if retrieve 
    disp('Retrieving data')
    getzarr = 'python \\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\get_zarr.py';
    changedirpath = 'C:\Nobackup\juditpcj\MATS\Images_J\All\';
    cd(changedirpath)
    ti = datetime(2023,03,01,00,00,00);
    te = datetime(2023,04,01,00,00,00);
    str_ti = append('2023 ',num2str(month(ti)), ' ',num2str(day(ti)), ' ', num2str(hour(ti)), ' ', num2str(minute(ti)), ' ', num2str(second(ti)));
    str_te = append('2023 ',num2str(month(te)), ' ',num2str(day(te)), ' ', num2str(hour(te)), ' ', num2str(minute(te)), ' ', num2str(second(te)));
    
    call_python = append(getzarr,' -c IR1 -b ',str_ti,' -e ', str_te)
    system(call_python);
    cd \\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\
    disp('Done')
end

%% Remove latitudes below 40º, separate data by hemispheres and save
if remove_lat_daylight
    minlat = 45;
    %source = 'C:\Nobackup\juditpcj\MATS\Images_J\All\L1b_IR1-2023_2_1_0_0_0-2023_4_30_0_0_0.nc';
    %source = '\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\L1b_IR1-2023_4_23_0_0_0-2023_4_24_0_0_0.nc'
    switch monthdo
        case 'feb'
            source = 'C:\Nobackup\juditpcj\MATS\Images_J\All\L1b_IR2-2023_2_8_0_0_0-2023_2_28_0_0_0.nc';
        case 'mar'
            source = 'C:\Nobackup\juditpcj\MATS\Images_J\All\L1b_IR2-2023_2_28_0_0_0-2023_4_1_0_0_0.nc';
        case 'apr'
            source = 'C:\Nobackup\juditpcj\MATS\Images_J\All\L1b_IR2-2023_4_1_0_0_0-2023_5_1_0_0_0.nc';
    end
    
    imcal = 'ImageCalibrated';
    imtime = 'time';
    imtplat = 'TPlat';
    imtplon = 'TPlon';
    imtpalt = 'TPheight';
    %imaltpix = 'TPheightPixel';
    vardata = ncread(source,imcal);
    vartime = ncread(source,imtime);
    varlat = ncread(source,imtplat);
    varlon = ncread(source,imtplon);
    %varaltpix = ncread(source,imaltpix);
    times = vartime*10^(-9);
    time = datetime(2000,01,01,00,00,00) + seconds(times);

    %Separate by hemisphere and cut anything more equatorward than 40
    %degrees
    NH.lat = varlat(varlat>minlat);
    NH.images = vardata(:,:,varlat>minlat);
    NH.time = time(varlat>minlat);
    NH.lon = varlon(varlat>minlat);
    %NH.altpix=varaltpix(:,:,varlat>minlat);
    SH.lat = varlat(varlat<-minlat);
    SH.images = vardata(:,:,varlat<-minlat);
    SH.time = time(varlat<-minlat);
    SH.lon = varlon(varlat<-minlat);
    %SH.altpix=varaltpix(:,:,varlat<-minlat);

    goodlatsNSH.lat = varlat(varlat>minlat | varlat<-minlat);
    goodlatsNSH.images = vardata(:,:,varlat>minlat | varlat<-minlat);
    goodlatsNSH.time = time(varlat>minlat | varlat<-minlat);
    goodlatsNSH.lon = varlon(varlat>minlat | varlat<-minlat);
    %goodlatsNSH.altpix=varaltpix(:,:,varlat>minlat & varlat<-minlat);

    save('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structures\goodlats.mat','NH','SH','goodlatsNSH')
end

if find_aurora
    %load('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structures\goodlats.mat')
    %Find cuts larger than 1 minute
    ind_cut = find(seconds(diff(NH.time))>60);
    ind_cut = [0 ind_cut' length(NH.time)];
    list_taurora = [];
    
    %For each pass
    for i=1:length(ind_cut)-1
        %Create keogram
        keogram = zeros([length(NH.images(1,:,1)),length(NH.images(1,1,ind_cut(i)+1:ind_cut(i+1)))]);
        times_aurora = [];
        t_keo = NH.time(ind_cut(i)+1:ind_cut(i+1));
        for j = ind_cut(i)+1:ind_cut(i+1)
            mid_col = NH.images(22,:,j)';
            keogram(:,j-ind_cut(i)) = mid_col';
        end

        %Plot keogram
        if plotyesno
            pxi = 100;pyi = 100;pxf = 1200;pyf = 600;  
            figure(Position= [pxi pyi pxf pyf]);
            tiledlayout(2,1)
            nexttile; hold on;
            imagesc(t_keo,1:length(NH.images(1,:,1)),keogram);
            set(gca,'YDir','normal');
        end
        %title(append('Latitude = ',num2str(goodlatsNSH.lat(ind_cut(i)+1))))
        if plotyesno 
            nexttile; hold on;
        end
        %FIND AURORA
        ic = 1;
        dev = 1.15;
        Conditions = zeros(10,1);
        %Evaluate some lines: two below pix 70, two between 70 and 130
        %and two between 130 and 180
        for rown = [140,150,160]%1:length(keogram(:,1))
            row = keogram(rown,:);
            tsecs = seconds(t_keo-t_keo(1));
           
            [fit1,~,mu] = polyfit(tsecs,row',1);
            line1 = polyval(fit1,tsecs,[],mu);
            %Find max deviation of the background
            clear1 = row-line1';
            imax = find(clear1 == max(clear1));
            
            if plotyesno
                scatter(t_keo, row,'*',MarkerEdgeColor=colors(ic,:))
                plot(t_keo,line1','-',color= colors(ic,:))
                plot(t_keo,line1'*dev,'--',color= colors(ic,:))
                scatter(t_keo(imax),row(imax),'o',MarkerEdgeColor=colors(ic,:))
                text(t_keo(end)+seconds(5),row(end),num2str(Conditions(ic)),color = colors(ic,:))
                pause
            end
            
            %Condition 1: If there is no big deviation from the average, 
            % there is no aurora.
            if row(imax) > line1(imax)*dev
                Conditions(ic)=1;
                times_aurora = [times_aurora t_keo(row'>line1*dev)'];
            end
            
            ic = ic+1;
        end
        if length(nonzeros(Conditions))<3
            continue % No aurora in strip, go to next iteration
        end
        title('YES')
        %Remove times not part of aurora
        ta_sorted = sort(times_aurora);
        ind_split = [0 find(seconds(diff(ta_sorted))>120) length(ta_sorted)];
        %ind_split_max = find(diff(ind_split) == max(diff(ind_split)));
        %group_t =ta_sorted(ind_split(ind_split_max):ind_split(ind_split_max+1));
        for splitn = 2:length(ind_split)
            group_t =ta_sorted(ind_split(splitn-1)+1:ind_split(splitn));                
            if length(group_t)>10
                average_splitn = mean(group_t);
                list_taurora = [list_taurora average_splitn];
            end
        end
    end
    switch monthdo
        case 'feb'
            list_taurora_feb = list_taurora;
            save('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\list_taurora_J.mat','list_taurora_feb')
        case 'mar'
            list_taurora_mar = list_taurora;
            save('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\list_taurora_J.mat','list_taurora_mar')
        case 'apr'
            list_taurora_apr = list_taurora;
            save('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\list_taurora_J.mat','list_taurora_apr')
    end
end

