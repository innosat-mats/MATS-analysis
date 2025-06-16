% Intensities Re-run
close all;

%Write 1 to the things you want to do and 0 if not:

find_aurora=1;
load_data = 0;
load_pixels = 0;

monthdo = 'mar';
hemisphere = 'N';

airg_top = 130;

colors= [0.9 0.2 0.8; 
         0.5 0.9 0.2; 
         0.2 0.2 0.2;
         0.9 0.5 0.8; 
         0.5 0.9 0.5; 
         0.2 0.5 0.5];

if load_pixels
    load('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\list_taurora_J.mat')
    switch monthdo
        case 'feb'
            if hemisphere == 'N'
                t_start = list_taurora_feb.NH.time(end);
                pixel_tim_list  = list_taurora_feb.NH.time;
                pixel_row_list  = list_taurora_feb.NH.row;
                altres_tim_list = list_taurora_feb.NH.altrestime;
                altres_row_list = list_taurora_feb.NH.altresrow;
            else
                t_start = list_taurora_feb.SH.time(end);
                pixel_tim_list  = list_taurora_feb.SH.time;
                pixel_row_list  = list_taurora_feb.SH.row;
                altres_tim_list = list_taurora_feb.SH.altrestime;
                altres_row_list = list_taurora_feb.SH.altresrow;

            end
        case 'mar'
            if hemisphere == 'N'
                t_start = list_taurora_mar.NH.time(end);
                pixel_tim_list  = list_taurora_mar.NH.time;
                pixel_row_list  = list_taurora_mar.NH.row;
                altres_tim_list = list_taurora_mar.NH.altrestime;
                altres_row_list = list_taurora_mar.NH.altresrow;
            else
                t_start = list_taurora_mar.SH.time(end);
                pixel_tim_list  = list_taurora_mar.SH.time;
                pixel_row_list  = list_taurora_mar.SH.row;
                altres_tim_list = list_taurora_mar.SH.altrestime;
                altres_row_list = list_taurora_mar.SH.altresrow;
            end
        case 'apr'
            if hemisphere == 'N'
                t_start = list_taurora_apr.NH.time(end);
                pixel_tim_list  = list_taurora_apr.NH.time;
                pixel_row_list  = list_taurora_apr.NH.row;
                altres_tim_list = list_taurora_apr.NH.altrestime;
                altres_row_list = list_taurora_apr.NH.altresrow;
            else
                t_start = list_taurora_apr.SH.time(end);
                pixel_tim_list  = list_taurora_apr.SH.time;
                pixel_row_list  = list_taurora_apr.SH.row;
                altres_tim_list = list_taurora_apr.SH.altrestime;
                altres_row_list = list_taurora_apr.SH.altresrow;
            end
    end
end

%% Remove latitudes below 45º, separate data by hemispheres
if load_data
    minlat = 45;
    switch monthdo
        case 'feb'
            source = 'C:\Nobackup\juditpcj\MATS\Images_J\All\L1b_IR2-2023_2_8_0_0_0-2023_2_28_0_0_0.nc';
        case 'mar'
            source = 'C:\Nobackup\juditpcj\MATS\Images_J\All\L1b_IR2-2023_2_28_0_0_0-2023_4_1_0_0_0.nc';
        case 'apr'
            source = 'C:\Nobackup\juditpcj\MATS\Images_J\All\L1b_IR2-2023_4_1_0_0_0-2023_5_1_0_0_0.nc';
    end
    %source = '\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\L1b_IR1-2023_4_23_0_0_0-2023_4_24_0_0_0.nc'
    
    imcal = 'ImageCalibrated';
    imtime = 'time';
    imtplat = 'TPlat';
    imtplon = 'TPlon';
    imtpalt = 'TPheight';
    %imaltpix = 'TPheightPixel';
    disp('Reading file')
    vardata = ncread(source,imcal);
    vartime = ncread(source,imtime);
    varlat = ncread(source,imtplat);
    varlon = ncread(source,imtplon);
    %varaltpix = ncread(source,imaltpix);
    times = vartime*10^(-9);
    time = datetime(2000,01,01,00,00,00) + seconds(times);

    %Separate by hemisphere and cut anything more equatorward than 45
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
    disp('Saving file')
    %save('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structures\goodlats.mat','NH','SH')
end
if find_aurora
    %load('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structures\goodlats.mat')   
    load('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\list_taurora_J.mat')
    
    switch hemisphere
        case 'N'
            disp('Starting analysis for NH')
            Htime = NH.time;
            Himages = NH.images;
            Hlat = NH.lat;
            
        case 'S'
            disp('Starting analysis for SH')
            Htime = SH.time;
            Himages = SH.images;
            Hlat = SH.lat;
    end
    if load_pixels
        [dif1,ind_start] = min(abs(Htime-t_start));
        Htime = Htime(ind_start:end);
        Himages = Himages(:,:,ind_start:end);
        Hlat = Hlat(ind_start:end);
    else
        pixel_tim_list = [];
        pixel_row_list = [];
        altres_tim_list = [];
        altres_row_list = [];
    end
    
    %Find cuts longer than 1 minute
    ind_cut = find(seconds(diff(Htime))>60);
    ind_cut = [0 ind_cut' length(Htime)];
    ind_cut2 = ind_cut;
    for i = 1:length(ind_cut)-1
        ind_cut2(i) = round((ind_cut(i+1)-ind_cut(i))*0.5)+ind_cut(i) ;
    end
    ind_cuts = sort([ind_cut ind_cut2]);
    ind_cut = ind_cuts';

    %For each half pass
    for i=1:length(ind_cut)-1
        close all
        %Create keogram
        keogram = zeros([length(Himages(1,:,1)),length(Himages(1,1,ind_cut(i)+1:ind_cut(i+1)))]);
        
        times_aurora = [];
        t_keo = Htime(ind_cut(i)+1:ind_cut(i+1));
        for j = ind_cut(i)+1:ind_cut(i+1)
            mid_col = Himages(22,:,j)';
            keogram(:,j-ind_cut(i)) = mid_col';            
        end

        if length(t_keo)<20
            continue
        end

        %Plot keogram
        pxi = 100;pyi = 100;pxf = 1600;pyf = 1000;  
        figure(Position= [pxi pyi pxf pyf]);
        tiledlayout(3,1,'TileSpacing','compact')
        ax1=nexttile; hold on;
        imagesc(ax1,t_keo,1:length(Himages(1,:,1)),keogram);
        %yline(120)
        set(gca,'YDir','normal');
        ylim(ax1,[1 length(keogram(:,1)) ])
        xlim(ax1,[t_keo(1) t_keo(end)])
        %title(append('Latitude = ',num2str(Hlat(ind_cut(i)+1))))
        

        keogram_wob = keogram;
        ax2 = nexttile; hold on;legend();
        ax3 = nexttile; hold on;
        %FIND AURORA
        ic = 0;
        %tmaxs = 
        conditions = zeros(length(keogram(:,1))-airg_top,1);
        ie_list = conditions;
        is_list = conditions;
        Aurora_yes = 1;
        for rown = airg_top:length(keogram(:,1))
            row_raw = keogram(rown,:);
            row = smoothdata(row_raw,'gaussian',4);
            tsecs = seconds(t_keo-t_keo(1));
           
            [fit1,~,mu] = polyfit(tsecs,row',3);
            line1 = polyval(fit1,tsecs,[],mu);
            %Find max deviation from background
            clear1 = row-line1';
            imax = find(clear1 == max(clear1));
            tmax = t_keo(imax);
            
            dev = mean(line1)*0.2;


            if ismember(rown,[135 140 160 180 186])
                ic = ic+1;
                scatter(ax2,t_keo, row,'*',MarkerEdgeColor=colors(ic,:),DisplayName=num2str(rown))
                plot(ax2,t_keo,line1','-',color= colors(ic,:),HandleVisibility='off')
                plot(ax2,t_keo,line1'+dev,'--',color= colors(ic,:),HandleVisibility='off')
                %plot(ax2,t_keo,avg_woaurora','.-',color= colors(ic,:))
                scatter(ax2,t_keo(imax),row(imax),'o',MarkerEdgeColor=colors(ic,:),HandleVisibility='off') 
                xlim(ax2,[t_keo(1) t_keo(end)])
            end
            
            if row(imax)<line1(imax)+dev
                row_wob = row-line1';
                keogram_wob(rown,:) = row_wob;
                continue
            end
            
            %Find where aurora starts and ends to separate the background
            is = find(row(1:imax)'<line1(1:imax),1,'last')-1;
            ie = find(row(imax:end)'<line1(imax:end),1,'first')+imax;
            if isempty(ie) || isempty(is) || is == 0 || ie>length(row)
                row_wob = row-line1';
                keogram_wob(rown,:) = row_wob;
                continue
            end
            %There might be aurora
            conditions(rown-airg_top+1)=1;
            y = [row(1:is) row(ie:end)];
            x = [tsecs(1:is)' tsecs(ie:end)'];
            [fit2,~,mu] = polyfit(x',y',1);
            avg_woaurora = polyval(fit2,tsecs,[],mu);
            %Row and keogram without background
            row_wob = row-avg_woaurora';
            keogram_wob(rown,:) = row_wob;   
            is_list(rown) = is;
            ie_list(rown) = ie;
        
            if ismember(rown,[135 140 160 180 186])
                xline(ax2,t_keo(is),color = colors(ic,:),HandleVisibility='off')
                xline(ax2,t_keo(ie),color = colors(ic,:),HandleVisibility='off')
            end

            
        end
        %pause(0.6)
        if length(nonzeros(conditions)) <length(conditions)/2
            continue
        end 

        is = floor(median(nonzeros(is_list)));
        ie = floor(median(nonzeros(ie_list)));
        xline(ax2,t_keo(is),color = 'k', linewidth = 2)
        xline(ax2,t_keo(ie),color = 'k', linewidth = 2)

        imagesc(ax3,t_keo,1:length(Himages(1,:,1)),keogram_wob,[0 max(keogram_wob,[],'all')/2]);
        %Find max intensity of aurora
        [max_above_nightglow,row_max] = max(keogram_wob(130:end-2,is:ie));
        scatter(ax3,t_keo(is:ie),row_max+129,'r*')
        set(gca,'YDir','normal');
        ylim(ax3,[1 length(keogram(:,1)) ])
        xlim(ax3,[t_keo(1) t_keo(end)]) 
        
        %[pixel_row,pixel_tim] = findpeaks(row_max,t_keo(is:ie),'MinPeakProminence',5);
        row_max_sort =sort(unique(row_max),'descend');
        if length(row_max_sort)<3
            Aurora_yes = 0;
            continue
        end

        for j = 1:length(row_max_sort)-1
            pixel_row = row_max_sort(j);
            if j == length(row_max_sort)-1 
                Aurora_yes = 0;
                break
            end
            
            if row_max_sort(j) > row_max_sort(j+1)+3
                continue
            end
            pixel_ind = find(row_max == pixel_row);
            for k = 1:length(pixel_ind)
                Aurora_yes =0;
                if pixel_ind(k) == 1 || pixel_ind(k) == length(row_max)
                    continue
                end
                if row_max(pixel_ind(k)) > row_max(pixel_ind(k)-1)+3 || row_max(pixel_ind(k)) > row_max(pixel_ind(k)+1)+3
                    continue
                end
                if row_max(pixel_ind(k)) == row_max(pixel_ind(k)-1) && row_max(pixel_ind(k)) == row_max(pixel_ind(k)+1)
                    continue
                end
                Aurora_yes = 1;
                pixel_tim = t_keo(pixel_ind+is-1);
                if row_max(pixel_ind(k)) <= min(row_max(pixel_ind(k):end))
                    continue
                end
                
                break            
            end
            if Aurora_yes
                break
            else
                continue
            end
            
        end
        if Aurora_yes == 0
            continue
        end
        
        if length(pixel_ind)>1
            [dif,index] = min(abs(pixel_tim-median(t_keo(is:ie))));
            pixel_tim = pixel_tim(index);
            %pixel_row = pixel_row(index);
        end

        plot(ax3,pixel_tim,pixel_row+129,'o',LineWidth=3)
        linkaxes([ax1,ax2,ax3],'x')
        

        waitforbuttonpress;
        key = double(get(gcf, 'CurrentCharacter'))
        switch key
            case 29 % Right arrow -> save
                pixel_tim_list = [pixel_tim_list pixel_tim];
                pixel_row_list = [pixel_row_list pixel_row+129];
                disp('El pixel sha guardat');
            case 28 % Left arrow -> save to list of doubtful
                altres_tim_list = [altres_tim_list pixel_tim];
                altres_row_list = [altres_row_list pixel_row+129];
                disp('El pixel sha guardat a la llista dolenta')
            case 31 % Down arrow -> not save
                disp('El pixel NO sha guardat');
            otherwise % Otherwise -> Save and stop
                break 
        end
         
    end
    switch hemisphere
        case 'N'
            list_taurora.NH.row = pixel_row_list;
            list_taurora.NH.time = pixel_tim_list;
            list_taurora.NH.altresrow = altres_row_list;
            list_taurora.NH.altrestime = altres_tim_list;
        case 'S'
            list_taurora.SH.row = pixel_row_list;
            list_taurora.SH.time = pixel_tim_list;
            list_taurora.SH.altresrow = altres_row_list;
            list_taurora.SH.altrestime = altres_tim_list;
    end
    
    disp('Guardant al document')
    switch monthdo
        case 'feb'
            list_taurora_feb = list_taurora;
            save('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\list_taurora_J.mat','list_taurora_feb','-append')
        case 'mar'
            list_taurora_mar = list_taurora;
            save('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\list_taurora_J.mat','list_taurora_mar','-append')
        case 'apr'
            list_taurora_apr = list_taurora;
            save('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\list_taurora_J.mat','list_taurora_apr','-append')
    end
end



