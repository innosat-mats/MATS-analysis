% Intensities Re-run
close all; clear all;
retrieve = 0;
just_update = 1;
already_saved=0;

month_chosen = 'February';
hemisphere = 'S';
filter = 'IR2';

colors= [0.9 0.2 0.8; 
         0.5 0.9 0.2; 
         0.2 0.2 0.2;
         0.9 0.5 0.8; 
         0.5 0.9 0.5; 
         0.2 0.5 0.5];


%Get times of all auroral events
switch month_chosen
    case 'February'
        addpath("C:\Nobackup\juditpcj\MATS\Monthdata\Februarymonth\")
        load("febpeaksNH.mat");
        load("febpeaksSH.mat") ;
        full_tlist = [febpeaksNH.time,febpeaksSH.time];
        switch hemisphere
            case 'N'
                tevent_list = febpeaksNH.time;
                row_list = febpeaksNH.row;
                peaksall = febpeaksNH;
                if already_saved
                    load(append('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structuresnew\',filter,'\febpeaksNH.mat'))

                    peakssaved = febNH;
                    tevent_list = tevent_list(tevent_list>febNH.time(end))
                    row_list = row_list(tevent_list>febNH.time(end))
                end
                
            case 'S'
                tevent_list = febpeaksSH.time;
                row_list = febpeaksSH.row;
                peaksall = febpeaksSH;
                if already_saved
                    load(append('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structuresnew\',filter,'\febpeaksSH.mat'))
                    peakssaved = febSH;
                    tevent_list = tevent_list(tevent_list>peakssaved.time(end))
                    row_list = row_list(tevent_list>peakssaved.time(end))
                end
        end
    case 'March'
        addpath("C:\Nobackup\juditpcj\MATS\Monthdata\Marchmonth\")
        load("marpeaksNH.mat");
        load("marpeaksSH.mat") ;
        full_tlist = [marpeaksNH.time,marpeaksSH.time];
        switch hemisphere
            case 'N'
                tevent_list = marpeaksNH.time;
                row_list = marpeaksNH.row;
                peaksall = marpeaksNH;
                if already_saved
                    load(append('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structuresnew\',filter,'\marpeaksNH.mat'))
                    peakssaved = marNH;
                    tevent_list = tevent_list(tevent_list>peakssaved.time(end))
                    row_list = row_list(tevent_list>peakssaved.time(end))
                end
            case 'S'
                tevent_list = marpeaksSH.time;
                row_list = marpeaksSH.row;
                peaksall = marpeaksSH;
                if already_saved
                    load(append('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structuresnew\',filter,'\marpeaksSH.mat'))
                    peakssaved = marSH;
                    tevent_list = tevent_list(tevent_list>peakssaved.time(end))
                    row_list = row_list(tevent_list>peakssaved.time(end))
                end
        end
    case 'April'
        addpath("C:\Nobackup\juditpcj\MATS\Monthdata\Aprilmonth\")
        load("aprpeaksNH.mat");
        load("aprpeaksSH.mat");
        full_tlist = [aprpeaksNH.time,aprpeaksSH.time];
        switch hemisphere
            case 'N'
                tevent_list = aprpeaksNH.time;
                row_list = aprpeaksNH.row;
                peaksall = aprpeaksNH;
                if already_saved
                    load(append('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structuresnew\',filter,'\aprpeaksNH.mat'))
                    peakssaved = aprNH;
                    tevent_list = tevent_list(tevent_list>peakssaved.time(end))
                    row_list = row_list(tevent_list>peakssaved.time(end))
                end
            case 'S'
                tevent_list = aprpeaksSH.time;
                row_list = aprpeaksSH.row;
                peaksall = aprpeaksSH;
                if already_saved
                    load(append('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structuresnew\',filter,'\aprpeaksSH.mat'))
                    peakssaved = aprSH;
                    tevent_list = tevent_list(tevent_list>peakssaved.time(end))
                    row_list = row_list(tevent_list>peakssaved.time(end))
                end
        end
    
end

%% Download images
if retrieve 
    getzarr = 'python \\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\get_zarr.py';
    changedirpath = append('C:\Nobackup\juditpcj\MATS\Images_J\',filter,'\', month_chosen, '\',hemisphere,'H\');
    cd(changedirpath)
    for i = 1:length(tevent_list)
        disp(tevent_list(i));
        ti = tevent_list(i)-minutes(4);
        te = tevent_list(i)+minutes(4);
        str_ti = append('2023 ',num2str(month(ti)), ' ',num2str(day(ti)), ' ', num2str(hour(ti)), ' ', num2str(minute(ti)), ' ', num2str(second(ti)));
        str_te = append('2023 ',num2str(month(te)), ' ',num2str(day(te)), ' ', num2str(hour(te)), ' ', num2str(minute(te)), ' ', num2str(second(te)));
        
        call_python = append(getzarr,' -c IR2 -b ',str_ti,' -e ', str_te)
        system(call_python);
    end
    cd \\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\
end

%% Select right events,fix intensity and save in list with only them
if just_update
    dirlist = dir(append('C:\Nobackup\juditpcj\MATS\Images_J\',filter,'\', month_chosen, '\',hemisphere,'H\'));
    imcal = 'ImageCalibrated';
    imtime = 'time';
    tevent_list = datetime(2023,month(tevent_list),day(tevent_list),hour(tevent_list),minute(tevent_list),...
        second(tevent_list));
    intensity_array = zeros([length(dirlist),1]);
    events_right = zeros([length(dirlist),1]);
    for i = 1:length(dirlist)-2
        namefile = dirlist(i+2).name
        source = append('C:\Nobackup\juditpcj\MATS\Images_J\',filter,'\',month_chosen,'\',hemisphere,'H\',namefile);
        vardata = ncread(source,imcal);
        vartime = ncread(source,imtime);
        times = vartime*10^(-9);
        time = datetime(2000,01,01,00,00,00) + seconds(times);
        if already_saved ==1 && time(1) < peakssaved.time(end)
            continue
        end
        ind_event = find(tevent_list> time(round(end/2)-10) & tevent_list< time(round(end/2)+10),1,"first");
        ind_row = row_list(ind_event);
        [val,ind_col] = min(abs(time-tevent_list(ind_event)));
        if isempty(ind_event)
            continue
        end

        %Generate and plot keogram
        keogram = zeros([length(vardata(1,:,1)),length(vardata(1,1,:))]);
        for j = 1:length(vardata(1,1,:))
            mid_col = vardata(22,:,j)';
            keogram(:,j) = mid_col';
        end
        Aurora_good =0;
        savenow = 0;
        while Aurora_good ==0
            %Figure
            pxi = 100; pyi = 100; pxf = 1500; pyf = 1200;
            figure(Position= [pxi pyi pxf pyf]);
            tiledlayout(2,1)
            nexttile; hold on;
            surf(time,1:1:length(mid_col),keogram,'edgecolor','none')
            view([0 90]);
            xline(time(ind_col),'r--')
            yline(ind_row,'r')
            xlim([time(1) time(end)])
            ylim([1 length(mid_col)])
            nexttile;
            
            keogram_wob = keogram;
            for rown = 130:length(keogram(:,1))
                row_raw = keogram(rown,:);
                row = smoothdata(row_raw,'gaussian',4);
                tsecs = seconds(time-time(1));
            
                [fit1,~,mu] = polyfit(tsecs,row',3);
                line1 = polyval(fit1,tsecs,[],mu);
                row_wob = row-line1';
                keogram_wob(rown,:) = row_wob;
            end
            imagesc(time,1:1:length(mid_col),keogram_wob,[0 1e14])
            xline(time(ind_col),'r--')
            yline(ind_row,'r')
            xlim([time(1) time(end)])
            ylim([1 length(mid_col)])
            set(gca,'YDir','normal');

            
            
            %Decide if it is aurora or not
            waitforbuttonpress;            
            key = double(get(gcf, 'CurrentCharacter'))
            switch key
                case 13 % Enter -> yes aurora
                    Aurora_good = 1;
                    disp('Good')
                case 29 % Right arrow -> col to right
                    ind_col = ind_col+1;                    
                case 28 % Left arrow -> col to left
                    ind_col = ind_col-2;                    
                case 31 % Down arrow
                    ind_row = ind_row-5;
                case 30 % Up arrow
                    ind_row = ind_row+2;
                case 8 % Erase -> not aurora
                    disp('NO');
                    break
                otherwise
                    savenow = 1;
                    break
            end
        end
        close all
        if savenow ==1
            break
        end
        if Aurora_good == 0
            continue
        end
        if ind_row >185
            ind_row = 185;
        end
        % Remove background for the wanted row and the ones before and
        % after
        ir = 1;
        ic = 1;
        intensity_matrix = zeros([5,3]);
        %figure(Position= [pxi pyi pxf pyf]);
        %tiledlayout(4,1)
        for coln = 21:23
            nexttile; hold on; grid;
            % We need to create 3 different keograms
            keogramb = zeros([length(vardata(1,:,1)),length(vardata(1,1,:))]);
            for j = 1:length(vardata(1,1,:))
                mid_col = vardata(coln,:,j)';
                keogramb(:,j) = mid_col';
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
                
                %Plot
                %scatter(time, row,'*',MarkerEdgeColor=colors(ic,:))
                %plot(time,line1','--',color= colors(ic,:))
                %plot(time,avg_woaurora', color = colors(ic,:))
                %xline(time(is),color = colors(ic,:))
                %xline(time(ie),color = colors(ic,:))
                %xline(time(ind_col),'--')
                %xline(time(imax),'r.-')
    
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
        

        %Get new value for intensity taking average of 9 values surrounding
        %pixel
        intensity_new = mean(intensity_matrix,'all');
        %Save in array
        intensity_array(ind_event)= intensity_new;
        events_right(ind_event)=1;
    end
    %doublecheck it worked and save in structure
    if already_saved
        skipped = zeros(length(peaksall.alt) - length(peakssaved.alt),1)';
        %events_add = [skipped events_right'];
        events_add = events_right;

        newsave.newI   = [peakssaved.newI   intensity_array(events_add==1)'];
        newsave.row    = [peakssaved.row    peaksall.row(events_add==1)   ];
        newsave.alt    = [peakssaved.alt    peaksall.alt(events_add==1)   ];
        newsave.maxlat = [peakssaved.maxlat peaksall.maxlat(events_add==1)];
        newsave.maxlon = [peakssaved.maxlon peaksall.maxlon(events_add==1)];
        newsave.lat    = [peakssaved.lat    peaksall.lat(events_add==1)   ];
        newsave.time   = [peakssaved.time   peaksall.time(events_add==1)  ];
        newsave.MLT    = [peakssaved.MLT    peaksall.MLT(events_add==1)   ];
        newsave.Mlat   = [peakssaved.Mlat   peaksall.Mlat(events_add==1)  ];
        newsave.Mlon   = [peakssaved.Mlon   peaksall.Mlon(events_add==1)  ];
        newsave.kp     = [peakssaved.kp     peaksall.kp(events_add==1)    ];
    else
        newsave.newI   = intensity_array(events_right==1)';
        newsave.row    = peaksall.row(events_right==1)   ;
        newsave.alt    = peaksall.alt(events_right==1)   ;
        newsave.maxlat = peaksall.maxlat(events_right==1);
        newsave.maxlon = peaksall.maxlon(events_right==1);
        newsave.lat    = peaksall.lat(events_right==1)   ;
        newsave.time   = peaksall.time(events_right==1)  ;
        newsave.MLT    = peaksall.MLT(events_right==1)   ;
        newsave.Mlat   = peaksall.Mlat(events_right==1)  ;
        newsave.Mlon   = peaksall.Mlon(events_right==1)  ;
        newsave.kp     = peaksall.kp(events_right==1)    ;
    end


    disp('Saving')
    %switch month_chosen
    %    case 'February'
    %        switch hemisphere
    %            case 'N'
    %                febNH = newsave;
    %                save(append('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structuresnew\',filter,'\febpeaksNH.mat'),'febNH')
    %            case 'S'
    %                febSH = newsave;
    %                save(append('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structuresnew\',filter,'\febpeaksSH.mat'),'febSH')
    %        end
    %    case 'March'
    %        switch hemisphere
    %            case 'N'
    %                marNH = newsave;
    %                save(append('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structuresnew\',filter,'\marpeaksNH.mat'),'marNH')
    %            case 'S'
    %                marSH = newsave; 
    %                save(append('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structuresnew\',filter,'\marpeaksSH.mat'),'marSH')
    %        end
    %    case 'April'
    %        switch hemisphere
    %            case 'N'
    %                aprNH = newsave; 
    %                save(append('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structuresnew\',filter,'\aprpeaksNH.mat'),'aprNH')
    %            case 'S'
    %                aprSH = newsave;
    %                save(append('\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\Structuresnew\',filter,'\aprpeaksSH.mat'),'aprSH')
    %        end
    %end
    if length(nonzeros(intensity_array)) == length(intensity_array)
        disp('All the events were found')
    else
        events_missing = length(intensity_array) - length(nonzeros(intensity_array))
        disp(append(num2str(events_missing),' out of ',num2str(length(intensity_array)),' events were not found, some intensities will be 0.'))
        
    end
       
end
