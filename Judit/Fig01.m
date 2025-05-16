clear all; close all;
source = '\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub\MATS-analysis\Judit\Downloadfiles\L1b_IR1-2023_4_23_0_0_0-2023_4_24_0_0_0.nc';
imcal = 'ImageCalibrated';
imtime = 'time';
savefig = '\\ug.kth.se\dfs\home\j\u\juditpcj\appdata\xp.V2\Documents\GitHub datafiles\MATS\Fig01.png';

%input
t_start = datetime(2023,04,23,11,28,09);
t_end = datetime(2023,04,23,11,45,04);
rown1 = 40 ; rown2 = 120; rown3 = 170;
c2= [0.9 0.2 0.8];
c3= [0.5 0.9 0.2];
c1= [0.2 0.2 0.2];

varinfo = ncinfo(source);
varimages = ncread(source,imcal);
vartime = ncread(source,imtime);
vartplat = ncread(source,'TPlat');
times = vartime*10^(-9);
time = datetime(2000,01,01,00,00,00) + seconds(times);
index_ts = find(time > t_start, 1, 'first');
index_te = find(time > t_end, 1, 'first');

TPlat = vartplat(index_ts:index_te);
timea = time(index_ts:index_te);

%keogram_i = zeros([187 14311]);
keogram_i = zeros([length(varimages(1,:,1)) index_te-index_ts]);
for j = index_ts:index_te
    i = j - index_ts+1;
    mid_row = varimages(round(end/2),:,j)';
    keogram_i(:,i) = mid_row';
end

%keogram_i=keogram_i(end:-1:1,:);

%%
im1t = datetime(2023,04,23,11,33,00);
im2t = datetime(2023,04,23,11,36,10);
im3t = datetime(2023,04,23,11,36,55);
im4t = datetime(2023,04,23,11,40,00);

index_im1t = find(time > im1t, 1, 'first');
index_im2t = find(time > im2t, 1, 'first');
index_im3t = find(time > im3t, 1, 'first');
index_im4t = find(time > im4t, 1, 'first');

im1 = varimages(:,:,index_im1t)';
im2 = varimages(:,:,index_im2t)';
im3 = varimages(:,:,index_im3t)';
im4 = varimages(:,:,index_im4t)';

%% Background subtraction

keogram_wb = zeros(size(keogram_i));
for rown = 1:length(keogram_i(:,1))
    row = keogram_i(rown,:);
    m = mean(row);
    is = find(row(10:end)>m,1,'first')-2;
    ie = find(row>m,1,'last')+2;
    
    avg_wa= mean([row(1:is) row(ie:end)]);
    
    switch rown
        case rown1 
            avg1 = avg_wa;
        case rown2
            avg2 = avg_wa;
        case rown3
            avg3 = avg_wa;        
    end

    keogram_wb(rown,:) = keogram_i(rown,:)-avg_wa;

    %if mod(rown,20)==0
    %    figure();hold on;
    %    scatter(timea,row)
    %    xline(timea(is))
    %    xline(timea(ie))
    %    yline(m)
    %    yline(avg_wa,'r')
    %end
end

row1 = keogram_i(rown1,:);
row2 = keogram_i(rown2,:);
row3 = keogram_i(rown3,:);

fig=figure(Position=[10 20 1200 850]);
%tiledlayout(3,4,'TileSpacing', 'compact', 'Padding', 'compact')

ax11 = subplot(3,4,1); hold on; 
title(datestr(im1t,'HH:MM:SS'))
imagesc(im1, [0 2e14])
xline(length(im1(1,:))/2, 'r', linewidth = 2)
xlim([0 length(im1(1,:))]);ylim([0 length(im1(:,1))]);
set(gca,'xtick',[])
set(gca,'ytick',[])
xlim([1 44])

ax12 = subplot(3,4,2); hold on; 
title(datestr(im2t,'HH:MM:SS'))
imagesc(im2, [0 2e14])
xline(length(im2(1,:))/2, 'r', linewidth = 2)
xlim([0 length(im2(1,:))]);ylim([0 length(im2(:,1))]);
set(gca,'xtick',[])
set(gca,'ytick',[])
xlim([1 44])

ax13 = subplot(3,4,3); hold on; 
title(datestr(im3t,'HH:MM:SS'))
imagesc(im3, [0 2e14])
xline(length(im1(1,:))/2, 'r', linewidth = 2)
xlim([0 length(im1(1,:))]);ylim([0 length(im1(:,1))]);
set(gca,'xtick',[])
set(gca,'ytick',[])
xlim([1 44])

ax14 = subplot(3,4,4); hold on; 
title(datestr(im4t,'HH:MM:SS'))
imagesc(im4, [0 2e14])
xline(length(im1(1,:))/2, 'r', linewidth = 2)
xlim([0 length(im1(1,:))]);ylim([0 length(im1(:,1))]);
set(gca,'xtick',[])
set(gca,'ytick',[])
xlim([1 44])

%keogram
ax1=subplot(3,4,[5 8]); hold on; 
imagesc(timea, 1:187, keogram_i, [0 2e14])
xline(im1t, 'r', linewidth = 2)
xline(im2t, 'r', linewidth = 2)
xline(im3t, 'r', linewidth = 2)
xline(im4t, 'r', linewidth = 2)
yline(rown1,'--', color = c1, linewidth = 1.5)
yline(rown2,'--', color = c2, linewidth = 1.5)
yline(rown3,'--', color = c3, linewidth = 1.5)
xlim([timea(1) timea(end)])
ylim([1 187])

ticklabelstim = timea(1):minutes(5):timea(end);
ticklabelslat = TPlat(ismember(timea,ticklabelstim));
labelArray = [compose('%5s',datestr(ticklabelstim,'HH:MM'))';compose('%.fº',ticklabelslat')];
labelArray = strjust(pad(labelArray),'center');
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
ax = gca(); 
ax.XTick = timea(1):minutes(3):timea(end); 
ax.XLim = [timea(1) timea(end)];
ax.XTickLabel = tickLabels; 
text(ax1,'Units', 'Normalized','Position', [0.95, -0.1],'string','Time', color = 'k')
text(ax1,'Units', 'Normalized','Position', [0.95, -0.25],'string','TPlat', color = 'k')

%xticks_current = get(gca, 'xtick');
%xticklabels_new = cellstr(datestr(xticks_current, 'HH:MM'));
%set(gca, 'xticklabel', xticklabels_new)
%ax2 = axes('Position', ax1.Position, 'Color', 'none');
%x2 = TPlat;
%ax2.XLim = [min(x2), max(x2)]; % Escala independent
%ax2.YLim = ax1.YLim; % Comparteix escala Y
%ax2.YTick = [];   % Eliminar marques Y
%ax2.YColor = 'none';
%ax2.XAxisLocation = 'bottom'; % Assegurar que estigui a baix
%ax2.Position(2) = ax1.Position(2) % Ajustar lleugerament cap avall
%xlabel(ax2, 'Latitude of Tangent Point');
ax1.Position(2) = ax1.Position(2) + 0.03;

ax3=subplot(3,4,[9 12]); hold on; grid; legend;
scatter(timea, row1, [], c1, '.', DisplayName=append('Row ',num2str(rown1)))
scatter(timea, row2, [], c2, '.', DisplayName=append('Row ',num2str(rown2)))
scatter(timea, row3, [], c3, '.', DisplayName=append('Row ',num2str(rown3)))
yline(avg1, color = c1, linewidth = 2, HandleVisibility='off')
yline(avg2, color = c2, linewidth = 2, HandleVisibility='off')
yline(avg3, color = c3, linewidth = 2, HandleVisibility='off')
xlim([timea(1) timea(end)])
ylim([0 2.1e14])
ylabel('ph \cdot nm^{-1} \cdot m^{-2} \cdot sr^{-1} \cdot s^{-1}')
%nexttile([1,4]); hold on;
%imagesc(timea, 1:187, keogram_wb, [0 2e14])
%xlim([timea(1) timea(end)])
%ylim([1 187])

text(ax11,'Units', 'Normalized','Position', [0.04, 0.1],'string','a)', color = 'w')
text(ax12,'Units', 'Normalized','Position', [0.04, 0.1],'string','b)', color = 'w')
text(ax13,'Units', 'Normalized','Position', [0.04, 0.1],'string','c)', color = 'w')
text(ax14,'Units', 'Normalized','Position', [0.04, 0.1],'string','d)', color = 'w')
text(ax1 ,'Units', 'Normalized','Position', [0.01, 0.1],'string','e)', color = 'w')
text(ax3 ,'Units', 'Normalized','Position', [0.01, 0.9],'string','f)')

set(findall(gcf,'-property','FontSize'),'FontSize',14)
saveas(fig,savefig)

