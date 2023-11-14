%Mean values
mean(maxalt15feb)
mean(maxlat15feb)
mean(maxlon15feb)
%% Plot a keogram and a corresponding position plot below (ex altitude)
clims = [min(min(keogramgrad)) max(max(keogramgrad))];
figure(1)
subplot(2,1,1) ;
imagesc(keogramgrad,clims)
ax1 = gca;
ax1.YDir = 'normal' ;
title('Keogram 3 Mars orbit 10')

dateStrings = time3mars; 
dateTimes = datetime(dateStrings, 'Format','dd/MM HH:mm:ss');  %dd/MM

subplot(2,1,2)
plot(alt3mars, '.') ;
title('Altitude 15 feb')
ylabel('Altitude (km)')
grid on
ylim([80,115])
xtickangle(45);
linkaxes([ax1 gca], 'x');  % Link x-axes of the two subplots

%%  Plotting the pixel values of different row of a keogram.
close all 
figure(1)
image = keogram15feborb11 ;
plot(image(140,:),'DisplayName','140')
hold on
plot(image(141,:),'DisplayName','141')
hold on
plot(image(150,:),'DisplayName','150')
hold on
plot(image(151,:),'DisplayName','151')
hold on
plot(image(152,:), 'DisplayName','152')
hold on
plot(image(153,:),'DisplayName','153')
hold on
plot(image(154,:),'DisplayName','154')
hold on
plot(image(155,:), 'DisplayName','155')
hold on
plot(image(160,:),'DisplayName','160')
title('Plot showing Keogram rows')

%% Plot keogram image and some of the polynomial regression rows
figure(1)
subplot(2,1,1)
clims = [min(min(keogramgrad)) max(max(keogramgrad))];
imagesc(keogramgrad,clims)
ax1 = gca;
ax1.YDir = 'normal' ;
title('Keogram 26 feb orbit 3')

subplot(2,1,2)
plot(keogram(180,:),'DisplayName','Row 180-diff')
hold on
plot(y, 'DisplayName','Row 140')
hold on
plot(y_reg, 'DisplayName','140-linearzied')
hold on
plot((y-y_reg), 'DisplayName','140-diff')
plot(keogram(150,:),'DisplayName','150')
hold on
plot(keogram(155,:), 'DisplayName','155')
hold on
plot(keogram(160,:),'DisplayName','160')

title('Plot showing Keogram rows')
linkaxes([ax1 gca], 'x');  % Link x-axes of the two subplots
grid on
legend
%% Check means and maximum values of specific keogram columns
 
col = keogram(151:end,106) ;
coltest = keogram(161:end,109);

sum(col)/length(col) 
max(col) 
sum(coltest)/length(coltest) 
max(coltest) 