%plot kp-index for whole period
% From https://kp.gfz-potsdam.de/en/data
clear all
[time, value, status] = getKpindex('2023-02-01', '2023-02-28', 'Kp') ;
days = 28 ;

% Reshape yourValues into a matrix
KPdata = reshape(value, 8, days);

figure(1)
bar(time,value)
title('Plot of the Kp-index vs time')
ylabel('Kp')
xlabel('Time')
grid on
%% Create kp_data matrix for a month 

%Loads kp data from website
[time, value, status] = getKpindex('2023-05-01', '2023-05-31', 'Kp') ;
days = 31 ;
% Reshape into matrix
KPdata = reshape(value, 8, days);
save('KPdataMay','KPdata')
% Different rows  in the matrix correspond to the measured time: 
% 1 = 00 , 2 = 03, 3 =06, 4 = 09, 5 = 12, 6 = 15, 7 = 18, 8 = 21

%% Function to create kp values corresponding to the peaks
clear all
clc

addpath("Weekdata\Maystrips\")
load("may1WpeaksSH.mat")
load("KPdataMay.mat")

timedata = may1WpeaksSH.time ;
kp = [] ;

%loops through list of strips
for i = 1:length(timedata)
    day = timedata(i).Day ;
    hour = timedata(i).Hour ;
    if (hour >= 0) && (hour <= 3)
        kp(end+1) = KPdata(1,day) ;

    elseif (hour > 3) && (hour <= 6)
        kp(end+1) = KPdata(2,day) ;    

    elseif (hour > 6) && (hour <= 9)
        kp(end+1) = KPdata(3,day) ; 

    elseif (hour > 9) && (hour <= 12)
        kp(end+1) = KPdata(4,day) ;
    
    elseif (hour > 12) && (hour <= 15)
        kp(end+1) = KPdata(5,day) ;

    elseif (hour >= 15) && (hour <= 18)
        kp(end+1) = KPdata(6,day) ;

    elseif (hour > 18) && (hour <= 21)
        kp(end+1) = KPdata(7,day) ;

    elseif (hour > 21) && (hour < 24)
        kp(end+1) = KPdata(8,day) ;
    end
end

may1WpeaksSH.kp = kp ;
save('Weekdata\may1WpeaksSH.mat', "may1WpeaksSH")