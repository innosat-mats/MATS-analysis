close all
%%% Plots altitude vs time
figure(1)
dateStrings = feb3Waurorastrips.time; 
dateTimes = datetime(dateStrings, 'Format', 'dd/MM HH:mm:ss');
plot(dateTimes,feb3Waurorastrips.alt,'.')
title('15 feb, limit, non-uniform time spacing of x-axis')
ylabel('Altitude (km)')
grid on
ylim([100,116]) ;
xtickangle(45);

% Select a subset of date strings to display as x-axis ticks, for
% non-uniform x-axis
numTicks = 12;  % Number of ticks to display
selectedIndices = round(linspace(1, length(dateTimes), numTicks));
xticks(dateTimes(selectedIndices));
xticklabels(datestr(dateTimes(selectedIndices), 'HH:MM:SS'));

figure(2)
plot(feb3Waurorastrips.alt,'.')
title('15 feb, limit, uniform spacing of x-axis')
ylabel('Altitude (km)')
grid on
ylim([100,116]) ;
xtickangle(45);

%% Extra scripts and functions
% This function splits strips to North and South hemisphere depending on TPlat 
function [MLTNH,MlatNH,MLTSH,MlatSH] = splitHemisphere(Mlat, MLT)
    L = length(Mlat) ;
    MLTNH =  [] ;
    MlatNH = [] ;
    MLTSH =  [] ;
    MlatSH = [] ;

    for m = 1:L
        if Mlat(m) > 0
            MLTNH(end+1) = MLT(m) ;
            MlatNH(end+1)= Mlat(m);
        end
    
        if Mlat(m) < 0
            MLTSH(end+1) = MLT(m) ;
            MlatSH(end+1) = Mlat(m);
        end
    end  
end