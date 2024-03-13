%% Replace time vector in struct with matlab datetime list
load('Weekdata\may1WpeaksSH.mat')

datechar = datetime(may1WpeaksSH.time,'Format','dd/MM HH:mm:ss');  %dd/MM
may1WpeaksSH.time = datechar' ;
save('Weekdata\may1WpeaksSH.mat', "may1WpeaksSH")


%% Remove certain indices
idx = [20,36] ; %indices to remove

load("Weekdata\Aprilstrips\apr3WpeaksNH.mat")
struct = apr3WpeaksNH ;

keys = fieldnames(struct) ;
for i = 1:length(keys)
    key = keys{i};
    field = getfield(struct,key) ; %#ok<*GFLD>
    field(idx) = [] ;
    struct = setfield(struct,key,field) ; %#ok<*SFLD>
    
end
apr3WpeaksNH = struct ;
save('Weekdata\Aprilstrips\apr3WpeaksNH.mat', "apr3WpeaksNH")


%% Combine peak files for a month
addpath("Weekdata/Aprilstrips/") ;
load("apr1WpeaksNH.mat");
load("apr1WpeaksSH.mat") ;
load("apr2WpeaksNH.mat") ;
load("apr2WpeaksSH.mat") ;
load("apr3WpeaksNH.mat") ;
load("apr3WpeaksSH.mat") ;
load("apr4WpeaksNH.mat");
load("apr4WpeaksSH.mat") ;

NH1 = apr1WpeaksNH;
NH2 = apr2WpeaksNH;
NH3 = apr3WpeaksNH;
NH4 = apr4WpeaksNH;

SH1 = apr1WpeaksSH;
SH2 = apr2WpeaksSH;
SH3 = apr3WpeaksSH;
SH4 = apr4WpeaksSH;

keys = fieldnames(NH4) ;
newSH = struct() ;
newNH = struct() ;
for j = 1:length(keys)
    field = keys{j} ;
    newSH.(field) = [SH1.(field), SH2.(field) , SH3.(field) , SH4.(field)] ;
    newNH.(field) = [NH1.(field), NH2.(field) , NH3.(field) , NH4.(field)] ;

end
aprpeaksSH = newSH ;
aprpeaksNH = newNH ;
save('Monthdata\aprpeaksSH.mat', "aprpeaksSH")
save('Monthdata\aprpeaksNH.mat', "aprpeaksNH")
%% Combine aurora strip files for a month
addpath("Weekdata/Marchstrips/") ;
load("mar1Waurorastrips.mat");
load("mar2Waurorastrips.mat") ;
load("mar3Waurorastrips.mat") ;
load("mar4Waurorastrips.mat");

S1 = mar1Waurorastrips;
S2 = mar2Waurorastrips;
S3 = mar3Waurorastrips;
S4 = mar4Waurorastrips;

keys = fieldnames(S2) ;
newS = struct() ;
for j = 1:length(keys)
    field = keys{j} ;
    newS.(field) = [S1.(field), S2.(field) , S3.(field) , S4.(field)] ;

end
maraurorastrips = newS ;
save('Monthdata\maraurorastrips.mat', "maraurorastrips")

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
        kp(end+1) = KPdata(2,day) ;    %#ok<*SAGROW>

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

%% This function splits list of strips to North and South hemisphere depending on TPlat 
function [MLTNH,MlatNH,MLTSH,MlatSH] = splitHemisphere(Mlat, MLT) %#ok<*DEFNU>
    L = length(Mlat) ;
    MLTNH =  [] ;
    MlatNH = [] ;
    MLTSH =  [] ;
    MlatSH = [] ;

    for m = 1:L
        if Mlat(m) > 0
            MLTNH(end+1) = MLT(m) ;
            MlatNH(end+1)= Mlat(m); %#ok<*AGROW>
        end
    
        if Mlat(m) < 0
            MLTSH(end+1) = MLT(m) ;
            MlatSH(end+1) = Mlat(m);
        end
    end  
end