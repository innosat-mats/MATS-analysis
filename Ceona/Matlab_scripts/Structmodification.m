%% Functions used to add point afterwards in structure
load("Weekdata\Maystrips\may1Waurorastrips.mat")
load("Weekdata\Maystrips\may1WpeaksSH.mat")
keys = fieldnames(may1Waurorastrips) ;
tostruct = may1WpeaksSH ;
fromstruct = may1Waurorastrips ;
fromidx = 946 ;
toidx = 26; %index which is before the location you want to place the new

for i = 1:length(keys)
    key = keys{i};
    fromfield = getfield(fromstruct,key) ; %gets the field list %#ok<*GFLD>
    val = fromfield(fromidx) ; %gets the specified value in field list 
    tofield = getfield(tostruct,key) ; %gets the field list of new structure
    newfieldlist = [tofield(1:toidx),val,tofield(toidx+1:end)] ; %inserts value in specified location ; 
    tostruct = setfield(tostruct,key,newfieldlist) ; %#ok<*SFLD>   
end

may1WpeaksSH = tostruct ;
save('Weekdata\Maystrips\may1WpeaksSH.mat', "may1WpeaksSH")

%% Replace time vector in struct with matlab datetime list
load('Weekdata\Marchstrips\mar3WpeaksSH.mat')

datechar = datetime(mar3WpeaksSH.time,'Format','dd/MM HH:mm:ss');  %dd/MM
mar3WpeaksSH.time = datechar' ;
save('Weekdata\mar3WpeaksSH.mat', "mar3WpeaksSH")

%% Remove certain indices
idx = [26] ; %ex indices to remove

load("Weekdata\Februarystrips\feb3Wpeaks.mat")
struct = feb3Wpeaks;

keys = fieldnames(struct) ;
for i = 1:length(keys)
    key = keys{i};                                                                                                          
    field = getfield(struct,key) ; %#ok<*GFLD>
    field(idx) = [] ;
    struct = setfield(struct,key,field) ; %#ok<*SFLD>
    
end
feb3Wpeaks = struct ;
save('Weekdata\Februarystrips\feb3Wpeaks.mat', "feb3Wpeaks")


%% Combine peak files for a month but also separately for the hemispheres.
addpath("Weekdata/Februarystrips/") ;
%load("apr1WpeaksNH.mat");
%load("apr1WpeaksSH.mat") ;
load("feb2WpeaksNH.mat") ;
load("feb2WpeaksSH.mat") ;
load("feb3WpeaksNH.mat") ;
load("feb3WpeaksSH.mat") ;
load("feb4WpeaksNH.mat");
load("feb4WpeaksSH.mat") ;

%NH1 = apr1WpeaksNH;
NH2 = feb2WpeaksNH;
NH3 = feb3WpeaksNH;
NH4 = feb4WpeaksNH;

%SH1 = apr1WpeaksSH;
SH2 = feb2WpeaksSH;
SH3 = feb3WpeaksSH;
SH4 = feb4WpeaksSH;

keys = fieldnames(NH4) ;
newSH = struct() ;
newNH = struct() ;
newH = struct() ;
for j = 1:length(keys)
    field = keys{j} ;
    newSH.(field) = [SH2.(field) , SH3.(field) , SH4.(field)] ; %SH1.(field), Remove SH1 and NH1 when February
    newNH.(field) = [NH2.(field) , NH3.(field) , NH4.(field)] ;  %NH1.(field),
    newH.(field) = [SH2.(field),NH2.(field) , SH3.(field),NH3.(field) , SH4.(field),  NH4.(field)] ; %SH1.(field),NH1.(field), 
end
febpeaksSH = newSH ;
febpeaksNH = newNH ;
febpeaks = newH ;
save('Monthdata\Februarymonth\febpeaks.mat', "febpeaks")
save('Monthdata\Februarymonth\febpeaksSH.mat', "febpeaksSH")
save('Monthdata\Februarymonth\febpeaksNH.mat', "febpeaksNH")
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
load("may1Wpeaks.mat")
%load("KPdataMarch.mat")

timedata = may1Wpeaks.time ;
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

may1Wpeaks.kp = kp ;
save('Weekdata\Maystrips\may1Wpeaks.mat', "may1Wpeaks")

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