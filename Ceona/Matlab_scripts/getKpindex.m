function [time, value, st] = getKpindex(starttime, endtime, index, status)
% getKpindex.m
% ===================================
% GFZ German Research Centre for Geosciences (CC BY 4.0)
% Author I. Wehner
% last modified on 25 May 2022
% -----------------------------------
% download 'Kp', 'ap', 'Ap', 'Cp', 'C9', 'Hp30', 'Hp60', 'ap30', 'ap60', 'SN', 'Fobs' or 'Fadj' index data from kp.gfz-potsdam.de
% date format for start and end is 'yyyy-mm-dd' or 'yyyy-mm-ddTHH:MM:SSZ'
% optional 'def' parameter to get only definitve values (only available for 'Kp', 'ap', 'Ap', 'Cp', 'C9', 'SN')
% Hpo indices ('Hp30', 'Hp60', 'ap30', 'ap60'), Fobs and Fadj do not have the status info
% example: [time, index, status] = getKpindex('2021-09-29', '2021-10-01','Ap','def');
% example: [time, index, status] = getKpindex('2021-09-29T00:00:00Z', '2021-10-01T12:00:00Z','Kp');

time=0; 
value=0; 
st=0;

% set date format
if strlength(starttime) == 10 && strlength(endtime) == 10
    starttime = strcat(starttime, 'T00:00:00Z');
    endtime = strcat(endtime, 'T23:59:00Z');
end

try 
    d1 = datetime(starttime,'InputFormat','yyyy-MM-dd''T''HH:mm:ssZ','TimeZone','UTC');
    d2 = datetime(endtime,'InputFormat','yyyy-MM-dd''T''HH:mm:ssZ','TimeZone','UTC');
    if d1>d2
        error('Start time must be before or equal to end time %s -%s',starttime,endtime);
    end
    if ~any(strcmp({'Kp','ap','Ap','Cp','C9','Hp30','Hp60','ap30','ap60','SN','Fobs','Fadj'},index))
        error('Wrong index parameter! \nAllowed are only the string parameter: Kp, ap, Ap, Cp, C9, Hp30, Hp60, ap30, ap60, SN, Fobs, Fadj');
    end

    % create url
    time_string = strcat("start=", starttime, "&end=", endtime);
    url = strcat('https://kp.gfz-potsdam.de/app/json/?',time_string,'&index=',index);

    if exist('status','var') && strcmp(status,'def')
        url = strcat(url,'&status=def');
    elseif exist('status','var') && ~strcmp(status,'def')
        error('Wrong option parameter: %s \nAllowed is only the string parameter: def', status);
    end
    
    tmp=webread(url);
    try
        time=datetime(tmp.datetime,'TimeZone','utc','InputFormat','yyyy-MM-dd''T''HH:mm:ssZ');
        value=tmp.(index);
        if (~strcmp(index,'Hp60')) && (~strcmp(index,'Hp30')) && (~strcmp(index,'ap60')) && (~strcmp(index,'ap30')) && (~strcmp(index,'Fobs')) && (~strcmp(index,'Fadj'))
            st=tmp.status;
        end
    catch 
        warning('No data availble for your selection!')
    end
    clear tmp;

catch ME
    switch ME.identifier
        case 'MATLAB:datetime:ParseErr'
            error('Wrong datetime string: %s - %s \nDatetime strings must be in format yyyy-mm-dd or yyyy-mm-ddTHH:MM:SSZ',starttime,endtime);
        case 'MATLAB:webservices:HTTP404StatusCodeError'
            error('Connection Error!!! Can not reach %s',url)
        otherwise
            ME.identifier
            rethrow(ME)
    end
end
