%author: Nicole Bienert
%date: 6/28/21
%Purpose: to read E312 data in different ways depending on formatting. 
%The E312 data with GPS stamping has the GPS info at the end of 8s of
%recording. 
%outputs lat and long in decimal degrees


function [data, gpsData] = readE312Data(directory, dataType, gpsYN);
gpsData = struct('string', [], 'valid', [], 'lat', [], 'lon', []);
gpsData.string = ['NoGPS'];
fID = fopen(directory); 
data = single(fread(fID, dataType)); 

if gpsYN == 1
    %100 bytes long=800 bits of GPS data 
    status = fseek(fID, -100, 'eof'); %moves to 100 bytes from end of file
    gpsstr = fread(fID, 'ubit8');     %read the last 100 bytes
    data = data(1:length(data)-50);   %The last 50 shorts are the GPS data
    
    %turn GPS data into a string
    gpsstr = char(gpsstr);
    gpsstr = string(reshape(gpsstr, [100,1])');
    gpsData.string = gpsstr;
    
    gpsstr = regexp(gpsData.string, ',', 'split');

    nan_check = true;
    
    %determine whether gpsstr is valid
    if (length(gpsstr) == 16 || length(gpsstr) == 15)    
        lat_check=strfind(gpsstr(3), 'nan');
        lon_check=strfind(gpsstr(5), 'nan');
        % nan_check = true;
        if length(lat_check) ~= 0 || length(lon_check) ~= 0
            nan_check=false;
        end
    end
    
    %GPS string has ~16 elements. If it doesn't, then this file was probably
    %stopped by the user and wasn't gps stamped.
	if (length(gpsstr) == 16 || length(gpsstr) == 15) && ...
            nan_check && (strcmp(gpsstr(4), 'N') ||  ...
            strcmp(gpsstr(4), 'S') ) & (strcmp(gpsstr(6), 'E') || ...
            strcmp(gpsstr(6), 'W') )
        
        gpsData.valid = 1;
        
        %find lat
        b = regexp(gpsstr(3), '\.', 'split');
        c = char(b(1));
        d = c(end-1:end);
        d = [d, '.', char(b(2))];
        gpsData.lat = str2double(d)/60;
        e = c(1:end-2);
        gpsData.lat = str2double(e)+gpsData.lat;
        if strcmp(gpsstr(4), 'S')
            gpsData.lat = gpsData.lat*-1;
        end
        
        %find lon
        b = regexp(gpsstr(5), '\.', 'split');
        c = char(b(1));
        d = c(end-1:end);
        d = [d, '.', char(b(2))];
        gpsData.lon = str2double(d)/60;
        e = c(1:end-2);
        gpsData.lon = str2double(e)+gpsData.lon;
        if strcmp(gpsstr(6), 'W')
            gpsData.lon = gpsData.lon*-1;
        end    
        
        %get GMT time
        b = regexp(gpsstr(2), '\.', 'split');
        c = char(gpsstr(2));
        gpsData.GMThr = c(1:2);
        gpsData.GMTmin = c(3:4);
        gpsData.GMTsec = c(5:end);
        gpsData.GMTtime = str2num(gpsData.GMThr)*3600+ ...
                          str2num(gpsData.GMTmin)*60+str2num(gpsData.GMTsec);
    else
        gpsData.valid = 0;
    end
    
end
fclose(fID);
%data is complex
data = complex(data(1:2:length(data)),data(2:2:length(data))); 

