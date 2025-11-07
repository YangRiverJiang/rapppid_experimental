function [obs_intv] = extractObsInterval(path_file)
% Extracts the observation interval from a RINEX observation file which has
% not observation interval information in the header
%
% INPUT:
%   path_file       string, path of RINEX observation file
% OUTPUT:
%   obs_intv        observation interval
%
% Revision:
%  	2025/02/18, MFWG: round interval to 3 fractional digits
%   2025/08/14, MFWG: switch to cal2gpstime
%
% This function belongs to raPPPid, Copyright (c) 2023, M.F. Glaner
% *************************************************************************

obs_intv = 0;
gps_time_1 = [];
gps_time_2 = [];
epochheader_v2 = '';
fid = fopen(path_file,'rt');         % open observation-file

% loop over header
while 1
    line = fgetl(fid);          % get next line
    if (line == -1)
        break                   % end of file reached
    end
    if contains(line,'END OF HEADER')
        break                   % end of header reached
    end
    if contains(line,'RINEX VERSION / TYPE')
        version = str2double(line(6));
    end
end

% loop to get time-stamp of the first two epochs
while 1
    line = fgetl(fid);          % get next line
    if version == 2 && contains(line, epochheader_v2)
        lvalues = textscan(line,'%f %f %f %f %f %f %d %2d%s','delimiter',',');
        epochheader_v2 = line(1:12);
        % convert date into gps-time [sow]
        [~, gps_time] = cal2gpstime([lvalues{1}, lvalues{2}, lvalues{3}, lvalues{4}, lvalues{5}, lvalues{6}]);
        % save and calculate observation interval
        if isempty(gps_time_1)
            gps_time_1 = gps_time;
        else
            gps_time_2 = gps_time;
            obs_intv = gps_time_2 - gps_time_1;
            break
        end
    end
    
    if version >= 3 && contains(line, '> ')
        lvalues = textscan(line,'%*c %f %f %f %f %f %f %d %2d %f');
        % convert date into gps-time [sow]
        [~, gps_time] = cal2gpstime([lvalues{1}, lvalues{2}, lvalues{3}, lvalues{4}, lvalues{5}, lvalues{6}]);
        % save and calculate observation interval
        if isempty(gps_time_1)
            gps_time_1 = gps_time;
        else
            gps_time_2 = gps_time;
            obs_intv = gps_time_2 - gps_time_1;
            break
        end
    end
    
end


% round to 3 fractional digits (according to RINEX specification and header)
obs_intv = round(obs_intv*1e3)/1e3;        

% %% print result of this function
% [~, obs_filename, ext] = fileparts(path_file);
% if obs_intv == 0
%     errordlg({ [obs_filename, ext],...
%         'No observation interval could be detected. Default: 1 sec'},...
%         'Error');
%     obs_intv = 1;
% else
%     msgbox({[obs_filename ext ':'],...
%         'No observation interval in the Rinex header. Observation interval was extracted from the first two epochs:', ...
%         [sprintf('%.1f', obs_intv) ' seconds']}, ...
%         'Information', 'help')
% end




