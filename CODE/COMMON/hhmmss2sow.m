function sow = hhmmss2sow(hh, mm, ss, startdate)
% Convert hours, minutes and seconds to seconds of week (GPS time).
% 
% INPUT:
% 	hh          hours
%   mm          minutes
%   ss          seconds
%   startdate   1x3, [year month day], startdate of observation file
% OUTPUT:
%	sow         seconds of week (GPS time)
%
% Revision:
% 2025/08/14, MFWG: switch to cal2gpstime
%
% This function belongs to raPPPid, Copyright (c) 2025, M.F. Wareyka-Glaner
% *************************************************************************

% convert to julian date and then seconds of week (GPS time)
[~, sow] = cal2gpstime([startdate(1:3) hh mm ss]);