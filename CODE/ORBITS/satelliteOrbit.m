function [X, V, dT_rel, exclude, status] = satelliteOrbit(prn, Ttr, preciseEph, settings, exclude, status)
% Calculate satellite position from precise satellite orbits (*.sp3).
%
% INPUT:
% 	prn         satellite vehicle number
% 	Ttr         transmission time, sow (GPS time)
% 	preciseEph  precise orbits for this GNSS (from sp3)
%   settings	struct with settings from GUI
%   exclude     true, if satellite has to be excluded
%   status      satellite status
%
% OUTPUT:
%   X           satellite position [m]
%   V           satellite velocity [m/s]
%   dT_rel      relativistic correction, calculated from nav message [s]
%   exclude     true, if satellite has to be excluded
%
% Revision:
%   2025/02/12, MFWG: added missing conversion [mm/s] to [m/s] (corr2brdc)
%   2025/02/12, MFWG: call of function SatPos_brdc changed
%   2025/09/13, MFWG: change structure
%   2025/10/02, MFWG: move brdc+corr2brdc to satelliteOrbitBrdc.m
%
% This function belongs to raPPPid, Copyright (c) 2023, M.F. Glaner
% *************************************************************************


%% Preparations
sv = mod(prn, 100);
bool_print = ~settings.INPUT.bool_parfor;
dT_rel = 0;


%% Calculation of Satellite Position
[X, V, exclude, status] ...
    = prec_satpos(preciseEph, prn, sv, Ttr, exclude, status, bool_print);

