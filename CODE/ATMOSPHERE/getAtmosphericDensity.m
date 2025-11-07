function rhoAtmo = getAtmosphericDensity(alt_m)
% Extracts the Earth's atmospheric density for a given altitude
% 
% INPUT:
%   altitude        altitude [m]
% OUTPUT:
%	rhoAtmo         Earth's atmospheric density  [kg/m^3] 
%
% Revision:
%   ...
%
% Created by Hoor Bano
% 
% This function belongs to raPPPid, Copyright (c) 2025, M.F. Wareyka-Glaner
% *************************************************************************

h = alt_m/1000;  % km

% Altitude [km]     Density [kg/m^3]  (order-of-magnitude, MSIS-like)
alt_km = [200  300    400     500     600     700     800     900     1000];
rho_tab = [1.9e-9 2.2e-11 3.5e-12 1.0e-12 3.0e-13 1.0e-13 5.0e-14 2.0e-14 1.0e-14];

% Clamp & log-linear interpolate
if h <= alt_km(1)
    rhoAtmo = rho_tab(1);
elseif h >= alt_km(end)
    rhoAtmo = rho_tab(end);
else
    rhoAtmo = exp(interp1(alt_km, log(rho_tab), h, 'linear'));
end