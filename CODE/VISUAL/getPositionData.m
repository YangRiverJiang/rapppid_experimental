function [pos_cart, pos_UTM] = getPositionData(storeData, obs, curr_label, MultiPlot)
% Used in MultiPlot.m and StationResultPlot.m
% Extract position data from storeData depending on selected solution to
% plot
%
% This function belongs to raPPPid, Copyright (c) 2023, M.F. Glaner
% *************************************************************************

if MultiPlot.float
    pos_cart = storeData.param(:,1:3);
    pos_UTM  = storeData.posFloat_utm;
elseif MultiPlot.fixed
    try
		try; pos_cart = storeData.param_fix(:,1:3); catch; pos_cart = storeData.xyz_fix(:,1:3); end
        pos_UTM  = storeData.posFixed_utm;
    catch
        pos_cart = NaN(1,3); pos_UTM = NaN(1,3);
        jd = cal2jd_GT(obs.startdate(1), obs.startdate(2), obs.startdate(3));
        [doy, year] = jd2doy_GT(jd);
        errordlg({'No fixed solution for:', curr_label}, ...
            sprintf('%s %04.0f/%03.0f', obs.stationname, year, floor(doy)));
        return
    end
end