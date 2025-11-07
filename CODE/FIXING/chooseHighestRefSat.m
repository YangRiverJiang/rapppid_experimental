function refSat = chooseHighestRefSat(Epoch, sats_gnss, elev_gnss, gnss, settings, bool_fixable)
% Function to find or change reference satellite for specific GNSS
%
% INPUT:
%   Epoch           struct, epoch-specific data for current epoch
%   sats_gnss       satellite numbers of current GNSS
%   elev_gps        [°], elevation of satellites of current GNSS
%   gnss            boolean vector for current GNSS
%   settings        struct, processing settings from GUI
%   bool_fixable    boolean, consider fixability of satellites
% OUTPUT:
%   refSat          selected reference satellite
%
% Revision:
%   2025/09/13, MFWG: check for full set of observations instead of
%                     preferring 3-frequency satellites
% 
% This function belongs to raPPPid, Copyright (c) 2023, M.F. Glaner
% *************************************************************************



%% Some conditions for reference satellite

if settings.AMBFIX.bool_AMBFIX          % check if integer ambiguity fixing is enabled
    
    if strcmp(settings.IONO.model, '2-Frequency-IF-LCs')	% PPPAR and IF
        % prefere satellites which have WL or NL fixed
        WL_fix = ~isnan(Epoch.WL_12(sats_gnss));
        NL_fix = ~isnan(Epoch.NL_12(sats_gnss));
        elev_gnss(WL_fix) = elev_gnss(WL_fix) + 90;         % increase their elevation
        elev_gnss(NL_fix) = elev_gnss(NL_fix) + 90;
        
        if settings.INPUT.proc_freqs == 2                   % PPPAR and 2xIF ?
            % prefere satellites which have EW or EN fixed
            EW_fix = ~isnan(Epoch.WL_23(sats_gnss));
            EN_fix = ~isnan(Epoch.NL_23(sats_gnss));
            elev_gnss(EW_fix) = abs(elev_gnss(EW_fix)) + 90;    % increase their elevation
            elev_gnss(EN_fix) = abs(elev_gnss(EN_fix)) + 90;
        end
    end
    
    if bool_fixable
        % reduce elevation of satellites which are not fixable
        unfixable = any(~Epoch.fixable(gnss, :), 2);
        elev_gnss(unfixable) = elev_gnss(unfixable) - 90;
    end
end

% prefere satellites with a full set of observations -> three frequency
% satellites are prefered
bool_all = check_all_obs(Epoch, settings.INPUT.n_gnss_freqs, settings.INPUT.use_GPS, settings.INPUT.use_GLO, settings.INPUT.use_GAL, settings.INPUT.use_BDS, settings.INPUT.use_QZSS);
missing_obs = ~bool_all(gnss);
elev_gnss = elev_gnss - 180 * missing_obs;



%% Find suitable reference satellite

if max(elev_gnss) > 0            
    % at least one satellite has "positive" elevation, take highest ascending satellite
    refSat = sats_gnss(elev_gnss == max(elev_gnss));
else
    % take "least lowest" satellite
    refSat = sats_gnss(elev_gnss == min(elev_gnss));
end

refSat = refSat(1);     % to be on the safe side


% Epoch.refSatGPS/GLO/GAL/BDS/QZS_idx is handled in change2refSat_IF/_DCM.m


