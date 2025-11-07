function [num_freqs, proc_freqs, n_gnss_freqs] = CountProcessedFrequencies(settings)
% This function determines the number of input and processed frequencies.
% Through building a LC the number of processed frequencies might be lower
% than the number of the input frequencies. 
% 
% INPUT:
%   settings        struct, processing settings from GUI
% OUTPUT:
%	proc_freqs      number of processed frequencies
%   num_freqs       number of input frequencies
%   n_gnss_freqs    number of input frequencies for each GNSS
%
% Revision:
%   ...
%
% This function belongs to raPPPid, Copyright (c) 2023, M.F. Glaner
% *************************************************************************


% number of input frequencies (1 or 2 or 3) 
num_freqs = max([ ... 
    settings.INPUT.use_GPS*numel(settings.INPUT.gps_freq(~strcmpi(settings.INPUT.gps_freq,'OFF'))), ...
    settings.INPUT.use_GLO*numel(settings.INPUT.glo_freq(~strcmpi(settings.INPUT.glo_freq,'OFF'))), ...
    settings.INPUT.use_GAL*numel(settings.INPUT.gal_freq(~strcmpi(settings.INPUT.gal_freq,'OFF'))), ...
    settings.INPUT.use_BDS*numel(settings.INPUT.bds_freq(~strcmpi(settings.INPUT.bds_freq,'OFF')))      ]);

% number of processed frequencies
proc_freqs = num_freqs;   
if strcmpi(settings.IONO.model,'2-Frequency-IF-LCs')
    % through building LC one frequency less is processed and the number
    % of input and processed frequencies is different
    proc_freqs = proc_freqs - 1;              
elseif strcmpi(settings.IONO.model,'3-Frequency-IF-LC')  
    proc_freqs = 1;
end  

% number of input frequencies for each GNSS
n_gnss_freqs = zeros(5,1);
n_gnss_freqs(1) = sum(settings.INPUT.gps_freq_idx  <= DEF.freq_GPS(end)) * settings.INPUT.use_GPS;
n_gnss_freqs(2) = sum(settings.INPUT.glo_freq_idx  <= DEF.freq_GLO(end)) * settings.INPUT.use_GLO;
n_gnss_freqs(3) = sum(settings.INPUT.gal_freq_idx  <= DEF.freq_GAL(end)) * settings.INPUT.use_GAL;
n_gnss_freqs(4) = sum(settings.INPUT.bds_freq_idx  <= DEF.freq_BDS(end)) * settings.INPUT.use_BDS;
n_gnss_freqs(5) = sum(settings.INPUT.qzss_freq_idx <= DEF.freq_QZSS(end))* settings.INPUT.use_QZSS;