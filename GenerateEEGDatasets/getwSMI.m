function [wsmi]=getwSMI(EEG_data, tau, srate, lowfreq)
cfg=struct();
cfg.chan_sel = 1:(size(EEG_data,1));  % compute for all pairs of channels
cfg.data_sel = 1:(size(EEG_data,2)); % compute using all samples
cfg.taus     = tau; % compute for taus 1 2 4 8
cfg.kernel   = 3; % kernel = 3 (3 samples per symbol)
cfg.sf       = srate;  % sampling frequency
cfg.lowfreq  = lowfreq;
cfg.over_trials       = 0;  
[~, ~, ~, wsmi] = smi_and_wsmi_zeropad(EEG_data, cfg);

end

