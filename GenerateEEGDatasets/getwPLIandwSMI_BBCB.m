function [wpli, wsmi]=getwPLIandwSMI_BBCB(EEG, freqrange, tau, chanlocs)

wpli=getwPLI(EEG.data,chanlocs, EEG.tr, freqrange, EEG.srate, 1);

[wsmi_cell]=getwSMI(EEG.data, tau, EEG.srate, freqrange(1));
wsmi=cell2mat(wsmi_cell);
wsmi=mean(wsmi,3);
                   
end
