function wpli=getwPLI(EEG_data,chanlocs, tr, freqrange, srate, trialsep)

labels=1:trialsep*size(EEG_data,2):(trialsep*size(EEG_data,2)*tr);
data.sampleinfo=[labels; labels+(size(EEG_data,2)-1)]';
   %%
for i=1:1:tr %(size(datavr,3)) 
    %count=floor((size(datavr,3)/2))+i
    data.time{i}=(1/srate:1/srate:((1/srate)*size(EEG_data,2)));
    data.trial{i}=EEG_data(:,:,i);
end
    data.label= chanlocs;    % cell-array containing strings, Nchan*1
    %channel = ft_channelselection({'all', '-Cz'}, data.label);

    cfg1           = [];
    cfg1.method    = 'mtmfft';
    cfg1.output    = 'powandcsd';
    cfg1.taper=      'hanning';
    cfg1.foilim = freqrange; %theta [1 4], alpha [8 12], beta [18 25], gamma [35 45]
    cfg1.keeptrials = 'yes';
    cfg1.feedback = 'none';
    
       
    cfg1.pad ='nextpow2';
    freq          = ft_freqanalysis(cfg1, data);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   [wpli, ~, ~] = ft_connectivity_wpli(freq.crsspctrm, 'dojack', 0,'feedback', 'none', 'debiased', 0);
%    wPLI=mean(abs(wpli),2);
%    wPLImatrix=transforming_sequence_into_matrix(abs(wPLI),length(chanlocs));
   
end
