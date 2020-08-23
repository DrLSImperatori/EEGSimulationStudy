function fildata=filterBandPass(EEG_data, srate, freqrange)
        data=EEG_data;
        eplength=size(data,2);
        ntrials = size(data,3);
        
        flippeddata=fliplr(data);
        
        fdata=NaN(size(data,1), eplength*2, ntrials);
        fdata(:,((eplength/2)+1):(3*eplength/2), :)=data;
        fdata(:,(1:(eplength/2)), :)=flippeddata(:,((eplength/2)+1):eplength,:);
        fdata(:,(3*eplength/2)+1:(2*eplength), :)=flippeddata(:,1:(eplength/2),:);

        fildata =NaN(size(fdata));
        for tr=1:ntrials
        fildata(:,:,tr) = ft_preproc_bandpassfilter(squeeze(fdata(:,:,tr)), srate, freqrange); 
        end
        fildata(:,1:(eplength/2), :)=[]; fildata(:,(eplength)+1:(3*eplength/2),:)=[];    
end
   
