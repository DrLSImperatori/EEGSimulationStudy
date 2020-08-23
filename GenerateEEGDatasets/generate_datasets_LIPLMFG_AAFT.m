function generate_datasets_LIPLMFG_AAFT(ndatasets, numWorkers) %, dataset_string, interaction)
% load head model and some miscallaneous data
load('/home/laurai/matlab/BBCBCode/data/sa')
load('/home/laurai/matlab/BBCBCode/data/miscdata')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saStruct=sa;
tau=14;
freqrange=[0.5 12];
eplength=2;
chanlocs=sa.EEG_clab_electrodes;
date_string=datestr(now,'ddmmyyyy');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataset_string={'true'; 'false'};
for connkind=1:2
% % create directory to store data in
mkdir(['/home/laurai/matlab/EEGSimulationStudy/SimulationResults/Conn_BBCB_Results_LIPLMFG_AAFT2_' date_string '/' dataset_string{connkind}])
end
%number of biological noise source
n_noise_sources = 500;
% sampling frequency
fs = 500;
% number of electrodes
EEG_M = length(sa.EEG_clab_electrodes);

parpool(numWorkers)
parfor idata = 1:ndatasets 
  Struct1=load(['/home/laurai/matlab/BBCB/truth_struct_liplmfg.mat']); 
  truth=Struct1.truth;
  Struct2=load(['/home/laurai/matlab/BBCB/timeseries.mat'], 'sources_int', 'x_h', 'y_h', 'x_i', 'y_i', 'x_r', 'y_r', 'z_r', 'x_l', 'y_l', 'z_l');
  Struct3=load(['/home/laurai/matlab/BBCB/noise_indices.mat']); 

  
for snr =linspace(0.01,1,100)
      snr_index=snr*100;
    % spatial standard deviation of the sources (along cortical manifold) in mm
% highpass filter coefficients
[b_high, a_high] = butter(3, 0.5/fs*2, 'high');
truth.snr=snr;
%truth.sigma_range = [10 30];
% length of recording in sec
truth.len = 120; % include additional epoch in beginning and end of data
% resulting number of samples
N = fs*truth.len;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interactionkind={'linear', 'henon', 'ikeda', 'lorenz_xy', 'lorenz_xz', 'lorenz_yz', 'rossler_xy', 'rossler_xz', 'rossler_yz'}; % 'henon', 'ikeda'
for idx=1:length(interactionkind) 
  
if ~exist(['/home/laurai/matlab/Simulations/Conn_BBCB_Results_LIPLMFG_AAFT2_' date_string '/' dataset_string{2} '/dataset_' char(interactionkind(idx))  '_' num2str(snr_index) '_' num2str(idata) '.mat'], 'file');

    % Load truth.sources_int and truth.sources_nonint.

if strcmp(interactionkind(idx),'linear')==1
    [truth.sources_int] = Struct2.sources_int; %generate_sources_ar_laura(fs, truth.len, truth.bandpass);
    truth.sources_nonint = [AAFT(truth.sources_int(1,:)),AAFT(truth.sources_int(2,:))].';
end

if strcmp(interactionkind(idx),'henon')==1
    truth.sources_int= [Struct2.x_h,Struct2.y_h].'; 
    truth.sources_nonint = [AAFT(truth.sources_int(1,:)),AAFT(truth.sources_int(2,:))].';
end

if strcmp(interactionkind(idx),'ikeda')==1 
    truth.sources_int= [Struct2.x_i,Struct2.y_i].'; 
    truth.sources_nonint = [AAFT(truth.sources_int(1,:)),AAFT(truth.sources_int(2,:))].';
end

if strcmp(interactionkind(idx), 'lorenz_xy')==1
    truth.sources_int= [Struct2.x_l,Struct2.y_l].'; 
    truth.sources_nonint = [AAFT(truth.sources_int(1,:)),AAFT(truth.sources_int(2,:))].';
end

if strcmp(interactionkind(idx), 'lorenz_yz')==1
    truth.sources_int= [Struct2.y_l,Struct2.z_l].'; 
    truth.sources_nonint = [AAFT(truth.sources_int(1,:)),AAFT(truth.sources_int(2,:))].';
end

if strcmp(interactionkind(idx), 'lorenz_xz')==1
    truth.sources_int= [Struct2.x_l,Struct2.z_l].'; 
    truth.sources_nonint = [AAFT(truth.sources_int(1,:)),AAFT(truth.sources_int(2,:))].';
end

if strcmp(interactionkind(idx), 'rossler_xy')==1
    truth.sources_int= [Struct2.x_r,Struct2.y_r].'; 
    truth.sources_nonint = [AAFT(truth.sources_int(1,:)),AAFT(truth.sources_int(2,:))].';
end

if strcmp(interactionkind(idx), 'rossler_yz')==1
    truth.sources_int= [Struct2.y_r,Struct2.z_r].'; 
    truth.sources_nonint = [AAFT(truth.sources_int(1,:)),AAFT(truth.sources_int(2,:))].';
end

if strcmp(interactionkind(idx), 'rossler_xz')==1
    truth.sources_int= [Struct2.x_r,Struct2.z_r].'; 
    truth.sources_nonint = [AAFT(truth.sources_int(1,:)),AAFT(truth.sources_int(2,:))].';
end
 
% sample noise source locations
noise_inds = Struct3.noise_inds_mat(:,idata);   % ceil(size(saStruct.cortex75K.EEG_V_fem_normal, 2)*rand(n_noise_sources, 1));
% generate pseudo-EEG/MEG with interacting sources     
no_signal = norm(truth.source_amp*truth.sources_int, 'fro');
EEG_signal = truth.EEG_field_pat*truth.sources_int;
% generate pseudo-EEG/MEG with non-interacting sources
no_signal_fake = norm(truth.source_amp*truth.sources_nonint, 'fro');
EEG_signal_fake = truth.EEG_field_pat*truth.sources_nonint;

EEG_signal = EEG_signal ./ no_signal;
EEG_signal_fake = EEG_signal_fake ./ no_signal_fake;


% biological noise, mixture of independent pink noise sources 
pn = mkpinknoise(N, n_noise_sources)';
EEG_brain_noise = saStruct.cortex75K.EEG_V_fem_normal(:, noise_inds)*pn;

if idx==1
  truth.bandpass = [8 12]; % bandpass filter coefficients
  [b_band, a_band] = butter(3, truth.bandpass/fs*2);
else
  truth.bandpass = [0.5 12]; % bandpass filter coefficients
  [b_band, a_band] = butter(3, truth.bandpass/fs*2);
end


norm_brain_noise = norm(filtfilt(b_band, a_band, pn')', 'fro');
% normalize noise by amplitude in band of interest
EEG_brain_noise = EEG_brain_noise ./ norm_brain_noise;

EEG_brain_signal_noise = truth.snr*EEG_signal + (1-truth.snr)*EEG_brain_noise;
EEG_brain_signal_noise = EEG_brain_signal_noise ./ norm(EEG_brain_signal_noise, 'fro');


EEG_brain_signal_noise_fake = truth.snr*EEG_signal_fake + (1-truth.snr)*EEG_brain_noise;
EEG_brain_signal_noise_fake = EEG_brain_signal_noise_fake ./ norm(EEG_brain_signal_noise_fake, 'fro');

% white sensor noise
EEG_sensor_noise = randn(EEG_M, N);
EEG_sensor_noise = EEG_sensor_noise ./ norm(EEG_sensor_noise, 'fro');

% overall noise is dominated by biological noise
EEG_data = 0.9*EEG_brain_signal_noise + 0.1*EEG_sensor_noise;
EEG_data_fake = 0.9*EEG_brain_signal_noise_fake + 0.1*EEG_sensor_noise;

% apply high-pass
EEG_data = filtfilt(b_high, a_high, EEG_data')';
EEG_data_fake = filtfilt(b_high, a_high, EEG_data_fake')';


Mats=load('/home/laurai/matlab/BBCB/BBCB_CSD_GH.mat');
G=Mats.G;
H=Mats.H;

CSD_EEG_data=CSD(EEG_data, G, H, 1.0e-5, 1);
CSD_EEG_data_fake=CSD(EEG_data_fake, G, H, 1.0e-5, 1);

EEG=struct();
EEG.srate=fs;
EEG.tr=truth.len/eplength;
EEG.data=reshape(CSD_EEG_data, [size(EEG_data,1), fs*eplength, EEG.tr]);
if idx==1
    [wpli, wsmi]=getwPLIandwSMI_BBCB(EEG,freqrange, tau, chanlocs);
    [wpli_alpha, wsmi_alpha]=getwPLIandwSMI_BBCB(EEG, truth.bandpass, tau, chanlocs);
else
    [wpli, wsmi]=getwPLIandwSMI_BBCB(EEG, freqrange, tau, chanlocs);
end


EEGfake=struct();
EEGfake.srate=fs;
EEGfake.tr=truth.len/eplength;
EEGfake.data=reshape(CSD_EEG_data_fake, [size(EEG_data,1), fs*eplength, EEGfake.tr]);
if idx==1
    [wpli_f, wsmi_f]=getwPLIandwSMI_BBCB(EEGfake, freqrange, tau, chanlocs);
    [wpli_f_alpha, wsmi_f_alpha]=getwPLIandwSMI_BBCB(EEGfake, truth.bandpass, tau, chanlocs);
else
    [wpli_f, wsmi_f]=getwPLIandwSMI_BBCB(EEGfake, freqrange, tau, chanlocs);
end

for connkind=1:2

if connkind==1
    truth.dataset = {'interacting'};
    truth.sources=truth.sources_int;
    if idx==1
        Res=struct();
        Res.wpli=wpli;
        Res.wsmi=wsmi;
        Res.wpli_alpha=wpli_alpha;
        Res.wsmi_alpha=wsmi_alpha; 
        Res.truth=truth;
    savetofile(['/home/laurai/matlab/Simulations/Conn_BBCB_Results_LIPLMFG_AAFT2_' date_string '/' dataset_string{1} '/dataset_' char(interactionkind(idx))  '_' num2str(snr_index) '_' num2str(idata) '.mat'], Res);
    else
        Res=struct();
        Res.wpli=wpli;
        Res.wsmi=wsmi;
        Res.truth=truth;
    savetofile(['/home/laurai/matlab/Simulations/Conn_BBCB_Results_LIPLMFG_AAFT2_' date_string '/' dataset_string{1} '/dataset_' char(interactionkind(idx))  '_' num2str(snr_index) '_' num2str(idata) '.mat'], Res);
    end
    
else
    truth.dataset = {'non-interacting'};
    truth.sources=truth.sources_nonint;
    if idx==1      
        Res=struct();
        Res.wpli_f=wpli_f;
        Res.wsmi_f=wsmi_f;
        Res.wpli_f_alpha=wpli_f_alpha;
        Res.wsmi_f_alpha=wsmi_f_alpha; 
    savetofile(['/home/laurai/matlab/Simulations/Conn_BBCB_Results_LIPLMFG_AAFT2_' date_string '/' dataset_string{2} '/dataset_' char(interactionkind(idx))  '_' num2str(snr_index) '_' num2str(idata) '.mat'], Res);
    else
        Res=struct();
        Res.wpli_f=wpli_f;
        Res.wsmi_f=wsmi_f;
    savetofile(['/home/laurai/matlab/Simulations/Conn_BBCB_Results_LIPLMFG_AAFT2_' date_string '/' dataset_string{2} '/dataset_' char(interactionkind(idx))  '_' num2str(snr_index) '_' num2str(idata) '.mat'], Res);
    end
end
end
end
end
end  
end
end