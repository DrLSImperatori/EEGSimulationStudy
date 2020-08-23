%% Determining Topographic Accuracy 100x100 Idea
close all;
clc

addpath(genpath('/home/laurai/matlab/QTWriter/'))
addpath('/home/laurai/matlab/BBCB/PlottingScripts/')
addpath('/home/laurai/matlab/fullfig/')
addpath(genpath('/home/laurai/matlab/notBoxPlot/'))
interactionkind={'linear', 'henon', 'ikeda', 'lorenz_xy', 'lorenz_xz', 'lorenz_yz', 'rossler_xy', 'rossler_xz', 'rossler_yz'}; % 'henon', 'ikeda'
datasets=100;
sourceloc={'LIPLMFG', 'RIPLMFG'};
interactionkind_des={'linear', 'H\''{e}non', 'Ikeda','Lorenz (x,y)', 'Lorenz (x,z)', 'Lorenz (y,z)', 'R{\"o}ssler (x,y)', 'R{\"o}ssler (x,z)', 'R{\"o}ssler (y,z)'} % 'henon', 'ikeda'

for idx=1:length(interactionkind)  


wSMIMeanSL=zeros(2, 100, 108, 108);
wSMISurrSL=zeros(2, 100, 108, 108);
LColl=zeros(100, 100, 2);
RColl=zeros(100, 100, 2);
LCollF=zeros(100, 100, 2);
RCollF=zeros(100, 100, 2);

for sl=1:2    
% 

wSMIMean=zeros(100, 108, 108);
wSMISurr=zeros(100, 108, 108);

parfor snr_index=1:100
BarChartCollLeft=zeros(datasets,2);
BarChartCollRight=zeros(datasets,2);
BarChartCollLFake=zeros(datasets,2);
BarChartCollRFake=zeros(datasets,2);

wSMI=zeros(108, 108, 100);
wSMIF=zeros(108, 108, 100);

for idata=1:datasets

Struct1=load(['/home/laurai/matlab/Simulations/Conn_BBCB_Results_' sourceloc{sl} '_shuffle2_09032019/true/dataset_' char(interactionkind(idx))  '_' num2str(snr_index) '_' num2str(idata) '.mat']); 

wsmi=Struct1.data.wsmi;
wSMI(:,:, idata)=wsmi;  

Struct2=load(['/home/laurai/matlab/Simulations/Conn_BBCB_Results_' sourceloc{sl} '_shuffle2_09032019/false/dataset_' char(interactionkind(idx))  '_' num2str(snr_index) '_' num2str(idata) '.mat']); 
wsmif=Struct2.data.wsmi_f;
wSMIF(:, :, idata)=wsmif;
wsmif(wsmif==0)=[];

if sl==1
l_wsmi=wsmi(51, 86);
wsmi(wsmi==0)=[];
end_point_wsmi_l=prctile(wsmi, 95);
BarChartCollLeft(idata,1)=l_wsmi;
BarChartCollLeft(idata,2)=end_point_wsmi_l;
BarChartCollLFake(idata,1)=median(wsmi);
BarChartCollLFake(idata,2)=median(wsmif);
elseif sl==2
r_wsmi=wsmi(51, 87);
wsmi(wsmi==0)=[];
end_point_wsmi_r=prctile(wsmi, 95);
BarChartCollRight(idata,1)=r_wsmi;
BarChartCollRight(idata,2)=end_point_wsmi_r;
BarChartCollRFake(idata,1)=median(wsmi);
BarChartCollRFake(idata,2)=median(wsmif);
end

end
wSMIMean(snr_index,:, :)=mean(wSMI, 3);
wSMISurr(snr_index, :, :)=mean(wSMIF, 3);

if sl==1
LColl(snr_index,:,:)=BarChartCollLeft;
LCollF(snr_index,:,:)=BarChartCollLFake;
elseif sl==2
RColl(snr_index,:,:)=BarChartCollRight;
RCollF(snr_index,:,:)=BarChartCollRFake;
end

end

wSMIMeanSL(sl, :, :, :)=wSMIMean;
wSMISurrSL(sl, :, :, :)=wSMISurr;
end


movObj = QTWriter(['/home/laurai/matlab/Simulations/wSMI_' interactionkind{idx} '.mov']);
t0=1;
% Create an animation
hf = fullfig; 
set(hf,'DoubleBuffer','on');
set(gca,'nextplot','replacechildren');    

for snr_index=1:100
suptitle({['wSMI for ' interactionkind_des{idx} ' dynamics at SNR = ' num2str(snr_index/100)]; ' '})  
LwSMI=squeeze(wSMIMeanSL(1, snr_index, :, :));
RwSMI=squeeze(wSMIMeanSL(2, snr_index, :, :));

LIM=max(abs(wSMIMeanSL(:)));
LwSMI=LwSMI+LwSMI';
RwSMI=RwSMI+RwSMI';

subplot(2,3,1)
imagesc(LwSMI); set(gca, 'YDir', 'normal');
xlabel('Channels'), ylabel('Channels')
caxis([-1.2 1.2]*LIM); title('LIPL-RMFG Real wSMI')

subplot(2,3,2)
imagesc(RwSMI); set(gca, 'YDir', 'normal');
xlabel('Channels'), ylabel('Channels')
caxis([-1.2 1.2]*LIM); title('RIPL-RMFG Real wSMI')

subplot(2,3,3)
marksize=2;
data=([squeeze((LColl(snr_index, :,:))).'; squeeze((RColl(snr_index,:, :))).'].');
colors=[0.4660, 0.6740, 0.1880; 1.0000,0.5490,0; 0.0039    0.4727    0.4336; 0.5020 0 0.5020];
H=notBoxPlot(data,[]); d=[H.data]; hold on;
set([H.mu],'Color','w');
for ii=1:length(H)
set([H(ii).data],'MarkerSize',marksize,'markerfacecolor',[0.2,1,0.2],'color',[0,0.2,0]);
set([H(ii).semPtch],'FaceColor',colors(ii,:),'EdgeColor','none');
set([H(ii).sdPtch],'FaceColor',colors(ii,:)*0.3,'EdgeColor','none');
end; hold on;
X=NaN(length(H),length(H(1).data.XData));
Y=NaN(length(H),length(H(1).data.YData));
for nd=1:length(H)
    X(nd,:)=[H(nd).data.XData];
    Y(nd,:)=[H(nd).data.YData];
end; clear nd;
plot(X, Y, 'o', 'MarkerSize', marksize, 'markerfacecolor',[1.0,0.2,0.2],'color',[0,0.2,0])
MIN=min([min(LColl(:)), min(RColl(:))]);
MAX=max([max(abs(LColl(:))), max(abs(RColl(:)))]);
ylim([MIN-abs(MIN/10), MAX+abs(MAX/10)])
set(gca,'FontName','Arial','FontSize',8); grid on; box on;
set(gca,'xticklabel', {'LIPL-RMFG', '95\% All', 'RIPL-RMFG', '95\% All'}) %Remove tick labels
xticklabel_rotate([],45,[],'Fontsize',10)
title({'wSMI for closest electrodes' ; 'to sources as compared to others'},'Fontsize',12, 'HorizontalAlignment', 'Center') 
LwSMIF=squeeze(wSMISurrSL(1, snr_index, :, :));
RwSMIF=squeeze(wSMISurrSL(2, snr_index, :, :));
axis square;

LwSMIF=LwSMIF+LwSMIF';
RwSMIF=RwSMIF+RwSMIF';

subplot(2,3,4)
imagesc(LwSMIF); set(gca, 'YDir', 'normal');
xlabel('Channels'), ylabel('Channels')
caxis([-1.2 1.2]*LIM); title('LIPL-RMFG Surrogate wSMI')

subplot(2,3,5)
imagesc(RwSMIF); set(gca, 'YDir', 'normal');
xlabel('Channels'), ylabel('Channels')
caxis([-1.2 1.2]*LIM); title('RIPL-RMFG Surrogate wSMI')

subplot(2,3,6)
marksize=2;
data=([squeeze((LCollF(snr_index, :,:))).'; squeeze((RCollF(snr_index,:, :))).'].');
colors=[0.4660, 0.6740, 0.1880; 1.0000,0.5490,0; 0.0039    0.4727    0.4336; 0.5020 0 0.5020];
H=notBoxPlot(data,[]); d=[H.data]; hold on;
set([H.mu],'Color','w');
for ii=1:length(H)
set([H(ii).data],'MarkerSize',marksize,'markerfacecolor',[0.2,1,0.2],'color',[0,0.2,0]);
set([H(ii).semPtch],'FaceColor',colors(ii,:),'EdgeColor','none');
set([H(ii).sdPtch],'FaceColor',colors(ii,:)*0.3,'EdgeColor','none');
end; hold on;
X=NaN(length(H),length(H(1).data.XData));
Y=NaN(length(H),length(H(1).data.YData));
for nd=1:length(H)
    X(nd,:)=[H(nd).data.XData];
    Y(nd,:)=[H(nd).data.YData];
end; clear nd;
plot(X, Y, 'o', 'MarkerSize', marksize, 'markerfacecolor',[1.0,0.2,0.2],'color',[0,0.2,0])
MIN=min([min(LCollF(:)), min(RCollF(:))]);
MAX=max([max(abs(LCollF(:))), max(abs(RCollF(:)))]);
ylim([MIN-abs(MIN/10), MAX+abs(MAX/10)])
set(gca,'FontName','Arial','FontSize',8); grid on; box on;
set(gca,'xticklabel', {'L-R Real', 'L-R Surr', 'R-R Real', 'R-R Surr'}) %Remove tick labels
xticklabel_rotate([],45,[],'Fontsize',10)
title({'Median wSMI across channels'; 'for real and surrogate data'},'Fontsize',12, 'HorizontalAlignment', 'Center')
colorbar('Position',[0.915 0.11 0.01 0.879715909090909]);
caxis([-1.2 1.2]*LIM);
 
movObj.FrameRate = 1;
        % Write each frame to the file
        writeMovie(movObj,getframe(hf));
        hp4 = get(subplot(2,2,4),'Position')   
        colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(2)+hp4(3)*3]); caxis([0 1])
        
    end
movObj.PlayAllFrames=true;
% Set palindromic looping flag
movObj.Loop = 'backandforth';
% Finish writing movie and close file
close(movObj);
delete(gcp('nocreate'));
end