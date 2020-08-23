%% Determining Topographic Accuracy 100x100 Idea
clear all; close all;
clc

addpath(genpath('/home/laurai/matlab/BBCBCode/tools/'))
addpath('/home/laurai/matlab/BBCB/PlottingScripts/')
addpath('/home/laurai/matlab/fullfig/')
interactionkind={'linear', 'henon', 'ikeda', 'lorenz_xy', 'lorenz_xz', 'lorenz_yz', 'rossler_xy', 'rossler_xz', 'rossler_yz'}; % 'henon', 'ikeda'
datasets=100;
wPLI_score_sum=zeros(2,9,100, 100);
wSMI_score_sum=zeros(2,9,100, 100);
sourceloc={'LIPLMFG', 'RIPLMFG'};

for sl=1:2    
    for idx=1:length(interactionkind)  
        for snr_index=1:100 
            for idata=1:datasets

            load(['/home/laurai/matlab/Simulations/Conn_BBCB_Results_' sourceloc{sl} '_shuffle2/true/dataset_' char(interactionkind(idx))  '_' num2str(snr_index) '_' num2str(idata) '.mat']); 

            wpli=data.wpli;
            wpli=mean(abs(wpli),2); %max(abs(wpli), [], 2);

            wpli=squareform(wpli);
            wsmi=data.wsmi;
            
            % Electrodes closest to imposed sources (Euclidean Distance).
            r_wpli=wpli(51, 87); 
            l_wpli=wpli(51, 86); 

            r_wsmi=wsmi(51, 87);
            l_wsmi=wsmi(51, 86);

            wpli(wpli==0)=[];
            wsmi(wsmi==0)=[];

            end_point_wpli=prctile(wpli, 95);
            end_point_wsmi=prctile(wsmi, 95);

              if  sl==1
                  score_wpli=l_wpli>end_point_wpli;
                  score_wsmi=l_wsmi>end_point_wsmi;
              else
                  score_wpli=r_wpli>end_point_wpli;
                  score_wsmi=r_wsmi>end_point_wsmi;
              end

              wPLI_score_sum(sl,idx,snr_index,idata)=score_wpli;
              wSMI_score_sum(sl,idx,snr_index,idata)=score_wsmi;

            end
        end
    end
end

%% Topographic accuracy
%load/save('/home/laurai/matlab/EEGSimulationStudy/TestPerformance/topographicacc.mat', 'wPLI_score_sum', 'wSMI_score_sum', 'p', 'p_snr')

idx_rev = [1 2 3 7 8 9 4 5 6];
Accuracy_wPLI=zeros(2,9);
Accuracy_wSMI=zeros(2,9);
p=NaN(2,9);
nPerm=10000;
for sl=1:2
    for idx=1:9
        i=idx_rev(idx);
        Accuracy_wPLI(sl,idx)=mean(sum(wPLI_score_sum(sl,i,:, :),4))/100;
        Accuracy_wSMI(sl,idx)=mean(sum(wSMI_score_sum(sl,i,:, :),4))/100;            
        p(sl, idx)=permtestcomplete(squeeze(sum(wPLI_score_sum(sl,i,:, :),4)), squeeze(sum(wSMI_score_sum(sl,i,:, :),4)), 'Unpaired', 'both', nPerm);
    end
end

interactionkind_des={'linear', 'H\''{e}non', 'Ikeda','Lorenz (x,y)', 'Lorenz (x,z)', 'Lorenz (y,z)', 'R{\"o}ssler (x,y)', 'R{\"o}ssler (x,z)', 'R{\"o}ssler (y,z)'} % 'henon', 'ikeda'
interactionkind_des_rev={'Linear', 'H\''{e}non map', 'Ikeda map', 'R{\"o}ssler (x,y)', 'R{\"o}ssler (x,z)', 'R{\"o}ssler (y,z)', 'Lorenz (x,y)', 'Lorenz (x,z)', 'Lorenz (y,z)'} % 'henon', 'ikeda'

%% Plot Topographic Accuracy over all SNRs.
set(0, 'defaultTextInterpreter', 'latex');
H=fullfig;
b_left=([squeeze(Accuracy_wPLI(1,:)); squeeze(Accuracy_wSMI(1,:))]);
b_right=([squeeze(Accuracy_wPLI(2,:)); squeeze(Accuracy_wSMI(2,:))]);
b_both=imagesc([b_left; b_right]);
title('wPLI and wSMI: Mean Topographic Accuracy', 'FontSize',16,'FontWeight', 'bold' )
%axis square
colormap(flipud(hot));
caxis([0,1]);
set(gca,'yticklabel',[], 'xticklabel', []) %Remove tick labels
%% Get tick mark positions
yTicks = get(gca,'ytick');
yTicklabel={' ', 'wPLI (L-R)',' ', 'wSMI (L-R)', ' ', 'wPLI (R-R)' , ' ', 'wSMI (R-R)', ' '}
xTicks = get(gca, 'xtick');
ax = axis; %Get left most x-position
HorizontalOffset = 0.1;
%% Reset the ytick labels in desired font
for i = 1:length(yTicks)
%Create text box and set appropriate properties
     text(ax(1) - HorizontalOffset,yTicks(i),[yTicklabel{i}],...
         'HorizontalAlignment','Right','FontSize',14,'interpreter', 'latex');   
end
%% Reset the xtick labels in desired font 
minY = max(yTicks);
verticalOffset = - 0.2;
for xx = 1:length(xTicks)
%Create text box and set appropriate properties
     text(xTicks(xx), minY - verticalOffset, [interactionkind_des_rev{xx}],...
         'HorizontalAlignment','Center','FontSize',12,'interpreter', 'latex');   
end

for sl=1:2
    for dyn=1:9
        if p(sl, dyn)<=0.05
            if sl==1
                line([dyn, dyn],[1,2], 'Color', ([102 255 0]./256), 'LineWidth', 10.0, 'LineStyle', '-');
            else 
                line([dyn, dyn],[3,4], 'Color', ([102 255 0]./256), 'LineWidth', 10.0, 'LineStyle', '-');
            end
        end
    end
end

colorbar('Ticks',[0.1:0.05:1],...
         'TickLabels',{'0.1','','0.2','','0.3','','0.4','','0.5','','0.6','','0.7','','0.8','','0.9',' ','1'})
     
export_fig(['/home/laurai/matlab/BBCB/Figures/2019_TOPO_FV'],...
    ['-r' num2str(150)], '-a2', '-nocrop', '-transparent');

%% Plot Topographic Accuracy w.r.t. SNR
p_snr=NaN(2,9,100);

for sl=1:2
    for idx=1:9
        for snr=1:100  
            i=idx_rev(idx);  
            p_snr(sl, idx, snr)=permtestcomplete(squeeze(wPLI_score_sum(sl,i,snr, :)), squeeze(wSMI_score_sum(sl,i,snr, :)), 'Unpaired', 'both', nPerm);
        end
    end
end

colors=[0.4660, 0.6740, 0.1880; 0.4940, 0.1840, 0.5560];

set(0, 'defaultTextInterpreter', 'latex');
fullfig;
suptitle({['wPLI and wSMI: Topographic Accuracy in Dependence of SNR']; ' '})
CollSigs=zeros(9,100);
for idx=1:9
    x=0.01:0.01:1;
    subplot(3,3,idx)
    i=idx_rev(idx);
    y1=squeeze(sum(wPLI_score_sum(1,i,:,:),4)/100);
    y2=squeeze(sum(wSMI_score_sum(1,i,:,:),4)/100);
    y3=squeeze(sum(wPLI_score_sum(2,i,:,:),4)/100);
    y4=squeeze(sum(wSMI_score_sum(2,i,:,:),4)/100);
    
    plot(x,y1,':','color',colors(1,:), 'LineWidth', 2)
    hold on
    plot(x,y2,':','color', colors(2,:), 'LineWidth', 2) 
    hold on
    plot(x,y3,'--','color',colors(1,:), 'LineWidth', 2)
    hold on
    plot(x,y4,'--','color', colors(2,:), 'LineWidth', 2) 
    hold on
    title([ char(interactionkind_des_rev(idx)) ' dynamics']) 
    xlabel('SNR')
    ylabel('Accuracy')
    ylim([0,1.3])
    amp_left=1.05;
    amp_right=0.05;
    pval_left=squeeze(p_snr(1,idx,:));
    pval_right=squeeze(p_snr(2,idx,:));
    signific_left = NaN(1, 100); signific_left(pval_left<0.05) = 1;
    signific_right = nan(1, 100); signific_right(pval_right<0.05) = 1;
    signific_both = (signific_left  + signific_right) - 1;
    signific_both(signific_both==0)=NaN;
    plot(x, signific_both *amp_left, '.k');
    text(0.6, amp_left+0.1, '$p_{LR} \, \& \, p_{RR} <0.05$'); 
    CollSigs(idx, :)=signific_both;

end
hp4 = get(subplot(3,3,9),'Position')   
h=legend('wPLI (L-R)', 'wSMI (L-R)', 'wPLI (R-R)','wSMI (R-R)', 'Location', 'SouthEast')
set(h, 'Position',[0.338906648507749 0.0209167143129407 0.361177728435142 0.0359133795903264],...
    'Orientation','horizontal');
export_fig(['/home/laurai/matlab/BBCB/Figures/2019_TOPO_FV_PvsSNR'],...
     ['-r' num2str(150)], '-a2', '-nocrop', '-transparent');

% load/save('/home/laurai/matlab/BBCB/PlottingScripts/topographicacc.mat', 'wPLI_score_sum', 'wSMI_score_sum', 'p', 'p_snr')