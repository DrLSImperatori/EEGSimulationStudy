%% Whole-Brain Accuracy Computation and Plot
close all;
addpath('/home/laurai/matlab/palm-alpha114/') % Palm Pareto Toolbox
addpath('/home/laurai/matlab/EEGSimulationStudy/TestPerformance/')
addpath('/home/laurai/matlab/fullfig/')

interactionkind={'linear', 'henon', 'ikeda', 'lorenz_xy', 'lorenz_xz', 'lorenz_yz', 'rossler_xy', 'rossler_xz', 'rossler_yz'};
SurrogateDataType = {'shuffle', 'AAFT'};
datasets=100;
sourceloc={'LIPLMFG', 'RIPLMFG'};

for surr = 1:2
    SurrType=SurrogateDataType{surr};
    for sl=1:2
        for idx=1:9
            P_wPLI=NaN(100,100);
            P_wSMI=NaN(100,100);
            parfor snr_index=1:100   
                wPLI_surrogate=zeros(100,1);
                wSMI_surrogate=zeros(100,1); 
                for idata=1:datasets
                    Struct1=load(['/home/laurai/matlab/Simulations/Conn_BBCB_Results_' sourceloc{sl} '_' SurrType '2/false/dataset_' char(interactionkind(idx))  '_' num2str(snr_index) '_' num2str(idata) '.mat']);
                    wplif=Struct1.data.wpli_f;
                    wplif=mean(abs(wplif),2);

                    wsmif=Struct1.data.wsmi_f;

                    wplif(wplif==0)=[];
                    wsmif(wsmif==0)=[];

                    wPLI_surrogate(idata)=median(wplif);
                    wSMI_surrogate(idata)=median(wsmif); 
                end
                P1=NaN(1,100);
                P2=NaN(1,100);
                for idata=1:datasets
                    Struct2=load(['/home/laurai/matlab/Simulations/Conn_BBCB_Results_' sourceloc{sl} '_' SurrType '2/true/dataset_' char(interactionkind(idx))  '_' num2str(snr_index) '_' num2str(idata) '.mat']);
                    wpli=Struct2.data.wpli;
                    wpli=mean(abs(wpli),2);

                    wsmi=Struct2.data.wsmi; 

                    wpli(wpli==0)=[];
                    wsmi(wsmi==0)=[];

                    P1(idata)= palm_pareto(median(wpli), wPLI_surrogate, false, 0.05, false); %sum(wPLI_surrogate(:)>=wpli_real)/datasets;
                    P2(idata)= palm_pareto(median(wsmi), wSMI_surrogate, false, 0.05, false); %sum(wSMI_surrogate(:)>=wsmi_real)/datasets; 
                end
                P_wPLI(snr_index, :)=P1;
                P_wSMI(snr_index, :)=P2;
            end
        Res=struct();
        Res.P_wPLI=P_wPLI;
        Res.P_wSMI=P_wSMI;   
        savetofile(['/home/laurai/matlab/BBCB/WB_PColl_' SurrType '/P_' num2str(sl) '_' num2str(idx) '.mat'], Res);
        delete(gcp('nocreate'));
        end
    end

P=zeros(2,9,2, 100,100);

for sl=1:2
    for idx=1:9
            load(['/home/laurai/matlab/BBCB/WB_PColl_' SurrType '/P_' num2str(sl) '_' num2str(idx) '.mat']);
            P(sl,idx,1, :,:)=data.P_wPLI;
            P(sl,idx,2, :,:)=data.P_wSMI;
    end
end

P_b=P;
P_b(P_b>=0.05)=NaN;
P_b(P_b<0.05)=1;
P_b(isnan(P_b))=0;

nPerm=10000;

%% Mean Accuracy

idx_rev = [1 2 3 7 8 9 4 5 6];
Accuracy_wPLI=zeros(2,9);
Accuracy_wSMI=zeros(2,9);

p=NaN(2,9);

for sl=1:2
    for idx=1:9
        i=idx_rev(idx);
        Accuracy_wPLI(sl,idx)=mean(squeeze(sum(P_b(sl,i,1,:,:), 5)/100)); %sum(P(sl,i,1,:)<=0.05)/100;
        Accuracy_wSMI(sl,idx)=mean(squeeze(sum(P_b(sl,i,2,:,:), 5)/100)); %sum(P(sl,i,2,:)<=0.05)/100;
        p(sl,idx)=permtestcomplete(squeeze(sum(P_b(sl,i,1,:,:), 5)),squeeze(sum(P_b(sl,i,2,:,:), 5)), 'Unpaired', 'both', nPerm);
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
title('wPLI and wSMI: Mean Whole-Brain Accuracy', 'FontSize',16,'FontWeight', 'bold' )
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
export_fig(['/home/laurai/matlab/BBCB/Figures/2019_WB_FV_' SurrType],...
     ['-r' num2str(150)], '-a2', '-nocrop', '-transparent');

%% Plot accuracy w.r.t. to SNR
colors=[0.4660, 0.6740, 0.1880; 0.4940, 0.1840, 0.5560];
interactionkind_des={'linear', 'H\''{e}non', 'Ikeda','Lorenz (x,y)', 'Lorenz (x,z)', 'Lorenz (y,z)', 'R{\"o}ssler (x,y)', 'R{\"o}ssler (x,z)', 'R{\"o}ssler (y,z)'} % 'henon', 'ikeda'
interactionkind_des_rev={'Linear', 'H\''{e}non map', 'Ikeda map', 'R{\"o}ssler (x,y)', 'R{\"o}ssler (x,z)', 'R{\"o}ssler (y,z)', 'Lorenz (x,y)', 'Lorenz (x,z)', 'Lorenz (y,z)'} % 'henon', 'ikeda'
ind=[1,2,3,7,8,9,4,5,6];
p_snr=zeros(2,9,100);
for sl=1:2
    for idx=1:9
        for snr=1:100  
        i=ind(idx);  
        p_snr(sl,idx, snr)=permtestcomplete(squeeze(P_b(sl,i,1,snr, :)), squeeze(P_b(sl,i,2,snr, :)), 'Unpaired', 'both', nPerm);
        end
    end
end


set(0, 'defaultTextInterpreter', 'latex');
fullfig;
suptitle({['wPLI and wSMI: Whole-Brain Accuracy in Dependence of SNR']; ' '})
CollSigs=zeros(9,100);
for idx=1:9
    x=0.01:0.01:1;
    subplot(3,3,idx)
    
    i=ind(idx);
    y1=squeeze((sum(P_b(1,i,1,:,:), 5)/100)); %smooth( , 20)
    y2=squeeze((sum(P_b(1,i,2,:,:), 5)/100)); %smooth(squeeze(1-P(1,i,2,:)), 20);   
    y3=squeeze((sum(P_b(2,i,1,:,:), 5)/100)); %smooth(squeeze(1-P(2,i,1,:)), 20);
    y4=squeeze((sum(P_b(2,i,2,:,:), 5)/100)); %smooth(squeeze(1-P(2,i,2,:)), 20);   
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
    CollSigs(idx, :)=signific_both;
    plot(x, signific_both *amp_left, '.k');
    text(0.6, amp_left+0.1, '$p_{LR} \, \& \, p_{RR} < 0.05$');    
end
hp4 = get(subplot(3,3,9),'Position')   
h=legend('wPLI (L-R)', 'wSMI (L-R)', 'wPLI (R-R)','wSMI (R-R)', 'Location', 'southeast')
set(h, 'Position',[0.338906648507749 0.0209167143129407 0.361177728435142 0.0359133795903264],...
    'Orientation','horizontal');
export_fig(['/home/laurai/matlab/BBCB/Figures/2019_WB_FV_PvsSNR_' SurrType],...
     ['-r' num2str(150)], '-a2', '-nocrop', '-transparent');

save(['/home/laurai/matlab/EEGSimulationStudy/TestPerformance/wholebrain2019_' SurrType '.mat'], 'P', 'p', 'p_snr')

end