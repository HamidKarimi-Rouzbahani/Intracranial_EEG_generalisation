% This code plots the AUCs obtained from the classifications
% of contacts in each of the 3 time sub-windows in
% Within-patient classification (interictal and ictal)

% INPUTS: data from C5_Permuting_3_classifciation_windows
% OUTPUTS: figures and numbers produced in command window
%%
clc;
clear all;
close all;
addpath(genpath('bayesFactor-master')) % Add the Bayes Factor analysis toolbox to Matlab path
load('random_permutations_within_time_wind.mat') %  Load the data and the permutation data
datas=[];
datas_random=[];
metric=8; %8=AUC
for tw=1:3 % time windows
    data_inter=squeeze(nanmean(Performance_within_inter(tw,:,:,metric),2));
    data_ictal=squeeze(nanmean(Performance_within_ictal(tw,:,:,metric),2));
    datas_random_inter=squeeze(nanmean(Performance_rand_within_inter(tw,:,:,metric,:),2));
    datas_random_ictal=squeeze(nanmean(Performance_rand_within_ictal(tw,:,:,metric,:),2));
    
    datas=horzcat(datas,[data_inter data_ictal]);
    datas_random=horzcat(datas_random,[datas_random_inter datas_random_ictal]);
end
colours={[0 0 0],[0.3 0.3 0.3],[0 0 0],[0.3 0.3 0.3],[0 0 0],[0.3 0.3 0.3]};
datas_tmp_outed=nan(55,6);
for i=1:6
    tmpp=rmoutliers(datas(:,i));
    datas_tmp_outed(1:length(tmpp),i)=tmpp;
    swarmchart([i*ones(length(tmpp),1)],tmpp,'MarkerFaceColor',colours{i},'MarkerEdgeColor',colours{i},'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
    hold on;
end
boxplot(datas_tmp_outed,'Whisker',inf,'color','k');
xticks([1:6])
ylabel('Area Under Curve (AUC)')
box off
grid on
set(gca,'TickDir','out','Fontsize',16)
title('Within-Patient Classification')
xticklabels({'Interictal (early)';'Ictal (early)';...
    'Interictal (middle)';'Ictal (middle)';...
    'Interictal (late)';'Ictal (late)'})
ylim([0.35 1])
plot([0.5 6.5],[0.5 0.5],'--k')
xlim([0.5 6.5])
for i=1:6
    bays_against_chance=bf.ttest2(datas(:,i),datas_random(:,i))
    [mean(datas(:,i)) std(datas(:,i))]
end
bays_diff=[bf.ttest(datas(:,1),datas(:,2)) bf.ttest(datas(:,3),datas(:,4)) bf.ttest(datas(:,5),datas(:,6))]
