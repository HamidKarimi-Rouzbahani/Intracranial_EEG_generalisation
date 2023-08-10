%% Plotting each feature
clc;
clear all;
close all;
inter_to_ictal=2;
if inter_to_ictal==1
    load(['Generalisation_performance_across_time_all_feats_comb_crcted_feats_ovrsmp_imp_100_10.mat'])
else
    load(['Generalisation_performance_across_time_ict_to_int_all_feats_comb_crcted_feats_ovrsmp_imp_100_10.mat'])
end
Measures={'accuracy','sensitivity','specificity','precision','recall','f1-measure','gmean','auc'};

classif=1; %
ff=1;
iter_rand=1;
for iter_equalis=1:size(Ground_truth,2)
    for p=1:size(Ground_truth,3)
        Performance(ff,iter_equalis,p,iter_rand,1:7)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif});
        [~,~,~,Performance(ff,iter_equalis,p,iter_rand,8)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif},1);
        for rep=1:1
            Predictions_tmp=randsample(Predictions{ff,iter_equalis,p,iter_rand,classif},length(Predictions{ff,iter_equalis,p,iter_rand,classif}));
            Performance_rand(ff,iter_equalis,p,iter_rand,1:7,rep)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions_tmp);
            [~,~,~,Performance_rand(ff,iter_equalis,p,iter_rand,8,rep)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions_tmp,1);
        end
    end
end
Performance=squeeze(Performance);
Performance_rand=squeeze(Performance_rand);

resect=bar([1:8],[nanmean(Performance)],0.6,'k');
hold on;
errorbar([1:8],[nanmean(Performance)],[nanstd(Performance)]/sqrt(size(Performance,1)),"LineStyle","none","Color","k")
grid on;
title(['Cross-time generalisation'])
xticks([1:8])
xticklabels(Measures)
ylim([0 1])

sig_level=6;
for metric=1:8
    bays=bf.ttest2(Performance(:,metric),squeeze(reshape(Performance_rand(:,metric,:),[size(Performance_rand(:,metric,:),1)*size(Performance_rand(:,metric,:),3)],[])))
    if bays>sig_level
        plot(metric,nanmean(Performance(:,metric))+0.1,'*r')
    end
end
%% Features importance
feats=[1:8 10:12 14:36];
Features_labels={'Mean','Median','Variance','Skewness','Kurtosis',...
    'LZ Comp','Higichi FD','Katz FD','Lyap Exp','Hurst Exp',...
    'Samp Ent','Apprx Ent','Autocor','Hjorth Comp','Hjorth Mob',...
    'Mean Freq','Med Freq','Avg Freq','SEF','Pow Med Freq',...
    'Phs Med Freq','Power','Energy ratio','Delta','Theta',...
    'Alpha','Beta','Gamma','H-Gamma','Correl','Delta Coh',...
    'Theta Coh','Alpha Coh','Beta Coh','Gamma Coh','HGamma Coh'};
iter_rand=1;
feat_imp=nan(1,34,1);
for iter_equalis=1
    for p=1:55
        feat_imp(p,included_feats{ff,iter_equalis,p,iter_rand},iter_equalis)=impCART{ff,iter_equalis,p,iter_rand};
    end
end
feat_imp(feat_imp==0)=nan;
figure;
bar([1:length(feats)],[nanmean(feat_imp,1)],0.6,'k');
hold on;
errorbar([1:length(feats)],[nanmean(feat_imp,1)],[nanstd(feat_imp)]/sqrt(size(feat_imp,1)),"LineStyle","none","Color","k")
grid on;
ylabel('Feature importance')
xticks([1:length(feats)])
xticklabels([Features_labels(feats)])
%% BF Anova
clc
clear all
participants_info=tdfread(['F:\RESEARCH\Hamid\Multicentre dataset\ds004100\participants.tsv']);
subjects_analysed=[1:25 27:39 41:45 47:58]; % already remove 40 becuase "data" does not have that
inter_to_ictal=1;
if inter_to_ictal==1
    load(['Generalisation_performance_across_time_all_feats_comb_crcted_feats_ovrsmp_imp_100_10.mat'])
else
    load(['Generalisation_performance_across_time_ict_to_int_all_feats_comb_crcted_feats_ovrsmp_imp_100_10.mat'])
end
classif=1; %
ff=1;
iter_rand=1;
for iter_equalis=1:size(Ground_truth,2)
    for p=1:size(Ground_truth,3)
        Performance(ff,iter_equalis,p,iter_rand,1:7)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif});
        [~,~,~,Performance(ff,iter_equalis,p,iter_rand,8)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif},1);
    end
end

feats=1;
metric=8;
Measures={'accuracy','sensitivity','specificity','precision','recall','f1-measure','gmean','auc'};

data=squeeze(nanmean(Performance(feats,:,:,1,metric),2));
pt=0;
for P=subjects_analysed
    pt=pt+1;
    g1{pt,1}=participants_info.engel(P,1);
    g2{pt,1}=participants_info.target(P,:);
    g3{pt,1}=participants_info.lesion_status(P,:);
    g4{pt,1}=participants_info.implant(P,:);
    participants_info.participant_id(P,:)
end
clearvars data_ready
nans=isnan(data);
nans([1])=1; % bad labeling of lesion (1)
data_ready=table(data,g1,g2,g3,g4);
data_ready(nans,:)=[];
data_ready.Properties.VariableNames = {'Effect' 'Outcome' 'ROI' 'Lesion' 'Recording'};
bfFull= bf.anova(data_ready,'Effect ~ Outcome+ROI+Lesion+Recording');
bfRestricted= bf.anova(data_ready,'Effect ~ ROI+Lesion+Recording');
effects(1) = bfFull/bfRestricted;
bfRestricted= bf.anova(data_ready,'Effect ~ Outcome+Lesion+Recording');
effects(2) = bfFull/bfRestricted;
bfRestricted= bf.anova(data_ready,'Effect ~ Outcome+ROI+Recording');
effects(3) = bfFull/bfRestricted;
bfRestricted= bf.anova(data_ready,'Effect ~ Lesion+Outcome+ROI');
effects(4) = bfFull/bfRestricted;

% bar(effects,0.2,'k');
% xticks([1:length(effects)])
% xticklabels({'Outcome','ROI','Lesion','Recording'})
% xtickangle(45)
% ylabel('Bayes factor')
% grid on;



% Plotting accuracies across variables
figure
clearvars data
for condition=1:4
    subplot(2,2,condition)
    % condition=4; %1= Engel levels; 2= roi; 3= Lesional?; 4= Recording
    feats=1;
    data(feats,:)=squeeze(nanmean(Performance(feats,:,:,1,metric),2));
    cond=zeros(size(data,2),1);
    
    if condition==1
        strings={'1','2','3','4'};
        col=cool;
        shiftee=0.1;
        s=0;
        for i=subjects_analysed
            s=s+1;
            for c=1:length(strings)
                if contains(participants_info.engel(s,1),strings{c})
                    cond(s)=(c);
                end
            end
        end
    elseif condition==2
        strings={'FRONTAL','TEMPORAL','FP','MTL','PARIETAL','MFL','INSULAR'};
        col=jet;
        shiftee=0.2;
        s=0;
        for i=subjects_analysed
            s=s+1;
            for c=1:length(strings)
                if contains(participants_info.target(s,:),strings{c})
                    cond(s)=(c);
                end
            end
        end
    elseif condition==3
        strings={'LESIONAL','NON-LESIONAL'};
        col=cool;
        shiftee=0.1;
        s=0;
        for i=subjects_analysed
            s=s+1;
            for c=1:length(strings)
                if contains(participants_info.lesion_status(s,:),strings{c})
                    cond(s)=(c);
                end
            end
        end
    elseif condition==4
        strings={'SEEG','ECOG'};
        col=cool;
        shiftee=0.1;
        s=0;
        for i=subjects_analysed
            s=s+1;
            for c=1:length(strings)
                if contains(participants_info.implant(s,:),strings{c})
                    cond(s)=(c);
                end
            end
        end
    end
    feats=1;
    bar_width=0.2;
    tmp=round(linspace(1,length(col),length(strings)));
    cols={col(tmp,:)};
    multip=floor(length(unique(cond))/2)*0.2;
    uniques=unique(cond);
    for i=uniques'
        cond(40)=0; % unused
        try
            plots{i}=bar(bar_width*(i-1),nanmean(data(feats,cond==i)),bar_width,'facecolor',cols{1,1}(i,:));
            hold on;
            errorbar(bar_width*(i-1),nanmean(data(feats,cond==i)),nanstd(data(feats,cond==i))/sqrt(size(data(feats,cond==i),2)),"LineStyle","none","Color","k")
        end
    end
    
    xes(feats)=bar_width*(ceil((length(uniques)-1)/2)-1)+feats*(1+multip)+shiftee;
    plot([-0.5 bar_width*(i-1)+0.5],[0.5 0.5],'--k');
    xticks([xes])
    xticklabels([])
    xlim([-0.5 bar_width*(i-1)+0.5])
    grid on;
    ylim([0 1])
    ylabel('Generalisation accuracy')
    
    if condition==1
        legend([plots{1} plots{2} plots{3} plots{4}],{['LEVEL ',strings{1}],['LEVEL ',strings{2}],['LEVEL ',strings{3}],['LEVEL ',strings{4}]},'location','northwest')
    elseif condition==2
        legend([plots{1} plots{2} plots{3} plots{4} plots{5} plots{6} plots{7}],{strings{1},strings{2},strings{3},strings{4},strings{5},strings{6},strings{7}},'location','northwest')
    elseif condition==3
        legend([plots{1} plots{2}],{strings{1},strings{2}},'location','northwest')
    elseif condition==4
        legend([plots{1} plots{2}],{strings{1},strings{2}},'location','northwest')
    end
    
    title(['Interictal to ictal generalisation; BF= ',num2str(round(effects(condition)))]);
    
    %% BFs across conditions    
    pt=0;
    for P=subjects_analysed
        pt=pt+1;
        if condition==1
            g{pt,1}=participants_info.engel(P,1);
        elseif condition==2
            g{pt,1}=participants_info.target(P,:);
        elseif condition==3
            g{pt,1}=participants_info.lesion_status(P,:);
        elseif condition==4
            g{pt,1}=participants_info.implant(P,:);
        end
    end
    
    data_bf=data;
    nans=isnan(data_bf);
    nans([1])=1; % bad labeling of lesion (1)
    data_bf(nans)=[];
    g(nans)=[];
    uniq_conds_tmp=unique(g);

    combs_tmp=nchoosek(1:length(uniq_conds_tmp),2);
    for comb=1:size(combs_tmp,1)
        conditions_BFtmp(comb,1)=bf.ttest2(data_bf(strcmp(g,uniq_conds_tmp(combs_tmp(comb,1)))),data_bf(strcmp(g,uniq_conds_tmp(combs_tmp(comb,2)))));        
        combinations_tmp{comb,1}=uniq_conds_tmp(combs_tmp(comb,1));
        combinations_tmp{comb,2}=uniq_conds_tmp(combs_tmp(comb,2));
    end
    conditions_BF{condition}=conditions_BFtmp;
    combinations{condition}=combinations_tmp;
    clearvars data_bf conditions_BFtmp combinations_tmp
end