% This code plots the AUCs obtained from the classifications
% of contacts in 6 classifications including:
% Within-patient classification (interictal and ictal)
% Across-patient generalisation (interictal and ictal)
% Across-time generalisation (interictal to ictal and vice versa)

% INPUTS: data from C5_Permuting_all_classifciation_windows
% OUTPUTS: figures and numbers produced in command window
%%
clc;
clear all;
close all;
addpath(genpath('bayesFactor-master'))
load('random_permutations.mat')
what=2;% 1=within patient; 2= across time; 3=across patients
Measures={'accuracy','sensitivity','specificity','precision','recall','f1-measure','gmean','auc'};
metric=8; %AUC
if what==1
    data_inter=Performance_within_inter(:,metric);
    data_ictal=Performance_within_ictal(:,metric);
    datas_random_inter=squeeze(nanmean(Performance_rand_within_inter(:,metric,:),3));
    datas_random_ictal=squeeze(nanmean(Performance_rand_within_ictal(:,metric,:),3));
elseif what==2
    data_inter=Performance_acrss_time_inter(:,metric);
    data_ictal=Performance_acrss_time_ictal(:,metric);
    datas_random_inter=squeeze(nanmean(Performance_rand_acrss_time_inter(:,metric,:),3));
    datas_random_ictal=squeeze(nanmean(Performance_rand_acrss_time_ictal(:,metric,:),3));
elseif what==3
    data_inter=Performance_acrss_subj_inter(:,metric);
    data_ictal=Performance_acrss_subj_ictal(:,metric);
    datas_random_inter=squeeze(nanmean(Performance_rand_acrss_subj_inter(:,metric,:),3));
    datas_random_ictal=squeeze(nanmean(Performance_rand_acrss_subj_ictal(:,metric,:),3));
end

datas=[data_inter data_ictal];
datas_random=[datas_random_inter datas_random_ictal];
colours={[0 0 0],[0.3 0.3 0.3]};
datas_tmp_outed=nan(55,2);
for i=1:2
    tmpp=rmoutliers(datas(:,i));
    datas_tmp_outed(1:length(tmpp),i)=tmpp;
    swarmchart([i*ones(length(tmpp),1)],tmpp,'MarkerFaceColor',colours{i},'MarkerEdgeColor',colours{i},'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
    hold on;
end
boxplot(datas_tmp_outed,'Whisker',inf,'color','k');
xticks([1 2])
ylabel('Area Under Curve (AUC)')
box off
grid on
set(gca,'TickDir','out','Fontsize',16)
if what==1
    title('Within-Patient Classification')
    xticklabels({'Interictal';'Ictal'})
    ylim([0.35 1])
    %     ylim([0.8 1])
    %     yticks([0.8 0.9 1])
    %     yticklabels({'0.8','0.95','1'})
elseif what==2
    title('Cross-Time Generalisation')
    xticklabels({'Interictal-to-Ictal';'Ictal-to-Interictal'})
    ylim([0.35 1])
    yticks([0.4:0.1:1])
elseif what==3
    title('Cross-Patient Generalisation')
    xticklabels({'Interictal';'Ictal'})
    ylim([0.35 1])
    yticks([0.4:0.1:1])
end

plot([0.5 2.5],[0.5 0.5],'--k')
bays_against_chance=[bf.ttest2(datas(:,1),datas_random(:,1)) bf.ttest2(datas(:,2),datas_random(:,2))]
bays_diff=bf.ttest(datas(:,1),datas(:,2))
[mean(datas(:,1)) std(datas(:,1)) mean(datas(:,2)) std(datas(:,2))]
%% Features importance
clc
close all
feats=[1:8 10:12 14:36];
Features_labels={'Mean','Median','Variance','Skewness','Kurtosis',...
    'LZ Comp','Higichi FD','Katz FD','Lyap Exp','Hurst Exp',...
    'Samp Ent','Apprx Ent','Autocorr','Hjorth Comp','Hjorth Mob',...
    'Mean Freq','Med Freq','Avg Freq','SEF','Pow Med Freq',...
    'Phs Med Freq','Power','Energy Ratio','Delta Pow','Theta Pow',...
    'Alpha Pow','Beta Pow','Gamma Pow','H-Gamma Pow','Correlation','Delta Coh',...
    'Theta Coh','Alpha Coh','Beta Coh','Gamma Coh','H-Gamma Coh'};

what=2;% 1=within patient; 2= across time; 3=across patients
ictal_or_inter='Interictal';
sorted=1;

if what==1
    load(['Within_subject_performance_sepnorm',ictal_or_inter,'_all_feats_crcted_feats_imp_100_10.mat'])
elseif what==2
    if strcmp(ictal_or_inter,'Interictal')
        load(['Generalisation_performance_across_time_all_feats_comb_crcted_feats_ovrsmp_imp_100_10.mat'])
    else
        load(['Generalisation_performance_across_time_ict_to_int_all_feats_comb_crcted_feats_ovrsmp_imp_100_10.mat'])
    end
elseif what==3
    load(['Generalisation_performance_across_subjects_',ictal_or_inter,'_all_feats_comb_crcted_feats_ovrsmp_imp_100_10.mat'])
end

ff=1;
iter_rand=1;
feat_imp=nan(1,34,1);
for fld=1:size(impCART,5)
    for iter_equalis=1
        for p=1:55
            feat_imp(p,included_feats{ff,iter_equalis,p,iter_rand},iter_equalis,fld)=impCART{ff,iter_equalis,p,iter_rand,fld};
        end
    end
end
feat_imp(feat_imp==0)=nan;
figure;
datas=nanmean(feat_imp,4);
datas_tmp_outed=nan(55,34);
orders=1:size(datas,2);
if sorted==1
    for i=1:size(datas,2)
        tmpp=rmoutliers(datas(:,i));
        datas_tmp_outed(1:length(tmpp),i)=tmpp;
    end
    [~,order]=sort(nanmedian(datas_tmp_outed),'descend');
    %     [~,order]=sort(nanmean(datas_tmp_outed),'descend');
    orders=order;
end
c=0;
datas_tmp_outed=nan(55,34);
for i=orders
    c=c+1;
    if i<6
        colours=[0.7 0.7 0.5];
    elseif i>5 && i<14
        colours=[0.9 0.5 0.9];
    elseif i>13 && i<28
        colours=[0.3 0.7 0.7];
    else
        colours=[0.5 0.5 0.9];
    end
    tmpp=rmoutliers(datas(:,i));
    datas_tmp_outed(1:length(tmpp),c)=tmpp;
    swarmchart([c*ones(size(tmpp,1),1)]',tmpp,20,'MarkerFaceColor',colours,'MarkerEdgeColor',colours,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.3);
    hold on;
    clearvars tmpp
end
boxplot(datas_tmp_outed,'Whisker',inf,'color','k');
box off
grid on
set(gca,'TickDir','out','Fontsize',16)
ylabel('Feature Contribution (A.U.)')
xticks([1:length(feats)])
if sorted==1
    xticklabels([Features_labels(feats(order))])
else
    xticklabels([Features_labels(feats)])
end
if what==1
    ylim([-0.5 2.5])
    title(['Within-Patient Classification ',ictal_or_inter])
elseif what==2
    ylim([-0.5 4])
    if strcmp(ictal_or_inter,'Interictal')
        title(['Cross-Time Generalisation Interictal to Ictal',])
    else
        title(['Cross-Time Generalisation Ictal to Interictal',])
    end
elseif what==3
    ylim([1 8])
    title(['Cross-Patient Generalisation ',ictal_or_inter])
end
%% BF Anova 4 factors
clc;
clear all;
close all;
addpath(genpath('bayesFactor-master'))
load('random_permutations.mat')
participants_info=tdfread(['F:\RESEARCH\Hamid\Multicentre dataset\ds004100\participants.tsv']);
subjects_analysed=[1:25 27:39 41:45 47:58]; % already remove 40 becuase "data" does not have that
what=1;% 1=within patient; 2= across time; 3=across patients
ictal_or_inter='Interictal'; % 1=ictal; non-1=inter
Measures={'accuracy','sensitivity','specificity','precision','recall','f1-measure','gmean','auc'};
metric=8; %AUC
if what==1
    if strcmp(ictal_or_inter,'Interictal')
        data_tmp=Performance_within_inter(:,metric);
        data_tmp_random=squeeze(nanmean(Performance_rand_within_inter(:,metric,:),3));
    else
        data_tmp=Performance_within_ictal(:,metric);
        data_tmp_random=squeeze(nanmean(Performance_rand_within_ictal(:,metric,:),3));
    end
elseif what==2
    if strcmp(ictal_or_inter,'Interictal')
        data_tmp=Performance_acrss_time_inter(:,metric);
        data_tmp_random=squeeze(nanmean(Performance_rand_acrss_time_inter(:,metric,:),3));
    else
        data_tmp=Performance_acrss_time_ictal(:,metric);
        data_tmp_random=squeeze(nanmean(Performance_rand_acrss_time_ictal(:,metric,:),3));
    end
elseif what==3
    if strcmp(ictal_or_inter,'Interictal')
        data_tmp=Performance_acrss_subj_inter(:,metric);
        data_tmp_random=squeeze(nanmean(Performance_rand_acrss_subj_inter(:,metric,:),3));
    else
        data_tmp=Performance_acrss_subj_ictal(:,metric);
        data_tmp_random=squeeze(nanmean(Performance_rand_acrss_subj_ictal(:,metric,:),3));
    end
end
pt=0;
for P=subjects_analysed
    pt=pt+1;
    g1{pt,1}=participants_info.engel(P,1);
    g2{pt,1}=participants_info.target(P,:);
    g3{pt,1}=participants_info.lesion_status(P,:);
    g4{pt,1}=participants_info.implant(P,:);
end
clearvars data_ready
nans=isnan(data_tmp);
nans([1])=1; % bad labeling of lesion (1)
nans([4 25 30 37 38 39])=1; % resection targets with fewer than 5 patients
data_ready=table(data_tmp,g1,g2,g3,g4);
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

% Plotting accuracies across variables
for condition=1:4
    clearvars data data_random
    figure
    feats=1;
    data(feats,:)=data_tmp;
    data_random(feats,:)=data_tmp_random;
    cond=nan(size(data,2),1);
    
    if condition==1
        % %         strings={'1','2','3','4'};
        ss={'S','F'};
        strings={'Engel I','Engel II-IV'};
        col=cool;
        shiftee=0.1;
        s=0;
        for i=subjects_analysed
            s=s+1;
            for c=1:length(ss)
                %                 if contains(participants_info.engel(s,1),strings{c})
                if contains(participants_info.outcome(s,1),ss{c})
                    cond(s)=(c);
                end
            end
        end
    elseif condition==2
        %         strings={'FRONTAL','TEMPORAL','FP','MTL','PARIETAL','MFL','INSULAR'};
        strings={'FRONTAL','TEMPORAL','MTL'};
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
        cond(isnan(cond))=[];
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
        cond(1)=[]; % unused
        data(1)=[];
        data_random(1)=[];
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
    datas=nan(55,length(uniques));
    c=0;
    datas=nan(55,length(uniques));
    for i=uniques'
        c=c+1;
        cond(40)=0; % unused
        data_tmp_tmp=rmoutliers(data(feats,cond==i));
        data_tmp_tmp_random=data_random(feats,cond==i);
        datas(1:length(data_tmp_tmp),c)=data_tmp_tmp;
        if strcmp(ictal_or_inter,'Interictal')
            colours=[0 0 0];
        else
            colours=[0.3 0.3 0.3];
        end
        swarmchart([c*ones(length(data_tmp_tmp),1)]',data_tmp_tmp,'MarkerFaceColor',colours,'MarkerEdgeColor',colours,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.3);
        hold on;
        bays_against_chance{condition,c}=bf.ttest2(data_tmp_tmp,data_tmp_tmp_random);
    end
    [nanmean(datas) nanstd(datas)]
    boxplot(datas,'Whisker',inf,'color','k');
    
    box off
    grid on
    set(gca,'TickDir','out','Fontsize',16)
    xticks([1:size(datas,2)])
    xlim([0.5 size(datas,2)+0.5])
    plot([0.5 size(datas,2)+0.5],[0.5 0.5],'--k')
    title(['ANOVA BF = ',num2str(effects(condition))])
    if what==1
        ylabel('Within-Patient Classification (AUC)')
        %         ylim([0.8 1.0])
        %         yticks([0.8 0.9 1])
        ylim([0.35 1.0])
        yticks([0.4:0.1:1])
    elseif what==2
        ylabel('Cross-Time Generalisation (AUC)')
        ylim([0.35 1])
        yticks([0.4:0.1:1])
    elseif what==3
        ylabel('Cross-Patient Generalisation (AUC)')
        ylim([0.35 1])
        yticks([0.4:0.1:1])
    end
    if condition==2
        %        strings={'FRT','TPR','FPL','MTL','PRT','MFL','INS'};
        strings={'FRT','TPR','MTL'};
    end
    xticklabels(strings)
    if condition==1
        %         xlabel('Outcome (Engel Level)')
        xlabel('Outcome')
    elseif condition==2
        xlabel('Region of Resection')
    elseif condition==3
        xlabel('Pathology')
    elseif condition==4
        xlabel('Recording')
    end
    
    
    clearvars datas
    % BFs across conditions
    pt=0;
    for P=subjects_analysed
        pt=pt+1;
        if condition==1
            %             g{pt,1}=participants_info.engel(P,1);
            g{pt,1}=participants_info.outcome(P,1);
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
    nans([4 25 30 37 38 39])=1; % resection targets with fewer than 5 patients
    data_bf(nans)=[];
    g(nans)=[];
    uniq_conds_tmp=unique(g);
    
    combs_tmp=nchoosek(1:length(uniq_conds_tmp),2);
    for comb=1:size(combs_tmp,1)
        try
            conditions_BFtmp(comb,1)=bf.ttest2(data_bf(strcmp(g,uniq_conds_tmp(combs_tmp(comb,1)))),data_bf(strcmp(g,uniq_conds_tmp(combs_tmp(comb,2)))));
        catch
            conditions_BFtmp(comb,1)=nan;
        end
        combinations_tmp{comb,1}=uniq_conds_tmp(combs_tmp(comb,1));
        combinations_tmp{comb,2}=uniq_conds_tmp(combs_tmp(comb,2));
    end
    combinations_BF{condition}=conditions_BFtmp;
    combinations{condition}=combinations_tmp;
    clearvars data_bf conditions_BFtmp combinations_tmp
end

%% Features importance-Outcome anova
clc
close all
feats=[1:8 10:12 14:36];
Features_labels={'Mean','Median','Variance','Skewness','Kurtosis',...
    'LZ Comp','Higichi FD','Katz FD','Lyap Exp','Hurst Exp',...
    'Samp Ent','Apprx Ent','Autocorr','Hjorth Comp','Hjorth Mob',...
    'Mean Freq','Med Freq','Avg Freq','SEF','Pow Med Freq',...
    'Phs Med Freq','Power','Energy Ratio','Delta Pow','Theta Pow',...
    'Alpha Pow','Beta Pow','Gamma Pow','H-Gamma Pow','Correlation','Delta Coh',...
    'Theta Coh','Alpha Coh','Beta Coh','Gamma Coh','H-Gamma Coh'};
participants_info=tdfread(['F:\RESEARCH\Hamid\Multicentre dataset\ds004100\participants.tsv']);
subjects_analysed=[1:25 27:39 41:45 47:58]; % already remove 40 becuase "data" does not have that
what=3;% 1=within patient; 2= across time; 3=across patients
ictal_or_inter='Ictal';
sorted=0;

if what==1
    load(['Within_subject_performance_sepnorm',ictal_or_inter,'_all_feats_crcted_feats_imp_100_10.mat'])
elseif what==2
    if strcmp(ictal_or_inter,'Interictal')
        load(['Generalisation_performance_across_time_all_feats_comb_crcted_feats_ovrsmp_imp_100_10.mat'])
    else
        load(['Generalisation_performance_across_time_ict_to_int_all_feats_comb_crcted_feats_ovrsmp_imp_100_10.mat'])
    end
elseif what==3
    load(['Generalisation_performance_across_subjects_',ictal_or_inter,'_all_feats_comb_crcted_feats_ovrsmp_imp_100_10.mat'])
end

ff=1;
iter_rand=1;
feat_imp=nan(1,34,1);
for fld=1:size(impCART,5)
    for iter_equalis=1
        for p=1:55
            feat_imp(p,included_feats{ff,iter_equalis,p,iter_rand},iter_equalis,fld)=impCART{ff,iter_equalis,p,iter_rand,fld};
        end
    end
end
feat_imp(feat_imp==0)=nan;
datas=nanmean(feat_imp,4);
data_plot=nan(1,34,55,2);
for condition=[1]
    for feat=1:34
        % BFs across conditions
        pt=0;
        for P=subjects_analysed
            pt=pt+1;
            if condition==1
                %                 if contains(participants_info.engel(P,1:2),'1A') || contains(participants_info.engel(P,1:2),'2A')
                %                     g{pt,1}='S';
                %                 else
                %                     g{pt,1}='F';
                %                 end
                %                 %             g{pt,1}=participants_info.engel(P,1);
                g{pt,1}=participants_info.outcome(P,1);
            elseif condition==2
                g{pt,1}=participants_info.target(P,:);
            elseif condition==3
                g{pt,1}=participants_info.lesion_status(P,:);
            elseif condition==4
                g{pt,1}=participants_info.implant(P,:);
            end
        end
        
        data_bf=datas(:,feat);
        nans=isnan(data_bf);
        nans([1])=1; % bad labeling of lesion (1)
        data_bf(nans)=[];
        g(nans)=[];
        uniq_conds_tmp=unique(g);
        
        combs_tmp=nchoosek(1:length(uniq_conds_tmp),2);
        for comb=1:size(combs_tmp,1)
            try
                conditions_BFtmp(comb,1)=bf.ttest2(data_bf(strcmp(g,uniq_conds_tmp(combs_tmp(comb,1)))),data_bf(strcmp(g,uniq_conds_tmp(combs_tmp(comb,2)))));
                dat_tmp=data_bf(strcmp(g,uniq_conds_tmp(combs_tmp(comb,1))));
                data_plot(comb,feat,1:length(dat_tmp),1)=dat_tmp;
                dat_tmp=data_bf(strcmp(g,uniq_conds_tmp(combs_tmp(comb,2))));
                data_plot(comb,feat,1:length(dat_tmp),2)=dat_tmp;
            catch
                conditions_BFtmp(comb,1)=nan;
            end
            combinations_tmp{comb,1}=uniq_conds_tmp(combs_tmp(comb,1));
            combinations_tmp{comb,2}=uniq_conds_tmp(combs_tmp(comb,2));
        end
        combinations_BF{condition,feat}=conditions_BFtmp;
        combinations{condition,feat}=combinations_tmp;
        clearvars g data_bf conditions_BFtmp combinations_tmp
    end
    
    figure
    c=0;
    for feat=1:34
        c=c+7;
        colours=[0.3 0.7 0.3];
        data_tmp2=squeeze(data_plot(1,feat,~isnan(data_plot(1,feat,:,2)),2));
        swarmchart([(c-1)*ones(sum(~isnan(data_tmp2)),1)]',data_tmp2,'MarkerFaceColor',colours,'MarkerEdgeColor',colours,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.3);
        hold on;
        colours=[0.7 0.3 0.3];
        data_tmp1=squeeze(data_plot(1,feat,~isnan(data_plot(1,feat,:,1)),1));
        swarmchart([(c+1)*ones(sum(~isnan(data_tmp1)),1)]',data_tmp1,'MarkerFaceColor',colours,'MarkerEdgeColor',colours,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.3);
        BFs=bf.ttest2(data_tmp1,data_tmp2);
        if what<3
            y_position=-0.45;
        elseif what==3
            y_position=1.05;
        end
        t=text(c,y_position,sprintf('%0.2f',BFs));
        t.Rotation=90;
        t.FontSize=14;
        if BFs>3
            t.FontWeight='Bold';
        end
        clearvars data_tmp1 data_tmp2
    end
    
    box off
    grid on
    set(gca,'TickDir','out','Fontsize',16)
    ylabel('Feature Contributon (A.U.)')
    feats=[1:8 10:12 14:36];
    xticks([7:7:7*length(feats)])
    xticklabels([Features_labels(feats)])
    
    if what==1
        ylim([-0.5 2.5])
        title(['Within-Patient Classification ',ictal_or_inter])
    elseif what==2
        ylim([-0.5 4])
        if strcmp(ictal_or_inter,'Interictal')
            title(['Cross-Time Generalisation Interictal to Ictal',])
        else
            title(['Cross-Time Generalisation Ictal to Interictal',])
        end
    elseif what==3
        ylim([1 8])
        title(['Cross-Patient Generalisation ',ictal_or_inter])
    end
end


