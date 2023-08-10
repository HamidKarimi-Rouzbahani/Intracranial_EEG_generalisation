%% Checking the effect of correlations
clc
clear all
ictal_or_inter='interictal';
data_to_plot='Normal';%'Absolute' or 'Normal'
differential_or_baselined_ictal=1;% 1=baselined; 2=not baselined

load(['Effect_size_data_',ictal_or_inter,'.mat'],'data_target_all','data_non_target_all')

for feats=1:21
    for Patient=[1:39 41:56]
        clearvars target abs_target nontarget abs_nontarget
        data_all{Patient,feats}=vertcat(data_target_all{Patient,feats},data_non_target_all{Patient,feats});
    end
end
iterations=101;
for iteration=1:iterations
    for feats=1:21
        for Patient=[1:39 41:56]
            clearvars target abs_target nontarget abs_nontarget
            samps_left=[1:size(data_all{Patient,feats},1)];
            if iteration==1
                samps=1:size(data_target_all{Patient,feats},1);
            else
                samps=randsample(samps_left,size(data_target_all{Patient,feats},1));
            end
            samps_left(samps)=[];
            data_target_all_new{Patient,feats}=data_all{Patient,feats}(samps,:,:);
            data_non_target_all_new{Patient,feats}=data_all{Patient,feats}(samps_left,:,:);
            
            for trial=1:size(data_target_all_new{Patient,feats},2)
                if strcmp(ictal_or_inter,'ictal')
                    abs_target(trial,:)=abs([(data_target_all_new{Patient,feats}(:,trial,2))-(data_target_all_new{Patient,feats}(:,trial,1))])./abs([(data_target_all_new{Patient,feats}(:,trial,2))+(data_target_all_new{Patient,feats}(:,trial,1))]);
                    if differential_or_baselined_ictal==1
                        target(trial,:)=([(data_target_all_new{Patient,feats}(:,trial,2))-(data_target_all_new{Patient,feats}(:,trial,1))])./abs([(data_target_all_new{Patient,feats}(:,trial,2))+(data_target_all_new{Patient,feats}(:,trial,1))]);
                    else
                        target(trial,:)=(data_target_all_new{Patient,feats}(:,trial,2));
                    end
                elseif strcmp(ictal_or_inter,'interictal')
                    target(trial,:)=(data_target_all_new{Patient,feats}(:,trial,1));
                end
            end
            for trial=1:size(data_non_target_all_new{Patient,feats},2)
                if strcmp(ictal_or_inter,'ictal')
                    abs_nontarget(trial,:)=abs([(data_non_target_all_new{Patient,feats}(:,trial,2))-(data_non_target_all_new{Patient,feats}(:,trial,1))])./abs([(data_non_target_all_new{Patient,feats}(:,trial,2))+(data_non_target_all_new{Patient,feats}(:,trial,1))]);
                    if differential_or_baselined_ictal==1
                        nontarget(trial,:)=([(data_non_target_all_new{Patient,feats}(:,trial,2))-(data_non_target_all_new{Patient,feats}(:,trial,1))])./abs([(data_non_target_all_new{Patient,feats}(:,trial,2))+(data_non_target_all_new{Patient,feats}(:,trial,1))]);
                    else
                        nontarget(trial,:)=(data_non_target_all_new{Patient,feats}(:,trial,2));
                    end
                elseif strcmp(ictal_or_inter,'interictal')
                    nontarget(trial,:)=(data_non_target_all_new{Patient,feats}(:,trial,1));
                end
            end
            if strcmp(data_to_plot,'Absolute')
                abs_target=nanmean(abs_target)';
                abs_nontarget=nanmean(abs_nontarget)';
                tmp=mes(abs_target,abs_nontarget,'hedgesg');
                data(feats,Patient)=abs(tmp.hedgesg);
            else
                target=nanmean(target)';
                nontarget=nanmean(nontarget)';
                tmp=mes(target,nontarget,'hedgesg');
                data(feats,Patient)=abs(tmp.hedgesg);
            end
        end
    end
    [iteration]
    data_iterations(:,:,iteration)=data;
end
save(['Data_for_correlation_',ictal_or_inter,'.mat'],'data_iterations');
cc
%% Correlation over all features
clc
clear all;
Features_labels={'Mean','Median','Variance','Skewness','Kurtosis',...
    'LZ Comp','Higichi FD','Katz FD','Lyap Exp','Hurst Exp',...
    'Samp Ent','Apprx Ent','Autocor','Hjorth Comp','Hjorth Mob',...
    'Mean Freq','Med Freq','Avg Freq','SEF','Pow Med Freq',...
    'Phs Med Freq','Power','Energy ratio'};
featuress=[1:8 10:12 14:23];
correlation_type=1; % 1 = over all features subjects avrgd; 2= on indiv feature over all subjects
iterations=101;
load(['Data_for_correlation_ictal.mat']);
data_iterations1=data_iterations;
load(['Data_for_correlation_interictal.mat']);
data_iterations2=data_iterations;
condition=ones(size(data_iterations1,2),1);
condition=logical(condition);
condition(40)=false; % unused

for iteration=1:iterations    
    [corrs(iteration),~]=corr(nanmean(data_iterations1(:,condition,iteration),2),nanmean(data_iterations2(:,condition,iteration),2));
    for i=1:21
        data_tmp=[data_iterations1(i,condition,iteration)' data_iterations2(i,condition,iteration)'];
        data_tmp=rmoutliers(data_tmp);     
        if iteration==1
            size(data_tmp,1)
        end
        [corrs_feat(i,iteration),~]=corr(data_tmp(:,1),data_tmp(:,2),'rows','complete');
    end
    [corrs_all(iteration),~]=corr(reshape(data_iterations1(:,condition,iteration),[21*55],[]),reshape(data_iterations2(:,condition,iteration),[21*55],[]),'rows','complete');
end

if correlation_type==1
    c=scatter(nanmean(data_iterations1(:,condition,1),2),nanmean(data_iterations2(:,condition,1),2),100,'filled');
    title(['Correlation over all features'])
    p=mean(corrs(1)<corrs(2:end));
    legend([c],['r=',sprintf('%.2f',corrs(1)),'; p=',sprintf('%.2f',p)],'location','northwest');
    xlabel('Ictal')
    ylabel('Interictal')
    title(['Correlation between Ictal and Interictal effects; Windows = 1000 and 2000 ms; ','Resected'])
elseif correlation_type==2
    for i=1:21
        subplot(5,5,i)
        c=scatter(data_iterations1(i,condition,1),data_iterations2(i,condition,2),50,'filled');
        p=mean(corrs_feat(i,1)<corrs_feat(i,2:end));        
        legend([c],['r=',sprintf('%.2f',corrs_feat(i,1)),'; p=',sprintf('%.2f',p)],'location','northwest');
        title([Features_labels(featuress(i))])
        xlabel('Ictal')
        ylabel('Interictal')
    end
    sgtitle(['Correlation between Ictal and Interictal effects; Windows = 1000 and 2000 ms; ','Resected'])
elseif correlation_type==3
    c=scatter(reshape(data_iterations1(:,condition,1),[21*55],[]),reshape(data_iterations2(:,condition,1),[21*55],[]),50,'filled');
end

%% BFs on Interictal-ictal correlations over all features
clc
clear all
Features_labels={'Mean','Median','Variance','Skewness','Kurtosis',...
    'LZ Comp','Higichi FD','Katz FD','Lyap Exp','Hurst Exp',...
    'Samp Ent','Apprx Ent','Autocor','Hjorth Comp','Hjorth Mob',...
    'Mean Freq','Med Freq','Avg Freq','SEF','Pow Med Freq',...
    'Phs Med Freq','Power','Energy ratio'};
featuress=[1:8 10:12 14:23];
load(['Data_for_correlation_ictal.mat']);
Effect_1=data_iterations(:,:,1);
load(['Data_for_correlation_interictal.mat']);
Effect_2=data_iterations(:,:,1);
participants_info=tdfread(['F:\RESEARCH\Hamid\Multicentre dataset\ds004100\participants.tsv']);

condition=1; %1= success/failure; 2= temporal roi/extra temporal; 3= lesional/non
subjects_analysed=[1:25 27:45 47:58];

pt=0;
for P=subjects_analysed
    pt=pt+1;
    g1{pt,1}=participants_info.engel(P,1);
    g2{pt,1}=participants_info.target(P,:);
    g3{pt,1}=participants_info.lesion_status(P,:);
    g4{pt,1}=participants_info.implant(P,:);
    y(pt)=corr(nanmean(Effect_1(:,pt),2),nanmean(Effect_2(:,pt),2));
end

clearvars data_ready
nans=isnan(y);
nans([1 40])=1; % bad labeling of lesion (1) and excluded patient (40)
data_ready=table(y',g1,g2,g3,g4);
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
bar(effects,0.2,'k');
xticks([1:length(effects)])
xticklabels({'Outcome','ROI','Lesion','Recording'})
xtickangle(45)
ylabel('Bayes factor')
grid on;
title(['BF ANOVA on Interictal-ictal correlations over all features; Windows = 1000 and 2000 ms; ','Resected'])

%% Single feature correlation
clc
clear all
Features_labels={'Mean','Median','Variance','Skewness','Kurtosis',...
    'LZ Comp','Higichi FD','Katz FD','Lyap Exp','Hurst Exp',...
    'Samp Ent','Apprx Ent','Autocor','Hjorth Comp','Hjorth Mob',...
    'Mean Freq','Med Freq','Avg Freq','SEF','Pow Med Freq',...
    'Phs Med Freq','Power','Energy ratio'};
featuress=[1:8 10:12 14:23];
load(['Data_for_correlation_ictal.mat']);
difference_1=data_iterations(:,:,1);
load(['Data_for_correlation_interictal.mat']);
difference_2=data_iterations(:,:,1);


artf_removal=0;
if artf_removal==1
    difference_1n=nan(21,56);
    difference_2n=nan(21,56);
    
    for feats=[1:21]
        differences=rmoutliers([difference_1(feats,:);difference_2(feats,:)]');
        difference_1n(feats,1:length(differences(:,1)))=differences(:,1);
        difference_2n(feats,1:length(differences(:,2)))=differences(:,2);
    end
else
    difference_1n=difference_1;
    difference_2n=difference_2;
end

for feats=[1:21]
    subplot(5,5,feats)
    c=scatter(difference_1n(feats,:),difference_2n(feats,:));
    [cors(feats) ps(feats)]=corr(difference_1n(feats,:)',difference_2n(feats,:)','row','complete');
    legend([c],['r=',sprintf('%.2f',cors(feats)),'; p=',sprintf('%.2f',ps(feats))]);
    xlabel('Ictal')
    ylabel('Interictal')
    title([Features_labels(featuress(feats))])
end
sgtitle(['Correlation between Ictal and Interictal effects; Windows = 1000 and 2000 ms; ','Resected'])
