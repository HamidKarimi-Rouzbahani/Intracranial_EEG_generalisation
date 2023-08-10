%% Plotting corrleations within patients
clc;
clear all;
close all;
ictal_or_inter='interictal';
load(['Within_subject_performance_',ictal_or_inter,'_all_feats_crcted_feats_imp_100_10.mat'])
ff=1;
iter_rand=1;
for fld=1:size(Ground_truth,5)
    for iter_equalis=1:size(Ground_truth,2)
        for p=1:size(Ground_truth,3)
            Performance_inter(ff,iter_equalis,p,iter_rand,fld,1:7)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand,fld},Predictions{ff,iter_equalis,p,iter_rand,fld});
            [~,~,~,Performance_inter(ff,iter_equalis,p,iter_rand,fld,8)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand,fld},Predictions{ff,iter_equalis,p,iter_rand,fld},1);
        end
    end
end
measure=8;
Performance_inter=squeeze(nanmean(nanmean(Performance_inter(:,:,:,:,measure),5),2));
for fld=1:size(impCART,5)
    for iter_equalis=1
        for p=1:55
            feat_imp_inter(p,included_feats{ff,iter_equalis,p,iter_rand},iter_equalis,fld)=impCART{ff,iter_equalis,p,iter_rand,fld};
        end
    end
end
feat_imp_inter(feat_imp_inter==0)=nan;
feat_imp_inter=squeeze(nanmean(feat_imp_inter,4));
ictal_or_inter='ictal';
load(['Within_subject_performance_',ictal_or_inter,'_all_feats_crcted_feats_imp_100_10.mat'])
ff=1;
iter_rand=1;
for fld=1:size(Ground_truth,5)
    for iter_equalis=1:size(Ground_truth,2)
        for p=1:size(Ground_truth,3)
            Performance_ictal(ff,iter_equalis,p,iter_rand,fld,1:7)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand,fld},Predictions{ff,iter_equalis,p,iter_rand,fld});
            [~,~,~,Performance_ictal(ff,iter_equalis,p,iter_rand,fld,8)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand,fld},Predictions{ff,iter_equalis,p,iter_rand,fld},1);
        end
    end
end
measure=8;
Performance_ictal=squeeze(nanmean(nanmean(Performance_ictal(:,:,:,:,:,measure),5),2));
for fld=1:size(impCART,5)
    for iter_equalis=1
        for p=1:55
            feat_imp_ictal(p,included_feats{ff,iter_equalis,p,iter_rand},iter_equalis,fld)=impCART{ff,iter_equalis,p,iter_rand,fld};
        end
    end
end

feat_imp_ictal(feat_imp_ictal==0)=nan;
feat_imp_ictal=squeeze(nanmean(feat_imp_ictal,4));

performances=rmoutliers([Performance_inter Performance_ictal]);
[cors,ps]=corr(performances(:,1),performances(:,2));
scatter(performances(:,1),performances(:,2),100,'MarkerEdgeColor',[0.3 0.3 0.95],'MarkerFaceColor',[0 0 0])
xlabel('Interictal Classification (AUC)')
ylabel('Ictal Classification (AUC)')
xlim([0.92 1])
ylim([0.92 1])
hold on;
plot([0.92 1],[0.92 1],'--k')
title({['r = ',sprintf('%.2f',cors), '; P = ',sprintf('%.2f',ps)]})
set(gca,'TickDir','out','Fontsize',20)

figure;

[cors2,ps2]=corr(nanmean(feat_imp_inter)',nanmean(feat_imp_ictal)');
scatter(nanmean(feat_imp_inter)',nanmean(feat_imp_ictal)',100,'MarkerEdgeColor',[0.3 0.3 0.95],'MarkerFaceColor',[0 0 0])
xlabel('Interictal Feature Contribution (A.U.)')
ylabel('Ictal Feature Contribution (A.U.)')
xlim([0.3 1])
ylim([0.3 1])
hold on;
plot([0.3 1],[0.3 1],'--k')
title({['r = ',sprintf('%.2f',cors2), '; P = ',sprintf('%.2f',ps2)]})
set(gca,'TickDir','out','Fontsize',20)

%% Plotting corrleations across time
 
clc;
clear all;
close all;
load(['Generalisation_performance_across_time_all_feats_comb_crcted_feats_ovrsmp_imp_100_10.mat'])

classif=1; %
ff=1;
iter_rand=1;
for iter_equalis=1:size(Ground_truth,2)
    for p=1:size(Ground_truth,3)
        Performance_inter(ff,iter_equalis,p,iter_rand,1:7)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif});
        [~,~,~,Performance_inter(ff,iter_equalis,p,iter_rand,8)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif},1);
    end
end
measure=8;

Performance_inter=squeeze(nanmean(Performance_inter(:,:,:,:,measure),2));
for iter_equalis=1
    for p=1:55
        feat_imp_inter(p,included_feats{ff,iter_equalis,p,iter_rand},iter_equalis)=impCART{ff,iter_equalis,p,iter_rand};
    end
end
feat_imp_inter(feat_imp_inter==0)=nan;

load(['Generalisation_performance_across_time_ict_to_int_all_feats_comb_crcted_feats_ovrsmp_imp_100_10.mat'])
endclassif=1; %
ff=1;
iter_rand=1;
for iter_equalis=1:size(Ground_truth,2)
    for p=1:size(Ground_truth,3)
        Performance_ictal(ff,iter_equalis,p,iter_rand,1:7)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif});
        [~,~,~,Performance_ictal(ff,iter_equalis,p,iter_rand,8)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif},1);
    end
end
measure=8;
Performance_ictal=squeeze(nanmean(Performance_ictal(:,:,:,:,measure),2));
for iter_equalis=1
    for p=1:55
        feat_imp_ictal(p,included_feats{ff,iter_equalis,p,iter_rand},iter_equalis)=impCART{ff,iter_equalis,p,iter_rand};
    end
end
feat_imp_ictal(feat_imp_ictal==0)=nan;

performances=rmoutliers([Performance_inter Performance_ictal]);
[cors,ps]=corr(performances(:,1),performances(:,2));
scatter(performances(:,1),performances(:,2),100,'MarkerEdgeColor',[0.3 0.3 0.95],'MarkerFaceColor',[0 0 0])
xlabel('Inter-to-Ict Generalisation (AUC)')
ylabel('Ict-to-Inter Generalisation (AUC)')
xlim([0.4 0.9])
ylim([0.4 0.9])
xticks([0.4:0.1:0.9])
yticks([0.4:0.1:0.9])
hold on;
plot([0.4 0.9],[0.4 0.9],'--k')
title({['r = ',sprintf('%.2f',cors), '; P = ',sprintf('%.2f',ps)]})
set(gca,'TickDir','out','Fontsize',20)

figure

[cors2,ps2]=corr(nanmean(feat_imp_inter)',nanmean(feat_imp_ictal)');
scatter(nanmean(feat_imp_inter)',nanmean(feat_imp_ictal)',100,'MarkerEdgeColor',[0.3 0.3 0.95],'MarkerFaceColor',[0 0 0])
xlabel('Inter-to-Ict Feature Contribution (A.U.)')
ylabel('Ict-to-Inter Feature Contribution (A.U.)')
xlim([0.5 2])
ylim([0.5 2])
xticks([0.5:0.5:2])
yticks([0.5:0.5:2])
hold on;
plot([0.5 2],[0.5 2],'--k')
title({['r = ',sprintf('%.2f',cors2), '; P = ',sprintf('%.2f',ps2)]})
set(gca,'TickDir','out','Fontsize',20)

%% Plotting corrleations across patients
clc;
clear all;
close all;
ictal_or_inter='interictal';
load(['Generalisation_performance_across_subjects_',ictal_or_inter,'_all_feats_comb_crcted_feats_ovrsmp_imp_100_10.mat'])
classif=1; %
ff=1;
iter_rand=1;
for iter_equalis=1:size(Ground_truth,2)
    for p=1:size(Ground_truth,3)
        Performance_inter(ff,iter_equalis,p,iter_rand,1:7)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif});
        [~,~,~,Performance_inter(ff,iter_equalis,p,iter_rand,8)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif},1);
    end
end
measure=8;

Performance_inter=squeeze(nanmean(Performance_inter(:,:,:,:,measure),2));
for iter_equalis=1
    for p=1:55
        feat_imp_inter(p,included_feats{ff,iter_equalis,p,iter_rand},iter_equalis)=impCART{ff,iter_equalis,p,iter_rand};
    end
end
feat_imp_inter(feat_imp_inter==0)=nan;

ictal_or_inter='ictal';
load(['Generalisation_performance_across_subjects_',ictal_or_inter,'_all_feats_comb_crcted_feats_ovrsmp_imp_100_10.mat'])
classif=1; %
ff=1;
iter_rand=1;
for iter_equalis=1:size(Ground_truth,2)
    for p=1:size(Ground_truth,3)
        Performance_ictal(ff,iter_equalis,p,iter_rand,1:7)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif});
        [~,~,~,Performance_ictal(ff,iter_equalis,p,iter_rand,8)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif},1);
    end
end
measure=8;
Performance_ictal=squeeze(nanmean(Performance_ictal(:,:,:,:,measure),2));
for iter_equalis=1
    for p=1:55
        feat_imp_ictal(p,included_feats{ff,iter_equalis,p,iter_rand},iter_equalis)=impCART{ff,iter_equalis,p,iter_rand};
    end
end
feat_imp_ictal(feat_imp_ictal==0)=nan;

performances=rmoutliers([Performance_inter Performance_ictal]);
[cors,ps]=corr(performances(:,1),performances(:,2));
scatter(performances(:,1),performances(:,2),100,'MarkerEdgeColor',[0.3 0.3 0.95],'MarkerFaceColor',[0 0 0])
xlabel('Interictal Classification (AUC)')
ylabel('Ictal Classification (AUC)')
xlim([0.3 0.9])
ylim([0.3 0.9])
hold on;
plot([0.3 0.9],[0.3 0.9],'--k')
title({['r = ',sprintf('%.2f',cors), '; P = ',sprintf('%.2f',ps)]})
set(gca,'TickDir','out','Fontsize',20)

figure

[cors2,ps2]=corr(nanmean(feat_imp_inter)',nanmean(feat_imp_ictal)');
scatter(nanmean(feat_imp_inter)',nanmean(feat_imp_ictal)',100,'MarkerEdgeColor',[0.3 0.3 0.95],'MarkerFaceColor',[0 0 0])
xlabel('Interictal Feature Contribution (A.U.)')
ylabel('Ictal Feature Contribution (A.U.)')
xlim([1.5 6])
ylim([1.5 6])
hold on;
plot([1.5 6],[1.5 6],'--k')
title({['r = ',sprintf('%.2f',cors2), '; P = ',sprintf('%.2f',ps2)]})
set(gca,'TickDir','out','Fontsize',20)



