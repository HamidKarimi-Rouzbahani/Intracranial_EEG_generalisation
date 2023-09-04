% This code classifies resected and non-resected contacts within each patient
% using Decision Tree Classifiers; It does it for an early, middle and late
% time window of data separately

% INPUTS: features data separated by contacts from C3_Separating_target_non_target_contacts_all_feats
% OUTPUTS: classification data to be permuted by C5_Permuting_3_classifciation_windows
%%
clc
clear all
%% Classification: Within ictal or interictal
% loading the data from C3
ictal_or_inter='interictal';
% time windows to analyse
if strcmp(ictal_or_inter,'interictal')
    time_windows=[1 14 28];
else
    time_windows=[1 7 14];
end
load(['Effect_size_data_',ictal_or_inter,'_crcted_all_feats_100_to_10_points.mat'],'data_target_all','data_non_target_all','positive_negative_modulation','num_samples_keep')
all_patients=[1:39 41:56]; % patient indices
ff=0;
features_used=[1:34]; % features
ff=ff+1;
tw=0;
for time_wind=time_windows
    tw=tw+1;
    for iter_equalis=1:10
        p=0;
        for Patient=all_patients
            % classificaiton is performed for each patient separately
            p=p+1;
            data_targ_res=[];
            data_non_res=[];
            f=0;
            num_samples_keep=time_wind; % how many samples per trial
            for feats=features_used
                f=f+1;
                % Downsample the data from the class with higher number
                % of contacts (usually non-resected) to equalise them with
                % the class with lower number of contacts (usually resected)
                if size(data_non_target_all{Patient,feats},1)>=size(data_target_all{Patient,feats},1)
                    if strcmp(ictal_or_inter,'ictal') % ictal data
                        % Normalising the feature values in the post onset by mean of values in the pre-onset
                        targ=squeeze(([(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                        data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                        
                        nont=squeeze(([(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                        data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                        
                        data_non=data_non(randsample([1:size(data_non,1)],size(data_targ,1)),:);
                    else % interictal data
                        targ=data_target_all{Patient,feats};
                        data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                        
                        nont=data_non_target_all{Patient,feats};
                        data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                        
                        data_non=data_non(randsample([1:size(data_non,1)],size(data_targ,1)),:);
                    end
                else
                    if strcmp(ictal_or_inter,'ictal')
                        targ=squeeze(([(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                        data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                        
                        nont=squeeze(([(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                        data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                        
                        data_targ=data_targ(randsample([1:size(data_targ,1)],size(data_non,1)),:);
                    else % interictal data
                        targ=data_target_all{Patient,feats};
                        data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                        
                        nont=data_non_target_all{Patient,feats};
                        data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                        
                        data_targ=data_targ(randsample([1:size(data_targ,1)],size(data_non,1)),:);
                    end
                end
                data_targ_res(:,f)=reshape(data_targ,[size(data_targ,1)*size(data_targ,2)],[]);
                data_non_res(:,f)=reshape(data_non,[size(data_non,1)*size(data_non,2)],[]);
            end
            X=[data_targ_res;data_non_res];
            y=[ones(size(data_targ_res,1),1);zeros(size(data_non_res,1),1)];
            c=0;
            for i=1:size(X,2)
                if mean(~isnan(X(:,i)))>0.5
                    c=c+1;
                    Xnew(:,c)=(X(:,i));
                    included_feats{ff,iter_equalis,p}(c)=i;
                end
            end
            if exist('Xnew','var')
                X=(Xnew);
                
                % randomising class labels (to generate null distribution for statistical
                % testing)?: No, we do it in another file; so iter_rand=1                
                for iter_rand=1:1
                    if iter_rand~=1
                        y_r=randsample(y,length(y));
                    else
                        y_r=y;
                    end
                    % classify the contacts using the decision tree
                    CVO = cvpartition(y_r,'k',10);
                    for fld = 1:CVO.NumTestSets
                        trIdx = CVO.training(fld);
                        teIdx = CVO.test(fld);
                        % normalising the data within the training ad testing
                        % sets separately to avoid circular analysis
                        Classifier_Model = TreeBagger(50,normalize(X(trIdx,:)),y_r(trIdx),...
                            Method="classification",...
                            OOBPrediction="on",OOBPredictorImportance="on");
                        impCART{tw,ff,iter_equalis,p,iter_rand,fld} = Classifier_Model.OOBPermutedPredictorDeltaError;
                        preds_tmp=predict(Classifier_Model,normalize(X(teIdx,:)));
                        for i=1:length(preds_tmp)
                            preds(i,1)=str2num(preds_tmp{i});
                        end
                        Predictions{tw,ff,iter_equalis,p,iter_rand,fld}=preds;
                        clearvars preds
                        Classifiers={'RF'};                        
                        Ground_truth{tw,ff,iter_equalis,p,iter_rand,fld}=y_r(teIdx);                        
                    end
                end
            else
                Predictions{tw,ff,iter_equalis,p,1,1,1}=nan;
                Ground_truth{tw,ff,iter_equalis,p,1,1}=nan;
                Classifiers={'Any'};
            end
            % saving the classification results + the ground truths data
            save(['Within_subject_performance_sepnorm_1_',ictal_or_inter,'_all_feats_crcted_feats_imp_100_10.mat'],'Classifiers','Predictions','Ground_truth','impCART','included_feats','-v7.3')
            [tw ff iter_equalis p]
            clearvars Xnew y_r X
        end
    end
end
