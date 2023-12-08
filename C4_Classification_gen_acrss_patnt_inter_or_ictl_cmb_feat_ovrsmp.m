% This code classifies resected and non-resected contacts, but does it by
% generalising across patients: trains the classifiers on all-minus-one
% patient and tests it on the left out patient using Decision Tree
% Classifiers

% INPUTS: features data separated by contacts from C3_Separating_target_non_target_contacts_all_feats
% OUTPUTS: classification data to be permuted by C5_Permuting_all_classifciation_windows
%%
clc
clear all
% loading the data from C3
ictal_or_inter='interictal';
load(['Effect_size_data_',ictal_or_inter,'_crcted_all_feats_100_to_10_points.mat'],'data_target_all','data_non_target_all','positive_negative_modulation','num_samples_keep')
%% Classification: Within ictal or interictal
all_patients=[1:39 41:56];
ff=0;
features_used=[1:34]; %% use all features
ff=ff+1;
for iter_equalis=1:1000 % how many iterations
    for p=1:length(all_patients)
        % Training and testing set indices
        all_patients_minus_one=all_patients;
        all_patients_minus_one(p)=[];
        
        Xtrain=[];
        ytrain=[];
        % Preparing the training set
        % Concatanating all-minus-one patient
        for Patient=all_patients_minus_one
            data_targ_res=[];
            data_non_res=[];
            f=0;
            for feats=features_used
                f=f+1;
                % Upsample the data from the class with lower number
                % of contacts (usually resected) to equalise them with
                % the class with higher number of contacts (usually non-resected)
                if size(data_non_target_all{Patient,feats},1)>=size(data_target_all{Patient,feats},1)
                    if strcmp(ictal_or_inter,'ictal') % ictal data
                        % Normalising the feature values in the post onset by mean of values in the pre-onset
                        targ=squeeze(([(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                        data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                        
                        nont=squeeze(([(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                        data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                        
                        samp=randsample([1:size(data_targ,1)],size(data_non,1)-size(data_targ,1),true);
                        data_targ=vertcat(data_targ,data_targ(samp,:));
                    else % interictal data
                        targ=squeeze(data_target_all{Patient,feats});
                        data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                        
                        nont=squeeze((data_non_target_all{Patient,feats}));
                        data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                        
                        samp=randsample([1:size(data_targ,1)],size(data_non,1)-size(data_targ,1),true);
                        data_targ=vertcat(data_targ,data_targ(samp,:));
                    end
                else
                    if strcmp(ictal_or_inter,'ictal') % ictal data
                        % Normalising the feature values in the post onset by mean of values in the pre-onset
                        nont=squeeze(([(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                        data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                        
                        targ=squeeze(([(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                        data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                        
                        samp=randsample([1:size(data_non,1)],size(data_targ,1)-size(data_non,1),true);
                        data_non=vertcat(data_non,data_non(samp,:));
                    else
                        nont=squeeze(data_non_target_all{Patient,feats});
                        data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                        
                        targ=squeeze((data_target_all{Patient,feats}));
                        data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                        samp=randsample([1:size(data_non,1)],size(data_targ,1)-size(data_non,1),true);
                        data_non=vertcat(data_non,data_non(samp,:));
                    end
                end
                data_targ_res(:,f)=reshape(data_targ,[size(data_targ,1)*size(data_targ,2)],[]);
                data_non_res(:,f)=reshape(data_non,[size(data_non,1)*size(data_non,2)],[]);
            end
            ytrain=vertcat(ytrain,[ones(size(data_targ_res,1),1);zeros(size(data_non_res,1),1)]);
            Xtrain=vertcat((Xtrain),[data_targ_res;data_non_res]);
        end
        
        % Preparing the testing set
        % Now do the same for the left-out patient
        % Similar to the above; with only one patinet; no concatanation
        ytest=[];
        Patient=all_patients(p);
        data_targ_res=[];
        data_non_res=[];
        f=0;
        for feats=features_used
            f=f+1;
            if size(data_non_target_all{Patient,feats},1)>=size(data_target_all{Patient,feats},1)
                if strcmp(ictal_or_inter,'ictal')
                    targ=squeeze(([(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                    data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                    
                    nont=squeeze(([(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                    data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                    
                    samp=randsample([1:size(data_targ,1)],size(data_non,1)-size(data_targ,1),true);
                    data_targ=vertcat(data_targ,data_targ(samp,:));
                else
                    targ=squeeze(data_target_all{Patient,feats});
                    data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                    
                    nont=squeeze((data_non_target_all{Patient,feats}));
                    data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                    
                    samp=randsample([1:size(data_targ,1)],size(data_non,1)-size(data_targ,1),true);
                    data_targ=vertcat(data_targ,data_targ(samp,:));
                end
            else
                if strcmp(ictal_or_inter,'ictal')
                    nont=squeeze(([(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                    data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                    
                    targ=squeeze(([(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                    data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                    
                    samp=randsample([1:size(data_non,1)],size(data_targ,1)-size(data_non,1),true);
                    data_non=vertcat(data_non,data_non(samp,:));
                else
                    nont=squeeze(data_non_target_all{Patient,feats});
                    data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                    
                    targ=squeeze((data_target_all{Patient,feats}));
                    data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                    
                    samp=randsample([1:size(data_non,1)],size(data_targ,1)-size(data_non,1),true);
                    data_non=vertcat(data_non,data_non(samp,:));
                end
            end
            data_targ_res(:,f)=reshape(data_targ,[size(data_targ,1)*size(data_targ,2)],[]);
            data_non_res(:,f)=reshape(data_non,[size(data_non,1)*size(data_non,2)],[]);
        end
        ytest=[ones(size(data_targ_res,1),1);zeros(size(data_non_res,1),1)];
        Xtest=([data_targ_res;data_non_res]);
        
        % Removing the nan features
        c=0;
        for i=1:size(Xtrain,2)
            if mean(~isnan(Xtrain(:,i)))>0.5 && mean(~isnan(Xtest(:,i)))>0.5
                c=c+1;
                Xnewtrain(:,c)=normalize(Xtrain(:,i));
                Xnewtest(:,c)=normalize(Xtest(:,i));
                included_feats{ff,iter_equalis,p}(c)=i;
            end
        end
        
        if exist('Xnewtrain','var') && exist('Xnewtest','var')
            Xtrain=(Xnewtrain);
            Xtest=(Xnewtest);
            
            % randomising class labels (to generate null distribution for statistical
            % testing)?: No, we do it in another file; so iter_rand=1
            for iter_rand=1:1
                if iter_rand~=1
                    ytrain_r=randsample(ytrain,length(ytrain));
                else
                    ytrain_r=ytrain;
                end
                
                % classify the contacts using the decision tree
                Classifier_Model = TreeBagger(30,Xtrain,ytrain_r,...
                    Method="classification",...
                    OOBPrediction="on",OOBPredictorImportance="on");
                impCART{ff,iter_equalis,p,iter_rand} = Classifier_Model.OOBPermutedPredictorDeltaError;
                preds_tmp=predict(Classifier_Model,Xtest);
                for i=1:length(preds_tmp)
                    preds(i,1)=str2num(preds_tmp{i});
                end
                Predictions{ff,iter_equalis,p,iter_rand,1}=preds;
                clearvars preds
                Classifiers={'RF'};
                
                Ground_truth{ff,iter_equalis,p,iter_rand}=ytest;
            end
        else
            Predictions{ff,iter_equalis,p,1,1}=nan;
            Ground_truth{ff,iter_equalis,p,1}=nan;
            Classifiers={'Any'};
        end
        [ff iter_equalis p]
        clearvars Xtest Xtrain ytrain ytrain_r ytest Xnewtrain Xnewtest
        % saving the classification results + the ground truths data
        save(['Generalisation_performance_across_subjects_',ictal_or_inter,'_all_feats_comb_crcted_feats_ovrsmp_imp_100_10_30bags.mat'],'Classifiers','Predictions','Ground_truth','impCART','included_feats','-v7.3')
    end
end
