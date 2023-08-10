clc;
clear all;
close all;
chars='_less_time_res'; %'' '_less_time_res' '_more_windows'
Features_labels={'Mean','Median','Variance','Skewness','Kurtosis',...
    'LZ Comp','Higichi FD','Katz FD','Lyap Exp','Hurst Exp',...
    'Samp Ent','Apprx Ent','Autocor','Hjorth Comp','Hjorth Mob',...
    'Mean Freq','Med Freq','Avg Freq','SEF','Pow Med Freq',...
    'Phs Med Freq','Power','Energy ratio'};
addpath(genpath('F:\RESEARCH\Hamid\Matlab\bayesFactor-master'))
Patient=4;
also_show_behaviour=0; % also behaviour?
ictal_or_inter='ictal'; % 'ictal' or 'interictal'
if strcmp(ictal_or_inter,'ictal')
    pre_siz=10000;  %*2ms
    post_siz=10000; %*2ms
    window_spans=[10 100 500]; %*2ms
    step_size=100; % *2ms
    min_max_baseline=[-10000 -5000];
elseif strcmp(ictal_or_inter,'interictal')
    pre_siz=1;  %*2ms
    post_siz=30000; %*2ms
    window_spans=[100 500 1000]; %*2ms
    step_size=100; % *2ms
end

shifts_to_onset=[0 0 0 0 0]; %EEG onset: ms

if Patient==1
    Patient_initials='060';
    modality='seeg';
    ictal_trials=[1:3];
    inter_ictal_trials=[1:2];
elseif Patient==2
    Patient_initials='064';
    ictal_trials=[1];
    modality='ecog';
    inter_ictal_trials=[2];
elseif Patient==3
    Patient_initials='065';
    modality='ecog';
    ictal_trials=[1:3];
    inter_ictal_trials=[1:2];
elseif Patient==4
    Patient_initials='070';
    ictal_trials=[1:3];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==5
    Patient_initials='074';
    ictal_trials=[1:3];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==6
    Patient_initials='075';
    ictal_trials=[1];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==7
    Patient_initials='080';
    ictal_trials=[1:4];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==8
    Patient_initials='082';
    ictal_trials=[1:5];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==9
    Patient_initials='086';
    ictal_trials=[1:2];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==10
    Patient_initials='087';
    ictal_trials=[1:2];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==11
    Patient_initials='088';
    ictal_trials=[1:3];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==12
    Patient_initials='089';
    ictal_trials=[1:4];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==13
    Patient_initials='094';
    ictal_trials=[1:3];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==14
    Patient_initials='097';
    ictal_trials=[1:5];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==15
    Patient_initials='105';
    ictal_trials=[1:2];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==16
    Patient_initials='106';
    ictal_trials=[1:3];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==17
    Patient_initials='107';
    ictal_trials=[1:5];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==18
    Patient_initials='111';
    ictal_trials=[1:5];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==19
    Patient_initials='112';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==20
    Patient_initials='114';
    ictal_trials=[1:4];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==21
    Patient_initials='116';
    ictal_trials=[1:3];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==22
    Patient_initials='117';
    ictal_trials=[1:3];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==23
    Patient_initials='123';
    ictal_trials=[1:4];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==24
    Patient_initials='126';
    ictal_trials=[1:4];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==25
    Patient_initials='130';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==26
    Patient_initials='133';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==27
    Patient_initials='134';
    ictal_trials=[1];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==28
    Patient_initials='135';
    ictal_trials=[1:2];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==29
    Patient_initials='138';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==30
    Patient_initials='139';
    ictal_trials=[1:3];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==31
    Patient_initials='140';
    ictal_trials=[1:3];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==32
    Patient_initials='141';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==33
    Patient_initials='142';
    ictal_trials=[1:3];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==34
    Patient_initials='144';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==35
    Patient_initials='146';
    ictal_trials=[1:3];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==36
    Patient_initials='148';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==37
    Patient_initials='150';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==38
    Patient_initials='151';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==39
    Patient_initials='157';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==40
    Patient_initials='158';
    ictal_trials=[4];
    modality='seeg';
    inter_ictal_trials=[];
elseif Patient==41
    Patient_initials='160';
    ictal_trials=[1:3];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==42
    Patient_initials='162';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==43
    Patient_initials='163';
    ictal_trials=[1:3];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==44
    Patient_initials='164';
    ictal_trials=[1:3];
    modality='seeg';
    inter_ictal_trials=[1:2];
    %     elseif Patient==45
    %         Patient_initials='165';
    %         seizures=[];
    %         modality='seeg';
    %         inter_ictal_trials=[1:2];
elseif Patient==45
    Patient_initials='166';
    ictal_trials=[1:2];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==46
    Patient_initials='171';
    ictal_trials=[1:4];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==47
    Patient_initials='172';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==48
    Patient_initials='173';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==49
    Patient_initials='177';
    ictal_trials=[1:3];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==50
    Patient_initials='179';
    ictal_trials=[1:2];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==51
    Patient_initials='180';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==52
    Patient_initials='181';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==53
    Patient_initials='185';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==54
    Patient_initials='187';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==55
    Patient_initials='188';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==56
    Patient_initials='190';
    ictal_trials=[1:3];
    modality='seeg';
    inter_ictal_trials=[1:2];
end


times=[-pre_siz:step_size:post_siz-max(window_spans)-1];
baseline_span=[find(times==min_max_baseline(1)):find(times==min_max_baseline(2))];
if also_show_behaviour==1
    Behavioural_timings=[behavioural_symps*(1000/step_size)]';
    Behav_tims=times([find(times==0)-1+Behavioural_timings(~isnan(Behavioural_timings))']);
end

clearvars tmp_data tmp_data_baseline
if strcmp(ictal_or_inter,'ictal')
%     load([Patient_initials,'_ds4100_Features_Seizure',num2str(1),chars,'.mat'])
    load([Patient_initials,'_All_Features_Seizure',num2str(1),chars,'.mat'])
elseif strcmp(ictal_or_inter,'interictal')
    load([Patient_initials,'_ds4100_Features_Interictal',num2str(1),chars,'.mat'])
end
chans_feats=repmat(1:length(channels),[length(featuress) 1]);


%% Plotting
if strcmp(ictal_or_inter,'ictal')
    trials=ictal_trials;
elseif strcmp(ictal_or_inter,'interictal')
    trials=inter_ictal_trials;
end
participants_info=tdfread(['F:\RESEARCH\Hamid\Multicentre dataset\ds004100\participants.tsv']);
% Plotting all features
ticks=[-20000 -15000 -10000 -5000 0 5000 10000 15000 20000];
subjects_analysed=[1:25 27:45 47:58];
w=[1:3];
for w=3
    figure;
    f=0;
    for feat=featuress
        f=f+1;
        tr=0;
        for trial=trials
            tr=tr+1;
            if strcmp(ictal_or_inter,'ictal')
                load([Patient_initials,'_ds4100_Features_Seizure',num2str(trial),chars,'.mat'])
            elseif strcmp(ictal_or_inter,'interictal')
                load([Patient_initials,'_ds4100_Features_Interictal',num2str(trial),chars,'.mat'])
            end
            patient_data_address=['F:\RESEARCH\Hamid\Multicentre dataset\ds004100\sub-HUP',Patient_initials,'\ses-presurgery\ieeg\'];
            folder_address=[patient_data_address,'sub-HUP',Patient_initials,'_ses-presurgery_task-',ictal_or_inter,'_acq-',modality,'_run-',sprintf('%02d',trial)];
            channel_info=tdfread([folder_address,'_channels.tsv']);
            c=0;
            for ch=1:size(channel_info.status_description,1)
                if contains(channel_info.status_description(ch,:),'resect') && ~isempty(find(contains(channels_id,strtrim(channel_info.name(ch,:)))))
                    c=c+1;
                    Specific_targ_chans{c}=channels_id{1,find(contains(channels_id,strtrim(channel_info.name(ch,:))))};
                end
            end
            chans_tmp=[];
            for ch=1:size(Specific_targ_chans,2)
                Actual_particip_id=find(strcmp(channels_id,Specific_targ_chans{ch}));
                chans_tmp=horzcat(chans_tmp,chans_feats(:,Actual_particip_id));
            end
            chans_feats_suspected=chans_tmp;
            included_chans=channels;
            included_chans(chans_tmp(1,:))=[];
            chans_feats_nonsuspected=chans_feats(:,included_chans);
            
            tmp_data_Suspected(:,:,tr)=circshift(squeeze(nanmean(features(feat,chans_feats_suspected(f,:),:,w),4))',-shifts_to_onset(trial))';
            tmp_data_baselineSuspected(:,:,tr)=circshift(squeeze(nanmean(features(feat,chans_feats_suspected(f,:),baseline_span,w),4))',-shifts_to_onset(trial))';
            tmp_data_NonSuspected(:,:,tr)=circshift(squeeze(nanmean(features(feat,chans_feats_nonsuspected(f,:),:,w),4))',-shifts_to_onset(trial))';
            tmp_data_baselineNonSuspected(:,:,tr)=circshift(squeeze(nanmean(features(feat,chans_feats_nonsuspected(f,:),baseline_span,w),4))',-shifts_to_onset(trial))';
        end
        tmp_data_Suspected=squeeze(nanmean(tmp_data_Suspected,3));
        tmp_data_baseline_Suspected=squeeze(nanmean(tmp_data_baselineSuspected,3));
        subplot(5,5,f)
        Suspected=shadedErrorBar(times,nanmean(tmp_data_Suspected),nanstd(tmp_data_Suspected)./sqrt(size(tmp_data_Suspected,1)),'lineprops',{'-r','markerfacecolor',[1,0.2,0.2]},'transparent',1);
        hold on;
        plot([times(1) times(end)],[nanmean(nanmean(tmp_data_baseline_Suspected)) nanmean(nanmean(tmp_data_baseline_Suspected))],'--r');
        plot([0 0],[nanmin(nanmean(tmp_data_Suspected)) nanmax(nanmean(tmp_data_Suspected))],'--k');
        
        tmp_data_allNonSuspected=squeeze(nanmean(tmp_data_NonSuspected,3));
        tmp_data_baseline_NonSuspected=squeeze(nanmean(tmp_data_baselineNonSuspected,3));
        NonSuspected=shadedErrorBar(times,nanmean(tmp_data_allNonSuspected),nanstd(tmp_data_allNonSuspected)./sqrt(size(tmp_data_allNonSuspected,1)),'lineprops',{'-b','markerfacecolor',[1,0.2,0.2]},'transparent',1);
        plot([times(1) times(end)],[nanmean(nanmean(tmp_data_baseline_NonSuspected)) nanmean(nanmean(tmp_data_baseline_NonSuspected))],'--b');
        data_for_margins=squeeze(nanmean([tmp_data_Suspected;tmp_data_allNonSuspected],2));
        if also_show_behaviour==1
            beh_line=fill([nanmin(Behav_tims) nanmin(Behav_tims) nanmax(Behav_tims) nanmax(Behav_tims)],[nanmin(data_for_margins) nanmax(data_for_margins) nanmax(data_for_margins) nanmin(data_for_margins)],'g','FaceAlpha',.3);
        end
        
        %     for t=1:size(tmp_data_Suspected,2)
        %         p(t)=bf.ttest2(tmp_data_Suspected(:,t),tmp_data_NonSuspected(:,t));
        %         %         p(t)=ttest2(tmp_data_Suspected(:,t),tmp_data_NonSuspected(:,t));
        %     end
        %     significance_threshold=10;
        %     p(p<significance_threshold)=nan;
        %     p(p>significance_threshold)=1;
        %     data_max=nanmin(data_for_margins);
        %     plot(times,p*data_max*0.98,'*k','markersize',3)
        
        xlabel('Time relative to onset [s]')
        xticks(ticks)
        xticklabels(num2str([ticks/1000]'))
        ylabel('Value [AU]')
        title([Features_labels(feat)])
        grid on;
        
        if f==1 && also_show_behaviour==1
            legend([Suspected.mainLine NonSuspected.mainLine beh_line],{'Resected','Other','Behaviour'},'location','northwest')
        elseif f==1 && also_show_behaviour==0
            legend([Suspected.mainLine NonSuspected.mainLine],{'Resected','Other'},'location','northwest')
        end
    end
%     sgtitle(['Patient #',num2str(Patient),'; Window = ',num2str(window_spans(w)*2),'ms; ','Outcome = Engel ',participants_info.engel(Patient,:)])

    for i=1:size(participants_info.participant_id,1)
        if contains(participants_info.participant_id(i,:),Patient_initials)        
           Actual_particip_id=i; 
        end
    end
    sgtitle(['Patient #',num2str(Patient),'; Outcome: Engel ',participants_info.engel(Actual_particip_id,:)])
end
