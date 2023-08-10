%% Feature extraction
clc
clear all;
close all;
addpath(genpath('F:\Toolbox\eeglab2021.1'))
p=0;
for Patient=[1:39 41:56]
    p=p+1;
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
    elseif Patient==57
        Patient_initials='jh103';
        modality='ecog';
        ictal_trials=[1:3];
        inter_ictal_trials=[1:2];
    elseif Patient==58
        Patient_initials='jh105';
        ictal_trials=[1:5];
        modality='ecog';
        inter_ictal_trials=[1:2];
    elseif Patient==59
        Patient_initials='pt1';
        modality='ecog';
        ictal_trials=[1:4];
        inter_ictal_trials=[1:4];
    elseif Patient==60
        Patient_initials='pt2';
        ictal_trials=[1:3];
        modality='ecog';
        inter_ictal_trials=[1:4];
    elseif Patient==61
        Patient_initials='pt3';
        ictal_trials=[1:2];
        modality='ecog';
        inter_ictal_trials=[1:3];
    elseif Patient==62
        Patient_initials='umf001';
        ictal_trials=[1];
        modality='ecog';
        inter_ictal_trials=[1:2];
    end
    inter_ictal_trials_all(p)=max(inter_ictal_trials);
    ictal_trials_all(p)=max(ictal_trials);
    [Patient]
end
