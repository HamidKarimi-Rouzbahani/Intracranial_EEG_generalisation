%% Feature extraction
% clc
% clear all;
% close all;
function Data_loading_mine(Patient)
% for Patient=1:56
    if Patient==1
        Patient_initials='060';
        modality='seeg';
        seizures=[1:3];
    elseif Patient==2
        Patient_initials='064';
        seizures=[1];
        modality='ecog';
    elseif Patient==3
        Patient_initials='065';
        modality='ecog';
        seizures=[1:3];
    elseif Patient==4
        Patient_initials='070';
        seizures=[1:3];
        modality='ecog';
    elseif Patient==5
        Patient_initials='074';
        seizures=[1:3];
        modality='ecog';
    elseif Patient==6
        Patient_initials='075';
        seizures=[1];
        modality='ecog';
    elseif Patient==7
        Patient_initials='080';
        seizures=[1:4];
        modality='ecog';
    elseif Patient==8
        Patient_initials='082';
        seizures=[1:5];
        modality='ecog';
    elseif Patient==9
        Patient_initials='086';
        seizures=[1:2];
        modality='ecog';
    elseif Patient==10
        Patient_initials='087';
        seizures=[1:2];
        modality='ecog';
    elseif Patient==11
        Patient_initials='088';
        seizures=[1:3];
        modality='ecog';
    elseif Patient==12
        Patient_initials='089';
        seizures=[1:4];
        modality='ecog';
    elseif Patient==13
        Patient_initials='094';
        seizures=[1:3];
        modality='ecog';
    elseif Patient==14
        Patient_initials='097';
        seizures=[1:5];
        modality='ecog';
    elseif Patient==15
        Patient_initials='105';
        seizures=[1:2];
        modality='ecog';
    elseif Patient==16
        Patient_initials='106';
        seizures=[1:3];
        modality='ecog';
    elseif Patient==17
        Patient_initials='107';
        seizures=[1:5];
        modality='ecog';
    elseif Patient==18
        Patient_initials='111';
        seizures=[1:5];
        modality='ecog';
    elseif Patient==19
        Patient_initials='112';
        seizures=[1:5];
        modality='seeg';
    elseif Patient==20
        Patient_initials='114';
        seizures=[1:4];
        modality='ecog';
    elseif Patient==21
        Patient_initials='116';
        seizures=[1:3];
        modality='seeg';
    elseif Patient==22
        Patient_initials='117';
        seizures=[1:3];
        modality='seeg';
    elseif Patient==23
        Patient_initials='123';
        seizures=[1:4];
        modality='ecog';
    elseif Patient==24
        Patient_initials='126';
        seizures=[1:4];
        modality='ecog';
    elseif Patient==25
        Patient_initials='130';
        seizures=[1:5];
        modality='seeg';
    elseif Patient==26
        Patient_initials='133';
        seizures=[1:5];
        modality='seeg';
    elseif Patient==27
        Patient_initials='134';
        seizures=[1];
        modality='seeg';
    elseif Patient==28
        Patient_initials='135';
        seizures=[1:2];
        modality='seeg';
    elseif Patient==29
        Patient_initials='138';
        seizures=[1:5];
        modality='seeg';
    elseif Patient==30
        Patient_initials='139';
        seizures=[1:3];
        modality='seeg';
    elseif Patient==31
        Patient_initials='140';
        seizures=[1:3];
        modality='seeg';
    elseif Patient==32
        Patient_initials='141';
        seizures=[1:5];
        modality='seeg';
    elseif Patient==33
        Patient_initials='142';
        seizures=[1:3];
        modality='seeg';
    elseif Patient==34
        Patient_initials='144';
        seizures=[1:5];
        modality='seeg';
    elseif Patient==35
        Patient_initials='146';
        seizures=[1:3];
        modality='seeg';
    elseif Patient==36
        Patient_initials='148';
        seizures=[1:5];
        modality='seeg';
    elseif Patient==37
        Patient_initials='150';
        seizures=[1:5];
        modality='seeg';
    elseif Patient==38
        Patient_initials='151';
        seizures=[1:5];
        modality='seeg';
    elseif Patient==39
        Patient_initials='157';
        seizures=[1:5];
        modality='seeg';
    elseif Patient==40
        Patient_initials='158';
        seizures=[4];
        modality='seeg';
    elseif Patient==41
        Patient_initials='160';
        seizures=[1:3];
        modality='seeg';
    elseif Patient==42
        Patient_initials='162';
        seizures=[1:5];
        modality='seeg';
    elseif Patient==43
        Patient_initials='163';
        seizures=[1:3];
        modality='seeg';
    elseif Patient==44
        Patient_initials='164';
        seizures=[1:3];
        modality='seeg';
        %     elseif Patient==45
        %         Patient_initials='165';
        %         seizures=[1:2];
        %         modality='seeg';
    elseif Patient==45
        Patient_initials='166';
        seizures=[1:2];
        modality='seeg';
    elseif Patient==46
        Patient_initials='171';
        seizures=[1:4];
        modality='seeg';
    elseif Patient==47
        Patient_initials='172';
        seizures=[1:5];
        modality='seeg';
    elseif Patient==48
        Patient_initials='173';
        seizures=[1:5];
        modality='seeg';
    elseif Patient==49
        Patient_initials='177';
        seizures=[1:3];
        modality='seeg';
    elseif Patient==50
        Patient_initials='179';
        seizures=[1:2];
        modality='seeg';
    elseif Patient==51
        Patient_initials='180';
        seizures=[1:5];
        modality='seeg';
    elseif Patient==52
        Patient_initials='181';
        seizures=[1:5];
        modality='seeg';
    elseif Patient==53
        Patient_initials='185';
        seizures=[1:5];
        modality='seeg';
    elseif Patient==54
        Patient_initials='187';
        seizures=[1:5];
        modality='seeg';
    elseif Patient==55
        Patient_initials='188';
        seizures=[1:5];
        modality='seeg';
    elseif Patient==56
        Patient_initials='190';
        seizures=[1:3];
        modality='seeg';
    end
    
    for seizure=seizures
        try
            load([Patient_initials,'_SEEG_Seizure',num2str(seizure),'.mat'],'signal','header','channels_id')
            [Patient seizure]
        catch
            display([num2str(Patient),' ',num2str(seizure),' bad'])
        end
    end
end