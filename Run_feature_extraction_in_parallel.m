

function Run_feature_extraction_in_parallel(poolSize,subjs)
parpool(poolSize); % Warning: you must call matlabpool with brackets
parfor subj=subjs
    [Subject] = Feature_extraction_parallel(subj);
    display(['Subj #',num2str(Subject),' Done!'])
end
delete(gcp('nocreate'))
end


