%% Script to move images from session to standard MNI space 

% Register masks from session to individual space
studyDir = '/mnt/sml_share/IDENT';
bidsDir = [studyDir '/derivatives/fpp'];
subjects = {'ident01','ident02','ident03','ident04','ident05','ident06',...
    'ident07','ident08','ident09','ident10'};

roiDescs = {'Hippocampus'};

for s=1:length(subjects)
    subject = subjects{s};
    subjDir = [bidsDir '/sub-' subjects{s}];
    anatDir = [subjDir '/anat'];
    roiDir = [subjDir '/roi'];
    xfmLinPath = [anatDir '/sub-' subject '_from-individual_to-session_mode-image_xfm.mat'];
    refPath = [anatDir '/sub-' subject '_space-session_desc-LHippocampus_aparcaseg_mask.nii.gz'];
    for r=1:length(roiDescs)
        % for h=1:2
            roiPath = [anatDir '/sub-' subject '_space-individual_res-2_desc-' roiDescs{r} '_mask.nii.gz'];
            outputPath = fpp.bids.changeName(roiPath,{'space','res'},{'session',[]});
            %if exist(outputPath,'file'), continue; end
            fpp.fsl.moveImage(roiPath,refPath,outputPath,xfmLinPath,'interp','nn');
            disp(['Wrote ' subject ' ' roiDescs{r}]);
        % end 
    end
end

% Register masks from individual to MNI space
hemis = {'L','R'};
for s=1:length(subjects)
    subject = subjects{s};
    subjDir = [bidsDir '/sub-' subjects{s}];
    anatDir = [subjDir '/anat'];
    roiDir = [subjDir '/roi'];
    %xfmLinPath = [anatDir '/sub-' subjects '_from-individual_to-MNI152NLin6Asym_mode-image_xfm.mat'];
    xfmPath = [anatDir '/sub-' subject '_from-individual_to-MNI152NLin6Asym_mode-image_xfm.nii.gz'];
    refPath = [anatDir '/sub-' subject '_space-individual_res-2_desc-Hippocampus_mask.nii.gz'];
    for r=1:length(roiDescs)
        for h=1:2
            roiPath = [roiDir '/sub-' subject '_space-individual_res-2_desc-' hemis{h} roiDescs{r} '_mask.nii.gz'];
            outputPath = fpp.bids.changeName(roiPath,{'space'},{'MNI152NLin6Asym'});
            %if exist(outputPath,'file'), continue; end
            fpp.fsl.moveImage(roiPath,refPath,outputPath,[],'warp',xfmPath,'interp','nn');
            disp(['Wrote ' subject ' ' roiDescs{r}]);
        end 
    end
end
