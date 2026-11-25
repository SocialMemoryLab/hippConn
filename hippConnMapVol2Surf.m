% Function to map volumetric ROIs in session space onto individual surface
% reconstructions

% Arguments:
% - roiDescs (cell array): BIDS descriptions of functional ROIs

function identMapVol2Surf(roiDescs)

studyDir = '/mnt/sml_share/IDENT';
bidsDir = [studyDir '/derivatives/fpp'];
subjects = {'ident01','ident02','ident03','ident04','ident05','ident06',...
    'ident07','ident08','ident09','ident10'};
hemis = {'L','R'};
surfaces = {'midthickness','white','pial'}; 
spaceStrVol = 'session';
spaceStrSurf = 'fsLR_res-2_den-32k';

for sub = 3:length(subjects)
    subject = subjects{sub};
    subjDir = [bidsDir '/sub-' subject];
    anatDir = [subjDir '/anat'];
    roiDir = [subjDir '/roi'];

    for r = 1:length(roiDescs)
        for h=1:2
            inputNiftiPath = [roiDir '/sub-' subject '_space-' spaceStrVol '_desc-' roiDescs{r} '_mask.nii.gz'];

            outputNiftiPath{h} = fpp.bids.changeName(inputNiftiPath,{'space','res'},{'individual','2'});
            outputCiftiPath{h} = [roiDir '/sub-' subject '_space-' spaceStrSurf '_desc-' roiDescs{r} '_mask.dscalar.nii'];

            % First convert to individual space 
            xfmLinPath = [anatDir '/sub-' subject '_from-session_to-individual_mode-image_xfm.mat'];
            refPath = [anatDir '/sub-' subject '_space-individual_res-2_desc-Hippocampus_mask.nii.gz'];
            % if exist(outputCiftiPath{h},'file')
            %     continue;
            % end
            % fpp.fsl.moveImage(inputNiftiPath,refPath,outputNiftiPath{h},xfmLinPath,'interp','nn');

            surfaceBasePath = [anatDir '/sub-' subject '_hemi-R_space-individual_den-native_midthickness.surf.gii'];  

            % Define surface paths
            for su=1:3
                inputSurfacePaths{h}{su} = fpp.bids.changeName(surfaceBasePath,'hemi',hemis{h},surfaces{su}); 
            end
            surfaceFsLRPaths{h} = [anatDir '/hemi-' hemis{h}...                
                '_space-fsLR_den-32k_desc-cortexAtlas_mask.shape.gii'];          
            sphereRegFsLRPaths{h} = [anatDir '/sub-' subject '_hemi-' hemis{h}...
                '_space-individual_den-native_desc-reg2fsLR_sphere.surf.gii'];
            midthickFsLRPaths{h} = [anatDir '/sub-' subject '_hemi-' hemis{h}...
                '_space-individual_den-32k_midthickness.surf.gii'];
        end 

        for h = 1:2
            fpp.func.surfaceResample(outputNiftiPath{h},inputSurfacePaths,surfaceFsLRPaths,outputCiftiPath{h}, ...
                'sphereRegFsLRPaths',sphereRegFsLRPaths,'midthickFsLRPaths',midthickFsLRPaths);
            fpp.wb.command('cifti-math', ' x > 0 ', [], outputCiftiPath{h}, [' -var x ' outputCiftiPath{h}]);
        end 
        disp(['Finished ' subject ' ' roiDescs{r}]);
    end     
end 

%% Correlation rstat files
studyDir = '/mnt/sml_share/IDENT';
bidsDir = [studyDir '/derivatives/fpp'];
subjects = {'ident01','ident02','ident03','ident04','ident05','ident06',...
    'ident07','ident08','ident09','ident10'};
roiDescs = {'HippocampuslanglocSentencesVsNonwordsTop5Pct', ...
    'HippocampusfamsemanticPersonVsPlaceTop5Pct', ...
    'HippocampusfamsemanticPersonVsPlaceInvertedTop5Pct'};
hemis = {'L','R'};
surfaces = {'midthickness','white','pial'}; 
spaceStrVol = 'session';
spaceStrSurf = 'individual_res-2_den-32k';

for s=1:length(subjects)
    for r = 1:length(roiDescs)
        for h = 1:2
            subject = subjects{s};
            subjDir = [bidsDir '/sub-' subject];
            anatDir = [subjDir '/anat'];
            roiDir = [subjDir '/roi'];
            analysisDir = [subjDir '/analysis'];

            inputNiftiPath = [analysisDir '/sub-' subject '_task-rest_space-session_funcconn/sub-' subject ...
                '_task-rest_space-session_desc-preproc-Seed_' hemis{h} roiDescs{r} '_rstat.nii.gz'];
            outputNiftiPath{h} = fpp.bids.changeName(inputNiftiPath,{'space','res'},{'individual','2'});
            outputCiftiPath{h} = [analysisDir '/sub-' subject '_task-rest_space-' spaceStrSurf '_funcconn/sub-' subject ...
                '_task-rest_space-' spaceStrSurf '_desc-preproc-Seed_' hemis{h} roiDescs{r} '_rstat.dscalar.nii'];

            % First convert to individual space 
            xfmLinPath = [anatDir '/sub-' subject '_from-session_to-individual_mode-image_xfm.mat'];
            refPath = [anatDir '/sub-' subject '_space-individual_res-2_desc-Hippocampus_mask.nii.gz'];
            if exist(outputNiftiPath,'file'), continue; end
            fpp.fsl.moveImage(inputNiftiPath,refPath,outputNiftiPath,xfmLinPath,'interp','nn');

            surfaceBasePath = [anatDir '/sub-' subject '_hemi-R_space-individual_den-native_midthickness.surf.gii'];  

            % Define surface paths
            for su=1:3
                inputSurfacePaths{h}{su} = fpp.bids.changeName(surfaceBasePath,'hemi',hemis{h},surfaces{su}); 
            end
            surfaceFsLRPaths{h} = [anatDir '/hemi-' hemis{h}...                
                '_space-fsLR_den-32k_desc-cortexAtlas_mask.shape.gii'];          
            sphereRegFsLRPaths{h} = [anatDir '/sub-' subject '_hemi-' hemis{h}...
                '_space-individual_den-native_desc-reg2fsLR_sphere.surf.gii'];
            midthickFsLRPaths{h} = [anatDir '/sub-' subject '_hemi-' hemis{h}...
                '_space-individual_den-32k_midthickness.surf.gii'];
        end 

        inputT1Path = [anatDir '/sub-' subject '_space-individual_res-p7_desc-preproc_T1w.nii.gz']; 

        for h = 1:2
            fpp.func.surfaceResample(outputNiftiPath{h},inputSurfacePaths,surfaceFsLRPaths,outputCiftiPath{h}, ...
                'sphereRegFsLRPaths',sphereRegFsLRPaths,'midthickFsLRPaths',midthickFsLRPaths); 
            disp(['Finished ' subject ' ' hemis{h} roiDescs{r}]);
        end 
    end 
end 