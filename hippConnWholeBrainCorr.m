%% Compute whole-brain resting state correlations between hippocampus ROIs and cortex
% Start in session/volumetric space, then convert to surface 

studyDir = '/mnt/sml_share/IDENT';
bidsDir = [studyDir '/derivatives/fpp'];
subjects = {'ident01','ident02','ident03','ident04','ident05','ident06',...
    'ident07','ident08','ident09','ident10'};
roiDescs = {'HippocampuslanglocSentencesVsNonwords', ...
    'HippocampusfamsemanticPersonVsPlace', ...
    'HippocampusfamsemanticPersonVsPlaceInverted'};
surfaces = {'midthickness','white','pial'};
spaceStr = 'session';
inputDesc = 'preprocGSR';
hemis = {'L','R'};

for s=1:length(subjects)
    for r = 1:length(roiDescs)
        for h=1:2
            subject = subjects{s};
            subjDir = [bidsDir '/sub-' subject];
            roiDir = [subjDir '/roi'];
        
            seedPath = [roiDir '/sub-' subject '_space-' spaceStr '_desc-' hemis{h} roiDescs{r} 'Top5Pct_mask.nii.gz'];
            identRestSeedCorrelation(seedPath); 
        end 
    end 
end 

spaceStrIndiv = 'individual_res-2_den-32k';

% Z transform before averaging
for s=1 %:length(subjects)
    for r = 1:length(roiDescs)
        for h=1:2
            subject = subjects{s};
            subjDir = [bidsDir '/sub-' subject];
            anatDir = [subjDir '/anat'];
            analysisDir = [subjDir '/analysis'];
            outputDirSession = [analysisDir '/sub-' subject '_task-rest_space-' spaceStr '_funcconn'];
            outputDirIndiv = [analysisDir '/sub-' subject '_task-rest_space-' spaceStrIndiv '_funcconn'];
    
            seedPath = [roiDir '/sub-' subject '_space-' spaceStr '_desc-' hemis{h} roiDescs{r} 'Top5Pct_mask.nii.gz'];
            seedDesc = fpp.bids.checkNameValue(seedPath,'desc');
    
            inputPath = [outputDirSession '/sub-' subject '_task-rest_space-' spaceStr '_desc-' inputDesc...
               '-Seed_' seedDesc '_rstat.nii.gz'];
    
            z_outputPath = [outputDirIndiv '/sub-' subject '_task-rest_space-' spaceStrIndiv '_desc-' inputDesc...
               'ZtransformSeed' seedDesc '_rstat.nii.gz'];
            z_outputPathIndiv{h} = fpp.bids.changeName(z_outputPath,{'space','res'},{'individual','2'});    
            z_outputPathSurf{s,h,r} = fpp.bids.changeName(z_outputPathIndiv{h},{'space','res','den'},{'fsLR',[],'32k'},'rstat.dscalar','nii'); 

            % Fisher z-transform
            [inputFile, hdr] = fpp.util.readDataMatrix(inputPath);
            z_inputFile = atanh(inputFile);
            fpp.util.writeDataMatrix(z_inputFile,hdr,z_outputPath);

            % Convert to individual space 
            xfmLinPath = [anatDir '/sub-' subject '_from-session_to-individual_mode-image_xfm.mat'];
            refPath = [anatDir '/sub-' subject '_space-individual_res-2_desc-Hippocampus_mask.nii.gz'];
            if exist(z_outputPathIndiv{h},'file'), continue; end
            fpp.fsl.moveImage(z_outputPath,refPath,z_outputPathIndiv{h},xfmLinPath,'interp','trilinear');
        end 

        surfaceBasePath = [anatDir '/sub-' subject '_hemi-R_space-individual_den-native_midthickness.surf.gii'];  

        % Define surface paths
        for h=1:2
            for su=1:3
                inputSurfacePaths{h}{su} = fpp.bids.changeName(surfaceBasePath,'hemi',hemis{h},surfaces{su}); 
            end
            surfaceFsLRPaths{h} = [anatDir '/hemi-' hemis{h}...                
                '_space-fsLR_den-32k_desc-cortexAtlas_mask.shape.gii'];          
            sphereRegFsLRPaths{h} = [anatDir '/sub-' subject '_hemi-' hemis{h}, ...
                '_space-individual_den-native_desc-reg2fsLR_sphere.surf.gii'];
            midthickFsLRPaths{h} = [anatDir '/sub-' subject '_hemi-' hemis{h}, ...
                '_space-individual_den-32k_midthickness.surf.gii'];
        end 

        inputT1Path = [anatDir '/sub-' subject '_space-individual_res-p7_desc-preproc_T1w.nii.gz'];

        for h=1:2
            fpp.func.surfaceResample(z_outputPathIndiv{h},inputSurfacePaths,surfaceFsLRPaths,...
                z_outputPathSurf{s,h,r},'sphereRegFsLRPaths',sphereRegFsLRPaths,'midthickFsLRPaths',midthickFsLRPaths); 
            disp(['Finished ' subject ' ' hemis{h} roiDescs{r}]);
        end 
    end 
end 

%% Average across ROIs to produce average whole-brain functional connectivity map
% TO DO: Make this into a function 

spaceStrSurf = 'fsLR_den-32k';

for h=1:2
    for r=1:length(roiDescs) 
        files = {};
        for s=1:length(subjects)
            subject = subjects{s};
            subjDir = [bidsDir '/sub-' subject];
            analysisDir = [subjDir '/analysis'];
            outputDir = [analysisDir '/sub-' subject '_task-rest_space-' spaceStrIndiv '_funcconn'];
    
            files(cellfun(@isempty,files)) = [];
            file = [outputDir '/sub-' subject '_task-rest_space-' spaceStrSurf '_desc-' ...
                inputDesc 'ZtransformSeed' hemis{h} roiDescs{r} 'Top5Pct_rstat.dscalar.nii'];
            files{end + 1} = file;
    
            %files = cellfun(@char, files, 'UniformOutput', false); 
        end
    
        disp("Getting all data for " + hemis{h} + roiDescs{r});
        try
            allData = cellfun(@fpp.util.readDataMatrix, files, 'UniformOutput', false);
        catch
            disp(['Problem using readDataMatrix for subject ', subject])
        end 
    
        disp("Converting into 3D array for " + hemis{h} + roiDescs{r});
        data3D = cat(3, allData{:});
    
        disp("Calculating mean for " + hemis{h} + roiDescs{r})
        meandata = mean(data3D,3);
        [temp,hdr] = fpp.util.readDataMatrix(files{1});
    
        disp("Writing output for " + hemis{h} + roiDescs{r});
        outputPath = [bidsDir '/group/langHippoc/task-rest_space-' spaceStrSurf '_desc-' ...
            inputDesc 'Seed' hemis{h} roiDescs{r} 'Top5Pct_avgrstat.dscalar.nii'];
        
        % Remove z-transform
        meandata = tanh(meandata);
        fpp.util.writeDataMatrix(meandata,hdr,outputPath); 
    end
end