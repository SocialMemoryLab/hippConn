%% Hippocampus subregion analysis
% Analyzes ROI distribution across head/body/tail manual segmentation and
% automatically generated subfield parcellation

%% Move to individualMTL space

studyDir = '/mnt/sml_share/IDENT';
bidsDir = [studyDir '/derivatives/fpp'];
subjects = {'ident01','ident02','ident03','ident04','ident05','ident06',...
    'ident07','ident08','ident09','ident10'};

templateSpace = 'individual_res-2';
templateSpaceNew = 'individualMTL_res-2';
regMatPath = '/mnt/sml_share/HCP/S1200AtlasBIDS32k/from-individual_to-individualMTL_mode-image_xfm.mat';
inputPathMTL = '/mnt/sml_share/IDENT/derivatives/fpp/group/langHippoc/sub-hcp100206_space-individualMTL_res-p7_desc-preproc_T1w.nii.gz';
hemis = {'L','R'};
roiDescs = {'CA1','CA23DG','Sub'};

for s=1:length(subjects)
    subject = subjects{s};
    anatDir = [bidsDir '/sub-' subject '/anat'];

    % Read in T1w
    % anatInputPath = [anatDir '/sub-' subject '_space-' templateSpace '_desc-preproc_T1w.nii.gz'];
    % anatOutputPath = [anatDir '/sub-' subject '_space-' templateSpaceNew '_desc-preproc_T1w.nii.gz'];
    % 
    % % Register individualMTL to individual space
    % fpp.fsl.moveImage(anatInputPath,inputPathMTL,anatOutputPath,regMatPath,'interp','trilinear');

    % Read in aparc segmentation
    for h=1:2
        for r=1:length(roiDescs)
            aparcInputPath = [anatDir '/sub-' subject '_space-' templateSpace '_desc-ashs' ...
                hemis{h} roiDescs{r} '_mask.nii.gz'];
            aparcOutputPath = [anatDir '/sub-' subject '_space-' templateSpaceNew '_desc-' ...
                hemis{h} roiDescs{r} '_mask.nii.gz'];
            fpp.fsl.moveImage(aparcInputPath,inputPathMTL,aparcOutputPath,regMatPath,'interp','nn');
        end 
    end 
    disp(['Finished ' subject]);
end 

%% Fix manual head/body/tail segmentations

studyDir = '/mnt/sml_share/IDENT';
bidsDir = [studyDir '/derivatives/fpp'];
subjects = {'ident01','ident02','ident03','ident04','ident05','ident06',...
    'ident07','ident08','ident09','ident10'};
hemis = {'L','R'};
segments = {'Head','Body','Tail'};
spaceName = 'individualMTL_res-2';

for s=1:length(subjects)
    subject = subjects{s};
    anatDir = [bidsDir '/sub-' subject '/anat'];

    % Read in aparc segmentation
    for h=1:2
        aparcPath = [anatDir '/sub-' subject '_space-' spaceName '_desc-' hemis{h} ...
            'Hippocampus_mask.nii.gz'];
        aparcFile = fpp.util.readDataMatrix(aparcPath);

        % Read in manual head/body/tail segmentations
        for seg=1:length(segments)
            segment = segments{seg};
            manualSegPath = [anatDir '/sub-' subject '_space-' spaceName '_desc-' hemis{h} ...
                'Hippocampus' segment '_mask.nii.gz'];
            [manualSegFile, hdr] = fpp.util.readDataMatrix(manualSegPath);

            % Remove the indices outside of aparcFile
            cleanedSeg = manualSegFile; 
            cleanedSeg(aparcFile == 0) = 0; 

            % Save corrected segmentation 
            outputPath = [anatDir '/sub-' subject '_space-' spaceName '_desc-' hemis{h} ... 
                'Hippocampus' segment '_mask.nii.gz']; 
            fpp.util.writeDataMatrix(cleanedSeg, hdr, outputPath);
            disp(['Finished ' subject ', ' segment])
        end 
    end

    % Convert manual segmentations back to individual space
    xfmLinPath = '/mnt/sml_share/HCP/S1200AtlasBIDS32k/from-individualMTL_to-individual_mode-image_xfm.mat';
    refPath = [anatDir '/sub-' subject '_space-individual_res-2_desc-LHippocampus_mask.nii.gz'];
    for seg=1:length(segments)
        segment = segments{seg};
        for h=1:2
            anatPath = [anatDir '/sub-' subject '_space-individualMTL_res-2_desc-' hemis{h} ...
                'Hippocampus' segment '_mask.nii.gz'];
            outputPath = fpp.bids.changeName(anatPath,'space','individual');
            % if exist(outputPath,'file'), continue; end
            fpp.fsl.moveImage(anatPath,refPath,outputPath,xfmLinPath,'interp','nn');
            disp(['Wrote ' subject ', ' segment]);
        end 
    end
    
    % Convert manual segmentations to session space 
    xfmLinPath = [anatDir '/sub-' subject '_from-individual_to-session_mode-image_xfm.mat'];
    refPath = [anatDir '/sub-' subject '_space-session_desc-LHippocampus_aparcaseg_mask.nii.gz'];
    for seg=1:length(segments)
        segment = segments{seg};
        for h=1:2
            anatPath = [anatDir '/sub-' subject '_space-individual_res-2_desc-' hemis{h} ...
                'Hippocampus' segment '_mask.nii.gz'];
            outputPath = fpp.bids.changeName(anatPath,'space','session');
            % if exist(outputPath,'file'), continue; end
            fpp.fsl.moveImage(anatPath,refPath,outputPath,xfmLinPath,'interp','nn');
            disp(['Wrote ' subject ', ' segment]);
        end 
    end
end

%% Analyze spatial overlap of different masks within and between individuals

% Read in langloc and fam semantic social masks in MNI space for each
% subject
studyDir = '/mnt/sml_share/IDENT';
bidsDir = [studyDir '/derivatives/fpp'];
subjects = {'ident01','ident02','ident03','ident04','ident05','ident06',...
    'ident07','ident08','ident09','ident10'};
roiDescs = {'LHippocampuslanglocSentencesVsNonwordsTop10Pct', 'RHippocampuslanglocSentencesVsNonwordsTop10Pct' ...
    'LHippocampusfamsemanticPersonVsPlaceTop10Pct', 'RHippocampusfamsemanticPersonVsPlaceTop10Pct'...
    'LHippocampusfamsemanticPersonVsPlaceInvertedTop10Pct', 'RHippocampusfamsemanticPersonVsPlaceInvertedTop10Pct'};
spaceStr = 'session';

% Initialize matrix of Dice coef values
diceMatrix = zeros(length(subjects), 6);
voxelCounts = zeros(length(roiDescs), length(subjects));
for s = 1:length(subjects)
    for r = 1:length(roiDescs)
        subject = subjects{s};
        subjDir = [bidsDir '/sub-' subject];
        anatDir = [subjDir '/anat'];
        roiDir = [subjDir '/roi'];
        roiPaths{r} = [roiDir '/sub-' subject '_space-' spaceStr '_desc-' roiDescs{r} '_mask.nii.gz'];
        try
            roiBWPaths{r} = logical(fpp.util.readDataMatrix(roiPaths{r})); 
        catch
            roiBWPaths{r} = logical(niftiread(roiPaths{r}));
        end 

        % Count voxels for each task
        voxelCount = length(find(roiBWPaths{r} == 1));
        voxelCounts(r,s) = voxelCount; 
    end 
    
    % Compute the Dice coefficient  
     % LH: langloc vs famsemantic (index 1 vs 3), langloc vs inverted (1 vs 5)
    diceMatrix(s, 1) = dice(roiBWPaths{1}, roiBWPaths{3});
    diceMatrix(s, 2) = dice(roiBWPaths{1}, roiBWPaths{5});
    diceMatrix(s, 3) = dice(roiBWPaths{3}, roiBWPaths{5});

    % RH: langloc vs famsemantic (2 vs 4), langloc vs inverted (2 vs 6)
    diceMatrix(s, 4) = dice(roiBWPaths{2}, roiBWPaths{4});
    diceMatrix(s, 5) = dice(roiBWPaths{2}, roiBWPaths{6});
    diceMatrix(s, 6) = dice(roiBWPaths{4}, roiBWPaths{6});
end 

diceAvg = zeros(1,6);
for c = 1:6
    diceAvg(c) = mean(diceMatrix(:,c));
end 
disp(diceAvg);

% Basically, there is little overlap between these regions

% Display averages
roiAvgs = zeros(length(roiDescs),1);
roiSTDs = zeros(length(roiDescs),1);
for r = 1:length(roiDescs)
    roiAvgs(r) = mean(voxelCounts(r,:));
    roiSTDs(r) = std(voxelCounts(r,:))/sqrt(length(subjects));  
end  

%% Determine the spatial distribution of lang, social, and spatial responses 

studyDir = '/mnt/sml_share/IDENT';
fppDir = [studyDir '/derivatives/fpp'];
subjects = {'ident01','ident02','ident03','ident04','ident05','ident06',...
    'ident07','ident08','ident09','ident10'};
hemis    = {'L','R'};
segments = {'Head','Body','Tail'};
roiDescs = {'langlocSentencesVsNonwordsTop5Pct', ...
            'famsemanticPersonVsPlaceTop5Pct', ...
            'famsemanticPersonVsPlaceInvertedTop5Pct'};
spaceName = 'session';

nSubs = length(subjects);
nROI  = length(roiDescs);
nSeg  = length(segments);
nHem  = length(hemis);

% Subject x ROI x Hemisphere x Segment
ROICounts = zeros(nSubs, nROI, nHem, nSeg);
segTotalCounts = zeros(nSubs, nHem, nSeg);

for s = 1:nSubs
    subject    = subjects{s};
    subjectDir = [fppDir '/sub-' subject];
    roiDir     = [subjectDir '/roi'];
    anatDir    = [subjectDir '/anat'];

    for h = 1:nHem
        for seg = 1:nSeg
            segment = segments{seg};
    
            % Read in segmentation file
            segPath = [anatDir '/sub-' subject '_space-' spaceName '_res-2_desc-' hemis{h} 'Hippocampus' segment '_mask.nii.gz'];
            segFile = fpp.util.readDataMatrix(segPath);
            segTotalVoxels = nnz(segFile);
            disp([subject ' ' hemis{h} ' ' segment ': ' num2str(segTotalVoxels) ' voxels'])
        
            for r = 1:nROI
                % Read in ROI file
                roiPath = [roiDir '/sub-' subject '_space-' spaceName '_desc-' hemis{h} 'Hippocampus' roiDescs{r} '_mask.nii.gz'];
                roiFile = fpp.util.readDataMatrix(roiPath);

                % Count overlapping labeled voxels (nonzero in both)
                overlapMask = (roiFile ~= 0) & (segFile ~= 0);
                segROIVoxels = nnz(overlapMask);

                segTotalCounts(s, h, seg) = segTotalVoxels;
                ROICounts(s, r, h, seg) = segROIVoxels;
            end
        end 
    end
    disp(['Finished ' subject]);
end

% Create index grids
[SubIdx, ROIIdx, HemIdx, SegIdx] = ndgrid(1:nSubs, 1:nROI, 1:nHem, 1:nSeg);

T = table;

T.Subject    = subjects(SubIdx(:))';
T.ROI        = roiDescs(ROIIdx(:))';
T.Hemisphere = hemis(HemIdx(:))';
T.Segment    = segments(SegIdx(:))';
T.VoxelCount = ROICounts(:);

% Stats

roiNames = {'Language','Social','Spatial'}; 

% Preallocate containers
SubjectCol   = {};
ROICol       = {};
HemCol       = {};
SegCol       = {};
PctCol       = [];

for s = 1:nSubs
    for h = 1:nHem
        for seg = 1:nSeg
            % total voxels in this ROI 
            totalROI = segTotalCounts(s,h,seg);

            for r = 1:nROI
                count = ROICounts(s,r,h,seg);
                pct   = 100 * count / totalROI;
    
                SubjectCol{end+1,1} = subjects{s};
                ROICol{end+1,1}     = roiNames{r};
                HemCol{end+1,1}     = hemis{h};
                SegCol{end+1,1}     = segments{seg};
                PctCol(end+1,1)     = pct;
            end
        end
    end
end

T = table(SubjectCol, ROICol, HemCol, SegCol, PctCol, ...
          'VariableNames', {'Subject','ROI','Hemisphere','Segment','Percent'});
T.Subject    = categorical(T.Subject);
T.ROI        = categorical(T.ROI);
T.Hemisphere = categorical(T.Hemisphere);
T.Segment    = categorical(T.Segment);

roiList = categories(T.ROI);

lmeModels = struct;

for i = 1:numel(roiList)
    thisROI = roiList{i};
    Troi = T(T.ROI == thisROI, :);

    % Mixed model
    formula = 'Percent ~ Hemisphere*Segment + (1|Subject)';
    lme = fitlme(Troi, formula);

    fprintf('\n===== ROI: %s =====\n', thisROI);
    lmeModels.([thisROI '_anova']) = anova(lme);
    disp(anova(lme))          % F-tests for Hemisphere, Segment, interaction
    disp(lme.Coefficients)    % fixed-effect estimates

    lmeModels.(char(thisROI)) = lme;
end

for i = 1:numel(roiList)
    roi = roiList{i};
    lme = lmeModels.(roi);

    fprintf('\n====================\nROI: %s\n====================\n', roi)
   
    p_results = lme_pairwise_pvals(lme, 'Segment');
    lmeModels.(['p_' roi]) = p_results;
    disp(p_results)
end

lmeOutputPath = [fppDir '/group/langHippoc/ROIHippocSegmentDistribution_LMEOutput.mat'];
%save(lmeOutputPath,"lmeModels")

%% Plot   
% --- Define ROI names ---
lmeModels = load(lmeOutputPath);
lmeModels = lmeModels.lmeModels;

roiGroups = {
    1, 'Language';
    2, 'Social';
    3, 'Spatial'
};

% --- Compute total ROI voxels per subject and ROI ---
totalPerROI = squeeze(sum(ROICounts,2));  
% size: Subjects × ROI

f = figure;
set(gcf,'Position',[200 200 1200 400])

for g = 1:size(roiGroups, 1)
    subplot(1,3,g)

    roiIdx  = roiGroups{g,1};
    roiName = roiGroups{g,2};

    % Extract voxel counts for this ROI
    % data = squeeze(ROICounts(:, roiIdx, :, :));  % size: Subjects × Hemisphere × Segment

    % Convert to percent of total ROI voxels
    pctData = zeros(nSubs,nHem,nSeg); 
    for s = 1:nSubs 
        for h = 1:nHem 
            for seg = 1:nSeg 
                rowMask = T.Subject == subjects{s} & ... 
                    T.ROI == roiName & ... 
                    T.Hemisphere == hemis{h} & ... 
                    T.Segment == segments{seg}; % Extract percent value 

                pctData(s,h,seg) = T.Percent(rowMask);
            end
        end 
    end 
    
    % Reshape
    pctData2D = reshape(pctData, nSubs, []);

    gapInd = [2 4];  
     
    % Set different plot colors 
    if isequal(g,1)
        colorArray = repmat({[1 1 .3],[1 1 .6]}, 1, 3);
    elseif isequal(g,2)
        colorArray = repmat({[217 40 34]/255,[250 95 65]/255}, 1, 3);
    else 
        colorArray = repmat({[43 119 181]/255,[96 164 208]/255}, 1, 3);
    end 
    
    % Plot
    [b,~,~,barInd] = fpp.util.barColor(pctData2D,colorArray,1,gapInd,0);
    hold on;

    % Add significance stars
    starGroups = {};
    starP = {};
    
    sigThresh = 0.05 / 3;

    pvals = lmeModels.(['p_' roiName]);

    % Head vs Body
    if pvals{1,"p"} < sigThresh
        starGroups{end+1} = [1.5 4.5];
        starP{end+1} = pvals{1,"p"};
    end

    % Body vs Tail
    if pvals{2,"p"} < sigThresh
        starGroups{end+1} = [4.5 7.5];
        starP{end+1} = pvals{2,"p"};
    end

    % Head vs Tail
    if pvals{3,"p"} < sigThresh
        starGroups{end+1} = [1.5 7.5];
        starP{end+1} = pvals{3,"p"};
    end

    % Plot stars
    starP = cell2mat(starP);
    hippConnSigStarSmallPlot(starGroups, starP);

    tickSpacing = 2;

    maxVal = max(mean(pctData2D(:,1:6))) + 1;
    yMax = ceil(maxVal/tickSpacing)*tickSpacing;
    
    ylim([0 yMax])
    yticks(0:tickSpacing:yMax)
    fontsize(gca,15,"points")
end

outputPath = [fppDir '/group/langHippoc/ROIHippocSegmentDistributionPlot.png'];
%exportgraphics(f,outputPath);
