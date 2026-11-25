%% Wrapper for identROIExtract
% Creates social, spatial, and language hippocampal ROIs and additionally creates group-average ROIs

%% Extract hippocampal ROIs from aparc+aseg

studyDir = '/mnt/sml_share/IDENT';
bidsDir = [studyDir '/derivatives/fpp'];
subjects = {'ident01','ident02','ident03','ident04','ident05','ident06',...
    'ident07','ident08','ident09','ident10'};

% Session space
for s=1:length(subjects)
    subject = subjects{s};
    anatDir = [bidsDir '/sub-' subject '/anat'];

    parcPathSession = [anatDir '/sub-' subject '_space-session_desc-aparcaseg_dseg.nii.gz'];
    
    parcIds = {17, 53}; % Left, right hippocampus
    hemis = {'L','R'};
    for h=1:2
        labelName = [hemis{h} 'Hippocampus'];
        outputPath = fpp.bids.changeName(parcPath, 'desc', labelName,'mask',[]);

        fpp.util.label2ROI(parcPath, parcIds{h}, outputPath); % NOT fslmaths

        if exist(outputPath,"file")
            disp(['Generated ' labelName ' mask file']);
        else
            disp(['Could not generate file for subject ' subjects{s}]);
        end 
    end 
end

%% Find tSNR in the hippocampus

studyDir = '/mnt/sml_share/IDENT';
bidsDir = [studyDir '/derivatives/fpp'];
subjects = {'ident01','ident02','ident03','ident04','ident05','ident06',...
    'ident07','ident08','ident09','ident10'};
roiDescs = {'LHippocampus', 'RHippocampus'};
nRuns = 1:6;
spaceStr = 'session';

tsnrValues = cell(length(subjects), length(nRuns));

for s=1:length(subjects)
    for reg=1:length(roiDescs)
        for r=1:length(nRuns)
            % Load resting file
            subject = subjects{s};
            roi = roiDescs{reg};
            run = num2str(nRuns(r));
            funcDir = [bidsDir '/sub-' subject '/func'];
            anatDir = [bidsDir '/sub-' subject '/anat'];
    
            rsPath = [funcDir '/sub-' subject '_task-rest_run-0' run '_space-session_desc-preproc_tsnr.nii.gz'];
            maskPath = [anatDir '/sub-' subject '_space-' spaceStr '_desc-' roi '_aparcaseg_mask.nii.gz'];
            maskFile = niftiread(maskPath);
    
            rsFile = fpp.util.readDataMatrix(rsPath, maskFile);
    
            meanROItsnr = mean(rsFile, 'all');
            fprintf('Avg tSNR: %.4f\n', meanROItsnr);
            tsnrValues(s,r) = meanROItsnr;
        end  
    end 
end

tsnrValues = cell2mat(tsnrValues);
tsnrMean = mean(tsnrValues, 'all');
tsnrStd = std(tsnrValues, 0, 'all');

disp([tsnrMean, tsnrStd]);

%%  Extract language responses in hippocampal ROIs

searchNames = {'LHippocampus_aparcaseg','RHippocampus_aparcaseg'}; 
for s=1:length(searchNames)
    identROIExtract('langloc','SentencesVsNonwords',searchNames{s},0,'roiSize', 5);
    identROIExtract('famsemantic','PersonVsPlace',searchNames{s},0,'roiSize', 5);
    identROIExtract('famsemantic','PersonVsPlace',searchNames{s},0,'roiSize', 5,'invertStats',1);
    identROIExtract('famsemantic','PersonVsObject',searchNames{s},0,'roiSize',5,'invertStats',1);
end

%% Display p-values for each contrast

fppDir = '/mnt/sml_share/IDENT/derivatives/fpp';
searchNames = {'LHippocampus', 'RHippocampus'};
taskNames = {'langloc','famsemantic'};
contrastNames = {'SentencesVsNonwords','PersonVsPlace','PersonVsPlaceInverted'};
pct = 5;
templateSpace = 'session';

statTables = cell(length(contrastNames), length(searchNames));
for s=1:length(searchNames)
    for t=1:length(taskNames)
        for c=1:length(contrastNames)
            groupDir = [fppDir '/group/langHippoc'];
            outputDesc = [searchNames{s} taskNames{t} contrastNames{c} 'Top' num2str(pct) 'PctN10'];
            outputPath = [groupDir '/space-' templateSpace '_desc-' outputDesc '_roiData.mat'];
            if exist(outputPath,'file')
                disp([newline '====== ' searchNames{s}, taskNames{t}, contrastNames{c} ' ======' newline]);
                statTable = identROIDispStats(taskNames{t},contrastNames{c},searchNames{s},0);

                statTable.Properties.VariableNames = {[searchNames{s} '_tStat'], [searchNames{s} '_pVal']};
                statTables{c,s} = statTable;
            else
                continue; 
            end 
        end
    end 
end

tables = {
    [statTables{1,1}, statTables{1,2}];
    [statTables{2,1}, statTables{2,2}];
    [statTables{3,1}, statTables{3,2}]
};

tableNames = {'langTable','socTable','spatTable'};

for t = 1:numel(tables)
    tableOutputPath = fullfile(groupDir, tableNames{t} + ".xlsx");
    writetable(tables{t}, tableOutputPath, 'WriteRowNames', true);
end


%% Extract responses at a range of ROI sizes

searchNames = {'LHippocampus_aparcaseg'}; 
for s=1:length(searchNames)
    identROIExtractBySize('SentencesVsNonwords',searchNames{s},0,'invertStats',0)
end 

%% Average ROIs across subjects to produce probability map for each hippocampal seed

studyDir = '/mnt/sml_share/IDENT';
bidsDir = [studyDir '/derivatives/fpp'];
subjects = {'ident01','ident02','ident03','ident04','ident05','ident06',...
    'ident07','ident08','ident09','ident10'};
roiDescs = {'HippocampuslanglocSentencesVsNonwords', ...
    'HippocampusfamsemanticPersonVsPlace', ...
    'HippocampusfamsemanticPersonVsPlaceInverted'};
spaceStr = 'MNI152NLin6Asym';
hemis = {'L','R'};

for r=1:length(roiDescs)
    for h=1:2
        cmd = '';
        for s=2:length(subjects)
            subject = subjects{s};
            subjDir = [bidsDir '/sub-' subject];
            roiDir = [subjDir '/roi'];
            inputPath = [roiDir '/sub-' subject '_space-' spaceStr '_res-2_desc-' hemis{h} roiDescs{r} 'Top5Pct_mask.nii.gz'];
            cmd = [cmd ' -add ' inputPath];
        end
        cmd = [cmd ' -div ' int2str(length(subjects))];
        subject = subjects{1};
        inputPath = [bidsDir '/sub-' subject '/roi/sub-' subject '_space-' spaceStr '_res-2_desc-' hemis{h} roiDescs{r} '_mask.nii.gz'];
        outputPath = [bidsDir '/group/langHippoc/space-' spaceStr '_res-2_desc-' hemis{h} roiDescs{r} 'CrossSub_probseg.nii.gz'];
        %fpp.fsl.maths(inputPath,cmd,outputPath);
        disp(['Wrote ' hemis{h} roiDescs{r}]);
    end
end

% Upsample to p5 resolution for visualization 
spaceStr = 'MNI152NLin6Asym';
for h=1:2
    for r=1:length(roiDescs)
        inputPath_res2 = [bidsDir '/group/langHippoc/space-' spaceStr '_res-2_desc-' ...
            hemis{h} roiDescs{r} 'Top5PctCrossSub_probseg.nii.gz'];
        outputPath_resp5 = [bidsDir '/group/langHippoc/sub-ident01_space-' spaceStr '_res-p5_desc-' ...
            hemis{h} roiDescs{r} 'Top5PctCrossSub_probseg.nii.gz'];
        refImg = '/mnt/sml_share/IDENT/derivatives/fpp/group/mni_icbm152_t1_tal_nlin_asym_09a.nii';
        % xfmLinPath = [anatDir '/sub-' subject '_from-individual-to-MNI152NLin6Asym_mode-image_xfm.mat'];
        
        fpp.wb.command('volume-resample',inputPath_res2,refImg,'TRILINEAR',outputPath_resp5);
        
        disp(['Moved ' hemis{h} roiDescs{r}]);
    end 
end 

%% Determine whether fROIs are intercorrelated

% Pearson correlations, small ROIs
corrMat = load('/mnt/sml_share/IDENT/derivatives/fpp/group/langHippoc/desc-HippocSeedsNoOverlap-to-mmpfROIsTop5PctUpdatedN10_ROICorrelationData.mat');
roiDescs = corrMat.roiDescs;
corrMat = corrMat.corrMat;

meanCorrs = mean(corrMat,3);

% Group ROIs by type
%hippROIs = contains(roiDescs,'Hipp');
langROIs    = contains(roiDescs, 'Sentences') & contains(roiDescs, 'Hipp');
socialROIs  = contains(roiDescs, 'PersonVsPlaceT') & contains(roiDescs, 'Hipp'); %& ~contains(roiDescs,'THfam');
spatialROIs = contains(roiDescs, 'PersonVsPlaceInverted') & contains(roiDescs, 'Hipp'); %& ~contains(roiDescs,'TPfam');

roiGroups = {langROIs, socialROIs, spatialROIs};
groupNames = {'Language', 'Social', 'Spatial'};

% Compute within-group and between-group correlations
results = struct();

for i = 1:length(groupNames)
    maskA = roiGroups{i};
    for j = 1:length(groupNames)
        maskB = roiGroups{j};

        % Skip if group has no ROIs
        if sum(maskA)==0 || sum(maskB)==0
            continue
        end

        % Extract correlations between group i and group j
        groupCorrs = meanCorrs(maskA, maskB);

        % For within-group, only take upper triangle (no duplicates)
        % if i == j
        %     groupCorrs = groupCorrs(triu(true(size(groupCorrs)),1));
        % end

        % Store results
        results(i,j).nameA = groupNames{i};
        results(i,j).nameB = groupNames{j};
        results(i,j).corrVals = groupCorrs(:);
        results(i,j).meanCorrs = mean(groupCorrs(:));
    end
end

% Print group-wise averages
fprintf('\n=== Mean Correlations Between Groups ===\n');
for i = 1:length(groupNames)
    for j = 1:length(groupNames)
        if isempty(results(i,j).corrVals)
            continue
        end
        fprintf('%s ↔ %s: mean r = %.3f\n', ...
            results(i,j).nameA, results(i,j).nameB, results(i,j).meanCorrs);
    end
end

% Statistical comparisons between groups
fprintf('\n=== Statistical Comparisons (t-tests) ===\n');

pairsToCompare = [
    1 1 1 2;  % Lang within vs Lang-Social
    1 1 1 3;  % Lang within vs Lang-Spatial
    2 2 2 3;  % Social within vs Social-Spatial
    3 3 3 1;  % Spatial within vs Spatial-Lang
    1 2 2 3;  % Lang-Social vs Social-Spatial
    1 3 2 3;  % Lang-Spatial vs Social-Spatial
    1 2 1 3;  % Lang-Social vs Lang-Spatial
];

for p = 1:size(pairsToCompare,1)
    i1 = pairsToCompare(p,1); j1 = pairsToCompare(p,2);
    i2 = pairsToCompare(p,3); j2 = pairsToCompare(p,4);

    vals1 = results(i1,j1).corrVals;
    vals2 = results(i2,j2).corrVals;

    [~, pval, ~, stats] = ttest2(vals1, vals2);
    fprintf('%s↔%s vs %s↔%s: t(%d)=%.2f, p=%.4f\n', ...
        results(i1,j1).nameA, results(i1,j1).nameB, ...
        results(i2,j2).nameA, results(i2,j2).nameB, ...
        stats.df, stats.tstat, pval);
end

% Optional: visualize as matrix
meanMat = nan(length(groupNames));
for i = 1:length(groupNames)
    for j = 1:length(groupNames)
        if ~isempty(results(i,j).corrVals)
            meanMat(i,j) = results(i,j).meanCorrs;
        end
    end
end

figure;
imagesc(meanMat, [-1 1]);
colorbar;
axis equal tight;
set(gca, 'XTick', 1:3, 'XTickLabel', groupNames, ...
         'YTick', 1:3, 'YTickLabel', groupNames);
title('Mean ROI Correlations Between Groups (Hippocampus)'); 
