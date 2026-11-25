%% Corticohippocampal connectivity analysis
% Extracts functional ROIs from cortical social, spatial, and language
% search spaces and computes partial correlations with
% hippocampal functional ROIs using resting-state data


%% ------ Extract ROIs from cortical search spaces -------
studyDir = '/mnt/sml_share/IDENT';
bidsDir = [studyDir '/derivatives/fpp'];
subjects = {'ident01','ident02','ident03','ident04','ident05','ident06',...
    'ident07','ident08','ident09','ident10'};

% Search space names for networks
searchNames{1} = {'PostTemp','AntTemp','MFG','IFG','TPJ'}; 
searchNames{2} = {'MPFC','MPC','TPJ','ASTS','SFG','TP'};
searchNames{3} = {'MPFC','MPC','TPJ','ASTS','SFG','TH'};

mmpNames = {'mmpLng','mmpApex','mmpApex'};
netNames = {'Lang','Social','Spatial'};
contrastNames = {'SentencesVsNonwords','PersonVsPlace','PersonVsPlaceInverted'};
taskNames = {'langloc','famsemantic','famsemantic'};
roiSize = 5;
spaceStr = 'session';
hemis = {'L','R'};

% ----- Extract ROIs ------
roiCortex = {};
for n = 1:length(netNames)
    for s = 1:length(searchNames{n})
        for h = 1:2
            if n < 3
                if n == 1 && contains(searchNames{n}{s},'TPJ')
                    %identROIExtract(taskNames{n}, contrastNames{n}, [mmpNames{2} hemis{h} searchNames{n}{s}], 0, 'roiSize', roiSize);
                else 
                    %identROIExtract(taskNames{n}, contrastNames{n}, [mmpNames{n} hemis{h} searchNames{n}{s}], 0, 'roiSize', roiSize);
                end 
            elseif n == 3
                %identROIExtract(taskNames{n}, contrastNames{n}, [hemis{h} searchNames{n}{s}], 0, 'invertStats', 1, 'roiSize', roiSize);
            end 
        end 
    end 
end 

%% ------ Merge into one search space for each task -------

nNets = length(netNames);
nROIs = max(cellfun(@length, searchNames));
searchNamesComb = cell(nNets, 1);

for n = 1:nNets
    for s = 1:length(searchNames{n})
        for h = 1:2
            if n == 1 && contains(searchNames{n}{s},'TPJ')
                searchNamesComb{n}{s, h} = [ ...
                    mmpNames{2} ...
                    hemis{h} ...
                    searchNames{n}{s} ...
                    taskNames{n} ...
                    contrastNames{n}, ...
                    'Top' num2str(roiSize) 'Pct'];
            else 
                searchNamesComb{n}{s, h} = [ ...
                    mmpNames{n} ...
                    hemis{h} ...
                    searchNames{n}{s} ...
                    taskNames{n} ...
                    contrastNames{n}, ...
                    'Top' num2str(roiSize) 'Pct'];
            end 
        end
    end
end

% ------ Combine ROIs into networks ------
for s = 1:length(subjects)
    subject = subjects{s};
    subjDir = [bidsDir '/sub-' subject];
    roiDir = [subjDir '/roi'];

    roiCortex = {};
    for n = 1:length(netNames)
        roiCortex{end+1} = [mmpNames{n} netNames{n} 'Top' num2str(roiSize) 'Pct'];
        outputPath = fullfile(roiDir,['sub-' subject '_space-' spaceStr '_desc-' roiCortex{n} '_mask.nii.gz']);

        cmdParts = {};
        for r = 1:numel(searchNamesComb{n})
            roiPath = fullfile(roiDir, ...
                ['sub-' subject '_space-' spaceStr '_desc-' searchNamesComb{n}{r} '_mask.nii.gz']);
            cmdParts{end+1} = roiPath;
        end
        cmd = strjoin(cmdParts, ' -add ');
        fslCmd = sprintf('fslmaths %s -add 0 %s', cmd, outputPath);
        %system(fslCmd);

        disp(['Finished ' subject ' ' netNames{n}]);
    end
end 

identMapVol2Surf(roiCortex);   

%% Generate correlation matrices and combine them in a struct

% Hippocampal ROIs  
roiHippoc = {'HippocampuslanglocSentencesVsNonwords', ...
    'HippocampusfamsemanticPersonVsPlace', ...
    'HippocampusfamsemanticPersonVsPlaceInverted'};

corrStruct = struct();

for h=1:2
    for r=1:length(roiHippoc)
        outputDesc = [hemis{h} roiHippoc{r} 'Top' num2str(roiSize) 'Pct-to-mmpNetworksPartialCorr' ];
        outputPath = ['/mnt/sml_share/IDENT/derivatives/fpp/group/langHippoc/desc-' outputDesc '_ROICorrelationData.mat'];
       
        % Run partial correlations for each ROI separately 
        if ~exist(outputPath,'file')
            identRestROIVolSurfPartialCorr([hemis{h} roiHippoc{r} 'Top' num2str(roiSize) 'Pct'], roiCortex, outputDesc); 
            disp(['Finished ' outputDesc]);
        else
            disp(['Corr mat already written for ' roiHippoc{r}]);
        end 
    
        corrMat = load(outputPath);
        corrMat = corrMat.corrMat;
        corrStruct.([hemis{h} roiHippoc{r}]) = corrMat;
        delete(outputPath);
    end 
end

structOutputPath = '/mnt/sml_share/IDENT/derivatives/fpp/group/langHippoc/desc-HippocSeedsTop5Pct-to-mmpNetworksPartialCorr_ROICorrelationData.mat';
save(structOutputPath, 'corrStruct')


%% Run stats for ROI-to-mmp

corrStruct = load(structOutputPath);
corrStruct = corrStruct.corrStruct;

studyDir = '/mnt/sml_share/IDENT';
bidsDir = [studyDir '/derivatives/fpp'];
subjects = {'ident01','ident02','ident03','ident04','ident05','ident06',...
    'ident07','ident08','ident09','ident10'};
roiHippoc = {'HippocampuslanglocSentencesVsNonwords', ...
    'HippocampusfamsemanticPersonVsPlace', ...
    'HippocampusfamsemanticPersonVsPlaceInverted'};
hemis = {'L','R'};

netNames = {'Lang','Social','Spatial'};
rows = {};

% Build LME table
for h=1:2
    for r=1:length(roiHippoc)
        for s = 1:length(subjects)

            rVals = corrStruct.([hemis{h} roiHippoc{r}])(:,:,s); 
            
            if isequal(h,1)
                hemi = hemis{1};
            else 
                hemi = hemis{2};
            end 

            % Fisher z-transform
            FCz = atanh(rVals);
            for n = 2:4
                rows{end+1,1} = FCz(1,n);         % correlation value
                rows{end,2}   = subjects{s};     % subject
                rows{end,3}   = hemi;
                rows{end,4}   = roiHippoc{r};    % hippocampal ROI
                rows{end,5}   = netNames{n-1};   % network name
            end
        end
    end 
end

lmeTable = cell2table(rows, ...
    'VariableNames', {'y','subj','hemis','roiHippoc','network'});

% Convert to categorical
lmeTable.subj      = categorical(lmeTable.subj);
lmeTable.hemis     = categorical(lmeTable.hemis);
lmeTable.roiHippoc = categorical(lmeTable.roiHippoc);
lmeTable.network   = categorical(lmeTable.network);

% Full mixed model across all hippocampal ROIs
lme_full = fitlme(lmeTable, 'y ~ hemis*roiHippoc*network + (1|subj)');

disp('=== ANOVA RESULT for LME ===');
corrStruct.anova_full = anova(lme_full);
disp(anova(lme_full));

%% POST HOC TESTS PER ROI

coefNames = lme_full.CoefficientNames;
beta = lme_full.Coefficients.Estimate;

results = {};
row = 1;

roiLevels = categories(lmeTable.roiHippoc);
refROI = roiLevels{1};

for r = 1:length(roiHippoc)

    roi = roiHippoc{r};

    comparisons = {
        'Social'  'Lang'
        'Spatial' 'Lang'
        'Spatial' 'Social'
        };

    for c = 1:size(comparisons,1)

        net1 = comparisons{c,1};
        net2 = comparisons{c,2};

        L = zeros(1,length(coefNames));

        % Main-effect portion
        if ~strcmp(net1,'Lang')
            idx = strcmp(coefNames,['network_' net1]);
            L(idx) = L(idx) + 1;
        end

        if ~strcmp(net2,'Lang')
            idx = strcmp(coefNames,['network_' net2]);
            L(idx) = L(idx) - 1;
        end

        % ROI x Network interaction
        if ~strcmp(roi,refROI)   

            if ~strcmp(net1,'Lang')

                intName = sprintf('roiHippoc_%s:network_%s', roi,net1);
                idx = find(strcmp(coefNames,intName));

                if any(idx)
                    L(idx) = L(idx) + 1;
                end
            end

            if ~strcmp(net2,'Lang')

                intName = sprintf('roiHippoc_%s:network_%s', roi,net2);
                idx = find(strcmp(coefNames,intName));

                if any(idx)
                    L(idx) = L(idx) - 1;
                end
            end

        end

        disp(roi)
        disp(L)
        [p,F,df1,df2] = coefTest(lme_full,L);

        estimate = L*beta;
        se = sqrt(L*lme_full.CoefficientCovariance*L');
        tstat = estimate/se;

        results(row,:) = { ...
            roi,...
            net1,...
            net2,...
            estimate,...
            tstat,...
            p,...
            df2};

        row = row + 1;

    end
end

networkPostHoc = cell2table(results,...
    'VariableNames',...
    {'ROI','Network1','Network2',...
     'Estimate','t','p','DF'});

disp(networkPostHoc)

%% POST HOC TESTS PER NETWORK

coefNames = lme_full.CoefficientNames;

results = {};
row = 1;

netLevels = categories(lmeTable.network);
refNet = netLevels{1};

for n = 1:length(netLevels)

    net = netLevels{n};

    comparisons = {
        roiHippoc{2} roiHippoc{1}
        roiHippoc{3} roiHippoc{1}
        roiHippoc{3} roiHippoc{2}
        };

    for c = 1:size(comparisons,1)

        roi1 = comparisons{c,1};
        roi2 = comparisons{c,2};

        L = zeros(1,length(coefNames));

        % ROI main effects

        coef1 = ['roiHippoc_' roi1];
        coef2 = ['roiHippoc_' roi2];

        if any(strcmp(coefNames,coef1))
            L(strcmp(coefNames,coef1)) = ...
                L(strcmp(coefNames,coef1)) + 1;
        end

        if any(strcmp(coefNames,coef2))
            L(strcmp(coefNames,coef2)) = ...
                L(strcmp(coefNames,coef2)) - 1;
        end

        % ROI x Network interactions

        if ~strcmp(net,refNet)

            int1 = ['roiHippoc_' roi1 ':network_' net];

            int2 = ['roiHippoc_' roi2 ':network_' net];

            if any(strcmp(coefNames,int1))
                L(strcmp(coefNames,int1)) = ...
                    L(strcmp(coefNames,int1)) + 1;
            end

            if any(strcmp(coefNames,int2))
                L(strcmp(coefNames,int2)) = ...
                    L(strcmp(coefNames,int2)) - 1;
            end
        end

        [p,F,df1,df2] = coefTest(lme_full,L);

        beta = lme_full.Coefficients.Estimate;

        estimate = L*beta;

        se = sqrt(L *lme_full.CoefficientCovariance * L');

        tstat = estimate/se;

        results(row,:) = { ...
            net,...
            roi1,...
            roi2,...
            estimate,...
            tstat,...
            p,...
            df2};

        row = row + 1;

    end
end

roiPostHoc = cell2table(results,...
    'VariableNames',...
    {'Network','ROI1','ROI2',...
     'Estimate','t','p','DF'});

disp(roiPostHoc)

%% OLD APPROACH

for r=1:length(roiHippoc)
    fprintf('\n=== POST-HOC TESTS FOR %s ===\n', roiHippoc{r});

    pairwisePvals_roi = lme_pairwise_pvals(lme_full,'network');
    disp(pairwisePvals_roi);
    corrStruct.(['p_' roiHippoc{r}]) = pairwisePvals_roi;
end 

% POST-HOC TESTS PER NETWORK
for r=1:length(netNames)
    fprintf('\n=== POST-HOC TESTS FOR %s ===\n', netNames{r});

    pairwisePvals_net = lme_pairwise_pvals(lme_full,'roiHippoc');
    disp(pairwisePvals_net);
    corrStruct.(['p_' netNames{n}]) = pairwisePvals_net;
end 

save(structOutputPath, 'corrStruct')

%% mmp Network plot

corrMat = load(structOutputPath);
corrMat = corrMat.corrStruct;

nSubs = length(subjects);
nHipp = numel(roiHippoc);
nNet = numel(netNames);

% Fisher transform and reshape 
for s = 1:nSubs
    colIdx = 1;
    for r = 1:nHipp
        for h=1:2
            % Correlation matrix for this ROI and subject
            rMat = corrMat.([hemis{h} roiHippoc{r}])(:,:,s);
    
            % Fisher z-transform
            FCz = atanh(rMat);
    
            % Networks are in columns 2:4 (Lang, Social, Spatial)
            for n = 1:nNet
                FCz_reshaped(s, colIdx) = FCz(1, n+1);
                colIdx = colIdx + 1;
            end
        end
    end
end 

% --- Plot ---
f = figure('Position',[100 100 1000 600]);

nConds = [3 3 3 3 3 3];  
lastConds = [0 cumsum(nConds)];
gapInd = lastConds(2:end-1);  
 
colorArray = repmat({[1 1 0],[250 95 65]/255,[96 164 208]/255}, 1, 6);

fpp.util.barColor(FCz_reshaped,colorArray,1,gapInd,0);
hold on;

barX = [1 2 3; 5 6 7; 9 10 11; 13 14 15; 17 18 19; 21 22 23];

starGroups = {};
starP = {};

sigThresh = 0.05 / 3;

for r = 1:nHipp
    for h = 1:2

        idx = (r-1)*2 + h;

        pvals = corrMat.(['p_' hemis{h} roiHippoc{r}]);

        % Lang vs Social
        if pvals{1,"p"} < sigThresh
            starGroups{end+1} = [barX(idx,1) barX(idx,2)];
            starP{end+1} = pvals{1,"p"};
        end

        % Lang vs Spatial
        if pvals{2,"p"} < sigThresh
            starGroups{end+1} = [barX(idx,1) barX(idx,3)];
            starP{end+1} = pvals{2,"p"};
        end

        % Social vs Spatial
        if pvals{3,"p"} < sigThresh
            starGroups{end+1} = [barX(idx,2) barX(idx,3)];
            starP{end+1} = pvals{3,"p"};
        end
    end
end

% Plot stars
starP = cell2mat(starP);
hippConnSigStarSmallPlot(starGroups, starP);
set(gca, 'YTick', -.1:.1:3, 'FontSize', 14)

%ylabel(ax, 'Mean correlation', 'FontSize', 15);
%title(ax, 'Correlation between hippocampal and cortical fROIs', 'FontSize', 20);
hold off;

outputDesc = 'desc-HippocSeedsTop5Pct-to-mmpNetworksPartialCorr';
exportgraphics(f, ['/mnt/sml_share/IDENT/derivatives/fpp/group/langHippoc/' outputDesc '.png'], 'BackgroundColor', 'none');
