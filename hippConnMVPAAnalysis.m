%% MVPA analysis
% Uses an SVM to classify language vs nonwords runs and social vs spatial
% runs within the semantic task and across tasks


% Basic information 
subjects = {'ident01','ident02','ident03','ident04','ident05','ident06','ident07','ident08','ident09','ident10'};
roiMaskNames = {'LHippocampus_aparcaseg','RHippocampus_aparcaseg'};
templateSpace = 'session';
fppDir = '/mnt/sml_share/IDENT/derivatives/fpp';

%% -----------------------------
%  Task label definitions
%% -----------------------------
% Number of runs
nLang = 4; nSoc = 5; nSpat = 5;

% taskNamesWithin = {'FamPerson', 'FamPlace'};
taskLabelsWithin = [1 1 1 1 1 2 2 2 2 2];
runLabelsWithin = [1 2 3 4 5 1 2 3 4 5]; 

% taskNamesBetween = {'VisualFamPerson', 'VisualFamPlace', ...
%     'SemanticFamPerson','SemanticFamPlace', ...
%     'EpisodicFamPerson','EpisodicFamPlace'};
taskLabelsBetween = [repmat(1,1,nSoc), repmat(2,1,nSpat), ...
    repmat(1,1,nSoc), repmat(2,1,nSpat), ...
    repmat(1,1,nSoc), repmat(2,1,nSpat)];
runLabelsBetween = repmat(1:5,1,6); 

% taskNamesLangloc = {'Sentence', 'NonwordList'};
taskLabelsLangloc = [1 1 1 1 2 2 2 2];
runLabelsLangloc  = [1 2 3 4 1 2 3 4]; 

nSubj = numel(subjects);
nROI  = numel(roiMaskNames);

withinAcc  = nan(nSubj,nROI);
betweenAcc = nan(nSubj,nROI);
langlocAcc = nan(nSubj,nROI);
withinConf = {};
betweenConf = {};
langlocConf = {};

% Permutation test set up
nPerm = 1000;

yPermBetween = nan(nPerm,numel(taskLabelsBetween));
yPermWithin = nan(nPerm,numel(taskLabelsWithin));
yPermLangloc = nan(nPerm,numel(taskLabelsLangloc));

permAccBetween = nan(nSubj,nROI,nPerm);
permAccWithin = nan(nSubj,nROI,nPerm);
permAccLangloc = nan(nSubj,nROI,nPerm);

for p=1:nPerm
    yPermBetween(p,:) = taskLabelsBetween(randperm(numel(taskLabelsBetween)));
    yPermWithin(p,:) = taskLabelsWithin(randperm(numel(taskLabelsWithin)));
    yPermLangloc(p,:) = taskLabelsLangloc(randperm(numel(taskLabelsLangloc)));
end 

%% -----------------------------
%  Loop over subjects & ROIs
%% -----------------------------

outputPath = [fppDir '/group/langHippoc/identMVPA_SVM_Results.mat'];

for s = 1:nSubj
    subject = subjects{s};
    subjDir = fullfile(fppDir, ['sub-' subject]);
    analysisDir = fullfile(subjDir, 'analysis');
    anatDir = fullfile(subjDir, 'anat');

    % Load social/spatial betas 
    beta = load(fullfile(analysisDir, ...
        ['sub-' subject '_space-' templateSpace '_desc-VisualSemanticEpisodicBrainNonZero_beta.mat']));
    betaBetween = beta.betaMat(:,:,[1 3 4 6 7 9]); % Corresponding to people and places
    betaWithin = beta.betaMat(:,:,[4 6]);

    betaLangloc = load(fullfile(analysisDir, ...
        ['sub-' subject '_space-' templateSpace '_desc-langlocBrainNonZero_beta.mat']));
    betaLangloc = betaLangloc.betaMat;

    % if exist(outputPath, 'file'), continue; end

    % Load brain mask
    maskPath = fullfile(subjDir, 'func', ...
        ['sub-' subject '_space-' templateSpace '_desc-brainNonZero_mask.nii.gz']);

    % brainNonZero voxel mapping
    brainMaskVec = fpp.util.readDataMatrix(maskPath); 
    brainIdx = (brainMaskVec == 1);                   % logical
    brainVoxIdx = find(brainIdx);

    for r = 1:nROI
        roiName = roiMaskNames{r};  

        % ROI mask
        roiPath = fullfile(anatDir, ...
            ['sub-' subject '_space-' templateSpace '_desc-' roiName '_mask.nii.gz']);
        roiVec = fpp.util.readDataMatrix(roiPath);
        roiInBrain = (roiVec == 1) & brainIdx;

        % ROI indices in beta space
        fullIdx = find(roiInBrain);
        [~, roiRows] = intersect(brainVoxIdx, fullIdx, 'stable');

        % ------ Run SVM ------
        X_Within = betaWithin(roiRows,:)';
        [withinAcc(s,r), withinConf{s,r}] = runLOROSVM(X_Within, taskLabelsWithin, runLabelsWithin);

        X_Between = betaBetween(roiRows,:)'; 
        [betweenAcc(s,r), betweenConf{s,r}] = runLOROSVM(X_Between, taskLabelsBetween, runLabelsBetween);

        X_Langloc = betaLangloc(roiRows,:)';
        [langlocAcc(s,r), langlocConf{s,r}] = runLOROSVM(X_Langloc, taskLabelsLangloc, runLabelsLangloc);

        fprintf('\n Sub %s ROI %s | Langloc: %.3f | Within: %.3f | Between: %.3f\n', ...
            subject, roiName, langlocAcc(s,r), withinAcc(s,r), betweenAcc(s,r));

        % ------ Build null distributions ------ 
        disp([newline subject ' ' roiName ': Starting permutation test'])
        parfor p=1:nPerm
            % Within
            [permAccWithin(s,r,p),~] = runLOROSVM(X_Within, yPermWithin(p,:), runLabelsWithin);
            % Between
            [permAccBetween(s,r,p),~] = runLOROSVM(X_Between, yPermBetween(p,:), runLabelsBetween);
            % Langloc
            [permAccLangloc(s,r,p),~] = runLOROSVM(X_Langloc, yPermLangloc(p,:), runLabelsLangloc);
        end 
    end
end

%% Stats
meanAccWithin = mean(withinAcc(:));
meanAccBetween = mean(betweenAcc(:));
meanAccLangloc = mean(langlocAcc(:));
nullWithin  = nan(nPerm,1);
nullBetween = nan(nPerm,1);
nullLangloc = nan(nPerm,1);

for p = 1:nPerm
    nullWithin(p)  = mean(permAccWithin(:,:,p),  'all');
    nullBetween(p) = mean(permAccBetween(:,:,p), 'all');
    nullLangloc(p) = mean(permAccLangloc(:,:,p), 'all');
end

nullWithin = vertcat(nullWithin, meanAccWithin);
nullBetween = vertcat(nullBetween, meanAccBetween);
nullLangloc = vertcat(nullLangloc, meanAccLangloc);

% one-sided permutation p-values
pValWithin  = mean(nullWithin  >= meanAccWithin);
pValBetween = mean(nullBetween >= meanAccBetween);
pValLangloc = mean(nullLangloc >= meanAccLangloc);

disp(['Langloc p: ', num2str(pValLangloc) ' | Within p: ', num2str(pValWithin) ...
    ' | Between p: ', num2str(pValBetween)])

%% -----------------------------
%  Return results
%% -----------------------------
results = struct();

results.withinAcc  = withinAcc;
results.betweenAcc = betweenAcc;
results.langlocAcc = langlocAcc;

results.meanWithin = meanAccWithin;
results.meanBetween = meanAccBetween;
results.meanLangloc = meanAccLangloc;

results.nullWithin = nullWithin;
results.nullBetween = nullBetween;
results.nullLangloc = nullLangloc;

results.withinConf = withinConf;
results.betweenConf = betweenConf;
results.langlocConf = langlocConf;

results.pValWithin = pValWithin;
results.pValBetween = pValBetween;
results.pValLangloc = pValLangloc;

save(outputPath, 'results')

%% -----------------------------
%  Visualization
%% -----------------------------

% Boxplot
box = figure('Position',[500 500 500 500]);

outputPath = [fppDir '/group/langHippoc/identMVPA_SVM_Results.mat'];
results = load(outputPath);
results = results.results;

N = numel(results.withinAcc);
dataMat = NaN(N, 3);

dataMat(1:numel(results.langlocAcc), 1) = results.langlocAcc(:);
dataMat(1:numel(results.withinAcc), 2)  = results.withinAcc(:);
dataMat(1:numel(results.betweenAcc), 3) = results.betweenAcc(:);

boxplot(dataMat, 'Labels', {'Langloc Accuracy', 'Within Task Accuracy','Between Task Accuracy'}, ...
    'Whisker', 1, 'Colors', 'k', 'BoxStyle', 'outline');
ylabel('Accuracy');
ylim([0 1]);
title('Decoding Accuracy');

% boxColors = [0.9 0.4 0.4; 0.3 0.6 0.95]; % blue red
% 
% hold on;
% 
% boxes = findobj(gca, 'Tag', 'Box'); 
% for i = 1:length(boxes) 
%     patch(get(boxes(i),'XData'), get(boxes(i),'YData'), ... 
%         boxColors(i,:), ... 
%         'FaceAlpha', 0.5, ... 
%         'EdgeColor', 'k'); 
% end
% 
% hold off;

boxOutputPath = '/mnt/sml_share/IDENT/derivatives/fpp/group/langHippoc/MVPA_Boxplot.png';
exportgraphics(box, boxOutputPath, 'BackgroundColor','white')

%% Confusion matrix
confMat = figure('Position',[500 500 900 500]);

% Create group confusion matrices
% Langloc
C = size(results.langlocConf{1,1}, 1);
groupConfLangloc = zeros(C, C);

for s = 1:nSubj
    for r = 1:nROI
        groupConfLangloc = groupConfLangloc + results.langlocConf{s,r};
    end
end
groupConfLangloc = groupConfLangloc / sum(groupConfLangloc(:));   % normalize

% Within
C = size(results.withinConf{1,1}, 1);
groupConfWithin = zeros(C, C);

for s = 1:nSubj
    for r = 1:nROI
        groupConfWithin = groupConfWithin + results.withinConf{s,r};
    end
end
groupConfWithin = groupConfWithin / sum(groupConfWithin(:));   % normalize

% Between
C = size(results.betweenConf{1,1}, 1);
groupConfBetween = zeros(C, C);

for s = 1:nSubj
    for r = 1:nROI
        groupConfBetween = groupConfBetween + results.betweenConf{s,r};
    end
end
groupConfBetween = groupConfBetween / sum(groupConfBetween(:));   % normalize

% Save group confusion matrices to results
save(outputPath, '-append', 'groupConfLangloc','groupConfWithin','groupConfBetween')

% ------ Langloc  -------
taskNamesLangloc = {'Sentences','Nonwords'};
nTasks = length(taskNamesLangloc);
allGroupConf = [groupConfLangloc; groupConfWithin; groupConfBetween];
climits = [min(allGroupConf(:)), max(allGroupConf(:))];

subplot(1,3,1);
imagesc(groupConfLangloc);
clim(climits);
colorbar;
axis square;
set(gca,'XTick',1:nTasks,'YTick',1:nTasks,...
    'XTickLabel',[],'YTickLabel',[]);
% title('Langloc Confusion Matrix');
% ylabel('True');
% xlabel('Predicted');
colormap(parula);

% ------ Within Task -------
taskNamesWithin = {'People','Places'};
nTasks = length(taskNamesWithin);

subplot(1,3,2);
imagesc(groupConfWithin);
clim(climits);
colorbar;
axis square;
set(gca,'XTick',1:nTasks,'YTick',1:nTasks,...
    'XTickLabel',[],'YTickLabel',[]);
% xlabel('Predicted');
% title('Within Task Confusion Matrix');
colormap(parula);

% ------ Between task -------
taskNamesBetween = {'People','Places'};
nTasks = length(taskNamesBetween);

subplot(1,3,3);
imagesc(groupConfBetween);
axis square;
clim(climits)
colorbar;
set(gca,'XTick',1:nTasks,'YTick',1:nTasks,...
    'XTickLabel',[],'YTickLabel',[]);
%title('Between Task Confusion Matrix');
% xlabel('Predicted');
colormap(parula);

confOutputPath = '/mnt/sml_share/IDENT/derivatives/fpp/group/langHippoc/MVPA_ConfMat.png';
exportgraphics(confMat, confOutputPath, 'BackgroundColor', 'white')

%% ========================================================================
%  Helper function: Leave-One-Out SVM
%% ========================================================================
function [acc,confMat] = runLOROSVM(X,y,runs)

uniqueRuns = unique(runs);
nRuns = numel(uniqueRuns);

yPred = nan(size(y));

for r = 1:nRuns

    testIdx  = runs == uniqueRuns(r);
    trainIdx = runs ~= uniqueRuns(r);

    Xtrain = X(trainIdx,:);
    Xtest  = X(testIdx,:);
    ytrain = y(trainIdx);

    % Normalize using training data only
    mu = mean(Xtrain);
    sigma = std(Xtrain);
    sigma(sigma==0) = 1;

    Xtrain = (Xtrain - mu) ./ sigma;
    Xtest  = (Xtest  - mu) ./ sigma;

    % Train with C tuning
    model = trainLinearSVM_Cgrid(Xtrain,ytrain);

    % Predict
    yPred(testIdx) = predict(model,Xtest);

end

acc = mean(yPred == y);
confMat = confusionmat(y,yPred);

end

% Multiclass SVM 
% for t = 1:nRuns
%     train = true(nRuns,1); train(t) = false;
%     model = fitcecoc(X(train,:), y(train), ...
%         'Learners', templateSVM('KernelFunction','linear'), ...
%         'Coding', 'onevsall',...
%         'ClassNames', classes);
% 
%     yPred(t) = predict(model, X(~train,:));
% end

%% C tuning function

function model = trainLinearSVM_Cgrid(Xtrain,ytrain)

Cs = logspace(-3,3,7);   % C values to test
bestC = Cs(1);
bestAcc = -inf;

% Use stratified CV to preserve class balance
cv = cvpartition(ytrain,'KFold',length(ytrain));

for c = 1:numel(Cs)
    Cval = Cs(c);
    accFold = nan(cv.NumTestSets,1);

    for k = 1:cv.NumTestSets

        tr = training(cv,k);
        te = test(cv,k);

        % Train
        mdl = fitcsvm(Xtrain(tr,:), ytrain(tr), ...
            'KernelFunction','linear', ...
            'BoxConstraint', Cval, ...
            'Standardize', false);

        % Predict
        yhat = predict(mdl, Xtrain(te,:));

        % Compute accuracy
        accFold(k) = mean(yhat == ytrain(te)');
    end

    mAcc = mean(accFold);

    if mAcc > bestAcc
        bestAcc = mAcc;
        bestC = Cval;
    end
end

% Retrain final model on all training data
model = fitcsvm(Xtrain,ytrain, ...
    'KernelFunction','linear', ...
    'BoxConstraint', bestC, ...
    'Standardize', false);

end