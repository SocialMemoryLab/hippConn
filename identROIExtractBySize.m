function identROIExtractBySize(contrastName, searchName, isCifti, varargin)

fppDir = '/mnt/sml_share/IDENT/derivatives/fpp';
subjects = {'ident01','ident02','ident03','ident04','ident05',...
            'ident06','ident07','ident08','ident09','ident10'};

% Variable argument definitions 
roiSize = 5:5:40;
sizeType = 'pct';
overwrite = 0;

varArgList = {'invertStats','isCifti','overwrite'};
for i=1:length(varArgList)
    argVal = fpp.util.optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end
end

% Handle input desc
switch contrastName
    case 'PersonVsPlace'
        defineROITask = 'famsemantic';
        contrast = [1 0 -1];
    case 'BeliefVsPhoto'
        defineROITask = 'tomloc2';
        contrast = [1 -1];
    case 'SentencesVsNonwords'
        defineROITask = 'langloc';
        contrast = [1 -1];
    otherwise
        error('Unknown contrast name');
end

if ~exist('isCifti','var') || isempty(isCifti)
    isCifti = 1;
end
if isCifti
    searchSpace = 'fsLR_den-32k';
    statSpace = 'individual_res-2_den-32k';
    imageExt = '.dscalar.nii';
else
    searchSpace = 'session';   
    statSpace = 'session';
    imageExt = '.nii.gz';
end
if invertStats
    invertStats = 1;
    contrast = -1 .* contrast;
    invertSuffix = 'Inverted';
else
    invertSuffix = '';
end 

% Define output paths 
groupDir = [fppDir '/group'];
if ~exist(groupDir,'dir'), mkdir(groupDir); end

outputDesc = [searchName contrastName invertSuffix 'BySizeN' num2str(length(subjects))];

outputPath = [groupDir '/langHippoc/space-' statSpace '_desc-' outputDesc '_roiData.mat'];
outputFigurePath = [groupDir '/langHippoc/space-' statSpace '_desc-' outputDesc '_responseplot.png'];

if exist(outputPath,'file') && ~overwrite
    return;
end

%% ---- Initialize ----
psc = [];
subNums = [];
runNums = [];
pscBySub = [];

% Extract psc
for s = 1:length(subjects)

    subject = subjects{s};
    subjDir = [fppDir '/sub-' subject];
    anatDir = [subjDir '/anat'];
    roiDir = [subjDir '/roi'];
    analysisDir = [subjDir '/analysis'];
    funcDir = [subjDir '/func'];

    maskPath = [funcDir '/sub-' subject '_space-' statSpace '_desc-brainNonZero_mask' imageExt];

    % if ~isCifti
    %     subjStr = ['sub-' subject '_'];
    % end
    subjStr = ['sub-' subject '_'];

    searchPath = [anatDir '/' subjStr ...
        'space-' searchSpace '_desc-' searchName '_mask' imageExt];

    defineROIDir = [analysisDir '/sub-' subject ...
        '_task-' defineROITask ...
        '_run-01_space-' statSpace '_modelarma'];

    extractResponseDir = defineROIDir;

    conTmp = [];

    for r = 1:length(roiSize)

        [pscTmp,~,runNames] = fpp.func.roiExtract( ...
            extractResponseDir,...
            defineROIDir,...
            contrastName,...
            searchPath,...
            'maskPath',maskPath,...
            'roiSize',roiSize(r),...
            'sizeType',sizeType,...
            'invertStats',invertStats,...
            'overwrite',overwrite);

        conTmp = [conTmp pscTmp*contrast'];
    end

    psc = [psc; conTmp];
    runNums = [runNums; cellfun(@str2num,runNames)'];
    subNums = [subNums; s*ones(length(runNames),1)];
    pscBySub(s,:) = mean(conTmp);

    disp(['Extracted: sub-' subject]);

end

%SE
pscRunwiseStdErr = std(psc)/sqrt(size(psc,1));

%Save 
save(outputPath,'psc','runNums','subNums','pscBySub',...
    'subjects','pscRunwiseStdErr');

%% PLOT 

addpath(genpath('/home/aagarwal2/Desktop/mseb'));
savepath; 

figure('Position',[200 200 800 400]);
xTicks = 5:5:40;
linProps.col = {'r'};

yMax = max(mean(psc) + pscRunwiseStdErr);
ylim([0 yMax]);

mseb(xTicks,mean(psc),pscRunwiseStdErr,linProps,1);
set(gca,'linewidth',2,'FontSize',16);
set(gca,'XTickLabel',xTicks,'XTick',xTicks);
set(gcf,'Color',[1 1 1]);
title(searchName, contrastName);
ax = gca;
ax.FontSize = 10;      % tick numbers
xlabel('ROI Size (%)', 'FontSize', 12);
ylabel('Percent Signal Change', 'FontSize',12);

saveas(gcf,outputFigurePath);

end