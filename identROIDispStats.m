
% identROIDispStats(defineROITask,contrastName,searchName,isCifti,varargin)
%
% Function to display p-values for 
%
% Arguments:
% - defineROITask (string): name of task used to define ROIs
% - contrastName (string): name of contrast used to define ROIs
% - searchName (string): name of search space used to define ROIs
% - isCifti (boolean): whether to use CIFTI inputs, or volumetric
% 
% Variable arguments:
% - roiSize (scalar): size of ROI, in % or # of coords in a search space
% - sizeType ('pct' or 'num'): whether size is measured in # or % of coords
% - invertStats (boolean, default = 0): whether to invert statistical map
% - inputSuffix (string): suffix for input modelperm directories
%
% NOTE: Conditions 4-5 (FamPlace, UnfamPlace) of famvisual are switched in
% order in the outputs, so that FamPlace is the last condition.

function [statTable] = identROIDispStats(defineROITask,contrastName,searchName,isCifti,varargin)

% Define constants
fppDir = '/mnt/sml_share/IDENT/derivatives/fpp';
groupDir = [fppDir '/group/langHippoc'];
subjects = {'ident01','ident02','ident03','ident04','ident05','ident06','ident07','ident08','ident09','ident10'};

lmeNames = {'famvisualPersonVsPlace','famvisualPersonVsObject','famvisualFamVsUnfamPerson',...
    'famvisualPersonVsPlaceInverted','famvisualPlaceVsObject','famvisualFamVsUnfamPlace',...
    'famsemanticPersonVsPlace','famsemanticPersonVsObject',...
    'famsemanticPersonVsPlaceInverted','famsemanticPlaceVsObject',...
    'famepisodicPersonVsPlace','famepisodicPersonVsObject',...
    'famepisodicPersonVsPlaceInverted','famepisodicPlaceVsObject',...
    'dylocPersonVsPlace','dylocPersonVsObject',...
    'dylocPersonVsPlaceInverted','dylocPlaceVsObject',...
    'tomloc2BeliefVsPhoto','langlocSentencesVsNonwords',...
    'famvisualPersonVsObjectInverted','famvisualPlaceVsObjectInverted',...
    'famsemanticPersonVsObjectInverted','famsemanticPlaceVsObjectInverted',...
    'famepisodicPersonVsObjectInverted','famepisodicPlaceVsObjectInverted',...
    'dylocPersonVsObjectInverted','dylocPlaceVsObjectInverted'};

if ~exist('isCifti','var') || isempty(isCifti)
    isCifti = 1;
end
if isCifti
    templateSpace = 'individual_res-2_den-32k';
else
    templateSpace = 'session';
end

% Variable argument defaults
roiSize = 5;
sizeType = 'pct';
inputSuffix = '';

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'roiSize','sizeType','invertStats','inputSuffix'};
for i=1:length(varArgList)
    argVal = fpp.util.optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end
end

% Define output suffices based on variable arguments
if strcmpi(sizeType,'num')
    roiSize = round(roiSize);
    numSuffix = ['Top' num2str(roiSize)];
else
    numSuffix = ['Top' num2str(roiSize) 'Pct'];
end

% Subject number suffix
subSuffix = ['N' num2str(length(subjects))];

% Define outputs
outputDesc = [searchName defineROITask inputSuffix contrastName numSuffix subSuffix];
outputPath = [groupDir '/space-' templateSpace '_desc-' outputDesc '_roiData.mat'];
load(outputPath);

statsToPlot = [];
if ismember(contrastName,{'PersonVsPlace','BeliefVsPhoto','FacesVsScenes',...
        'SentencesVsNonwords'})
    statsToPlot = [1:3 7 8 11 12 15 16 19 20];
elseif ismember(contrastName,{'PersonVsPlaceInverted','FacesVsScenesInverted'})
    statsToPlot = [4:6 9 10 13 14 17 18 19 20];
elseif ismember(contrastName,{'PlaceVsObjectInverted','PersonVsObjectInverted','FacesVsObjectsInverted',...
        'ScenesVsObjectsInverted'})
    statsToPlot = 19:28;
end

statArray = nan([length(lmeNames), 2]);

for i=statsToPlot
    for p=1:length(lmeOneTailedP)
        if lmeOneTailedP{p}==0
            lmeOneTailedP{p} = lmeResults{2}.Coefficients.pValue/2;
        end
    end 

    tStat = lmeResults{i}.Coefficients.tStat;
    df = lmeResults{i}.Coefficients.DF;
    pVal = lmeOneTailedP{i};
    statArray(i,1) = tStat;
    statArray(i,2) = pVal;
    % disp([lmeNames{i} ': t(' num2str(df) ') = ', num2str(tStat) ', P = ' num2str(pVal)]);  
end  

statTable = array2table(statArray);
statTable = rmmissing(statTable);
statTable.Properties.RowNames = lmeNames(statsToPlot);
statTable.Properties.Description = [searchName defineROITask contrastName];

end