
% Function to compute resting-state partial correlations between a set of ROIs in
% IDENT data.
%
% VolSurf version: allow ROI files to be NIFTI volume files OR CIFTI
% surface/combined files. Check first for surface, then volume, then
% combined data.
%
% Arguments:
% - roiDescs (cell array): BIDS descriptions of ROIs to correlation
% - outputDesc (string): BIDS description for output file

function hippConnRestROIVolSurfPartialCorr(roiHippoc,roiCortex,outputDesc)

studyDir = '/mnt/sml_share/IDENT';
bidsDir = [studyDir '/derivatives/fpp'];
subjects = {'ident01','ident02','ident03','ident04','ident05','ident06',...
    'ident07','ident08','ident09','ident10'};
spaceStr = 'fsLR_den-32k';
spaceStrVol = 'session';
spaceStrComb = 'individual_res-2_den-32k';
nRuns = 6;
nVols = 300;
roiDescs = [roiHippoc,roiCortex];
nROIs = length(roiDescs);
inputDesc = 'preprocGSR';
outputPath = [bidsDir '/group/langHippoc/desc-' outputDesc '_ROICorrelationData.mat'];
corrMat = zeros(nROIs,nROIs,length(subjects));
isVol = zeros(length(roiDescs),1);

for s=1:length(subjects)
    
    % Define paths
    subject = subjects{s};
    subjDir = [bidsDir '/sub-' subject];
    roiDir = [subjDir '/roi'];
    anatDir = [subjDir '/anat'];
    funcDir = [subjDir '/func'];
    
    % Load ROIs
    roiMat = {};
    roiAvgMat = {};
    for r=1:nROIs   % set path to anat for LANG mask
        roiPath = [roiDir '/sub-' subject '_space-' spaceStr '_desc-' roiDescs{r} '_mask.dscalar.nii'];
        if ~exist(roiPath,'file')
            roiPath = [roiDir '/sub-' subject '_space-' spaceStrVol '_desc-' roiDescs{r} '_mask.nii.gz'];
            isVol(r) = 1;
            if ~exist(roiPath,'file')   
                roiPath = [anatDir '/sub-' subject '_space-' spaceStrVol '_desc-' roiDescs{r} '_mask.nii.gz'];
                isVol(r) = 1;
                if ~exist(roiPath,'file')
                    roiPath = [roiDir '/sub-' subject '_space-' spaceStrComb '_desc-' roiDescs{r} '_mask.dscalar.nii'];
                    isVol(r) = 2;
                    if ~exist(roiPath,'file')
                        error(['ROI file ' roiPath ' does not exist.']);
                    end
                end 
            end 
        end
        roiMat{r} = fpp.util.readDataMatrix(roiPath);
        roiAvgMat{r} = roiMat{r}/sum(roiMat{r});  % Dot product with this vector = mean across ROI
    end

    for r=1:nRuns
        restPath = [funcDir '/sub-' subject '_task-rest_run-' fpp.util.numPad(r,2)...
            '_space-' spaceStr '_desc-' inputDesc '_bold.dtseries.nii'];
        restPathVol = [funcDir '/sub-' subject '_task-rest_run-' fpp.util.numPad(r,2)...
            '_space-' spaceStrVol '_desc-' inputDesc '_bold.nii.gz'];
        restPathComb = [funcDir '/sub-' subject '_task-rest_run-' fpp.util.numPad(r,2)...
            '_space-' spaceStrComb '_desc-' inputDesc '_bold.dtseries.nii'];
        
        % Define outlier volumes
        confoundPath = fpp.bids.changeName(restPath,'desc','','confounds','.tsv');
        outlierPath = fpp.bids.changeName(restPath,{'space','res','den','desc'},...
            {'','','',''},'outliers','.tsv');
        outlierInd = fpp.util.readOutlierTSV(outlierPath);
        goodInd = setdiff(1:nVols,outlierInd);
        
        % Load resting data
        if sum(isVol==0)>0
            restMat = fpp.util.readDataMatrix(restPath);
            restMat = restMat(:,goodInd);
        end
        if sum(isVol==1)>0
            restMatVol = fpp.util.readDataMatrix(restPathVol);
            restMatVol = restMatVol(:,goodInd);
        end
        if sum(isVol==2)>0
            restMatComb = fpp.util.readDataMatrix(restPathComb);
            restMatComb = restMatComb(:,goodInd);
        end
        
        % Compute ROI-averaged time series
        roiAvgSeries = [];
        for r2=1:nROIs
            if isVol(r2)==0
                roiAvgSeries(:,r2) = restMat'*roiAvgMat{r2};
            elseif isVol(r2)==1
                roiAvgSeries(:,r2) = restMatVol'*roiAvgMat{r2};
            else
                roiAvgSeries(:,r2) = restMatComb'*roiAvgMat{r2};
            end
        end
        
        % Compute correlations
        corrMat(:,:,s) = corrMat(:,:,s) + partialcorr(roiAvgSeries)/nRuns;
        
        disp(['Computed correlation matrix for ' subject ', run ' int2str(r)]);
    end
end

save(outputPath,'corrMat','roiDescs');
end