function [currentPoint,data] = optimizeAllAdjProjectParameters(currentPoint, opts, repeatGlobalEstimationlp)
%[currentPoint, data] = optimizeAllAdjProjectParameters( currentPoint, opts, repeatGlobalEstimationlp);
% Helper function for fit_MRI that optimizes all alignment parameters 
% (projectParameters) to improve the fit. Uses adjProject.
% (use optimizeAllProjectParameters to optimize project parameters)
%
% INPUTS
%   currentPoint : the full 'current point' specifying structure. 
%   opts         : the option structure
%   repeatGlobalEstimationlp: global loop counter. Optimization should be 
%                           skipped in the first loop.
% OUTPUT
%   currentPoint : the new current point, with the projectParameters
%                  (and projectGrad) updated.
%   data  : currently predicted data.
%
% Created by Dirk Poot, Erasmus MC, 5-4-2012

% adjProject data 
if (repeatGlobalEstimationlp <= opts.projectOptimization.skip)
    % do not align in first iteration (since we only have the initial values for theta)
    data = cell(1,numel(opts.data_in));
    for k=1:numel(opts.data_in);
        data{k} = opts.adjProject{k}( opts.data_in{k}, opts.projectParameters{k} );
    end;
else
    if opts.numPDFoptpar>0
        optsnoiselevelbackup = opts.noiseLevel;
        noiseLevel = currentPoint.theta(end-opts.numPDFoptpar+1:end,:,:,:,:);
    else
        noiseLevel = opts.noiseLevel;
    end;
    if isempty(currentPoint.predictedImages)
        currentPoint.predictedImages = predictAllImagesfun( currentPoint.theta, opts );
    end;
    % align after first round of ML fit:
%     logLikBeforeAlign = opts.logPDFfun(data, currentPoint.predictedImages, noiseLevel);
%     logLikBeforeAlign = sum(logLikBeforeAlign(:,:),2);
    logLikAfterAlign = zeros(1,numel(opts.data_in));
    szData = [opts.numImages opts.spatialSize];
    predictedImagesp = permute( currentPoint.predictedImages, [2:numel(szData) 1]);
    sel = repmat({':'},1,numel(szData));
    data = cell(1,numel(opts.data_in));
    progressbar('start',numel(opts.data_in));
    for k=1:numel(opts.data_in);
        sel{end} = opts.projectSelectSource{k};
        predictedImagesps = predictedImagesp(sel{:});
        if size(opts.projectParameterPrior,1)==1
            projectParameterPrior_k = opts.projectParameterPrior;
        else
            projectParameterPrior_k = opts.projectParameterPrior(k,:);
        end;
        noiseLevel = opts.noiseLevel;
        if size(noiseLevel,1)==opts.numImages
            noiseLevel = permute(noiseLevel(sel{[end 1:end-1]}),[2:numel(szData) 1]);
        end;
        [alignpar , fval, data{k}] = optimizeAdjProjectParameters( predictedImagesps , opts.data_in{k}, opts.adjProject(k,:) ,opts.projectParameters{k}, opts, noiseLevel, projectParameterPrior_k);
        currentPoint.projectParameters{k} = alignpar;
        logLikAfterAlign(k) = fval;
        progressbar(k);
    end;            
    progressbar('ready');
%     logLikAfterAlign = sum(logLikAfterAlign);
%     logLikBeforeAlign = sum(logLikBeforeAlign);
end;
data = permute(cat(numel(opts.spatialSize)+1, data{:}),[numel(opts.spatialSize)+1, 1:numel(opts.spatialSize)]);
