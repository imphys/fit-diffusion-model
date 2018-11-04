function [currentPoint] = optimizeAllProjectParameters(currentPoint, opts)
%[currentPoint] = optimizeAllProjectParameters( currentPoint, opts);
% Helper function for fit_MRI that optimizes all alignment parameters 
% (projectParameters) to improve the fit. Uses project. 
% (use optimizeAllAdjProjectParameters to optimize adjProject parameters)
%
% INPUTS
%   currentPoint : the full 'current point' specifying structure. 
%   opts         : the option structure
% OUTPUT
%   currentPoint : the new current point, with the projectParameters
%                  (and projectGrad) updated.
%
% Created by Dirk Poot, Erasmus MC, 5-4-2012


% if projectGrad is provided or not needed and project optimization should be skipped, we trust the projectGrad content and dont optimize.
if ~opts.projectOptimization.skip || (isempty(currentPoint.projectGrad) && opts.optimizeBlocks) 
    if isempty(currentPoint.predictedImages)
        currentPoint.predictedImages = predictAllImagesfun( currentPoint.theta, opts );
    end;
    szData = [opts.numImages opts.spatialSize];
    predictedImagesp = permute( currentPoint.predictedImages, [2:numel(szData) 1]); % TODO: can we avoid to permute this large(st) matrix?
    projectGrad = cell(numel(opts.data_in),1);
    sel = repmat({':'},1,numel(szData));
    % determine number of projection steps, for progressbar:
    cnt = 0;
    for k = 1 : size(opts.project,1)
        if iscell( opts.project{k,1} )
            cnt = cnt + numel(opts.project{k,1});
        else
            cnt = cnt + 1;
        end;
    end;
    progressbar('start', cnt);
    done = 0;
    % optimize all projection's:
    for k=1:size(opts.project,1)
        sel{end} = opts.projectSelectSource{k};
        predictedImagesps = predictedImagesp(sel{:});
        
        % Noise level specified per project?:
        if size(opts.noiseLevel,1)==1
            noiseLevel_k = opts.noiseLevel;
        else
            if isempty(opts.noiseLevel)
                error('An initial noise level has to be provided when fitting with projection.')
            end;
            sznoiselvl = size(opts.noiseLevel);
            noiseLevel_k = reshape(opts.noiseLevel(k,:), [1 sznoiselvl(2:end)]);
        end;
        if numel(noiseLevel_k)==1 && iscell(noiseLevel_k)
            noiseLevel_k = noiseLevel_k{1};
        end;
        % one or multiple project parameters provided?
        if size(opts.projectParameterPrior,1)==1
            projectParameterPrior_k = opts.projectParameterPrior;
        else
            projectParameterPrior_k = opts.projectParameterPrior(k,:);
        end;
        % one logPDFfun or for each project separately?
        if iscell( opts.logPDFfun )
            if numel(opts.logPDFfun)==1
                logPDFfun_k = opts.logPDFfun{1};
            else
                logPDFfun_k = opts.logPDFfun{k};
            end;
        else
            logPDFfun_k = opts.logPDFfun;
        end;        
        % Multiple projects per (block of) predicted image(s)?
        if iscell( opts.project{k,1} )
            nprojk = numel(opts.project{k,1});
            projectGrad{k} = struct('projectGrad',{cell(nprojk,1)}, 'logPDFhess',{cell(nprojk,1)},'logPDFgrad',{cell(nprojk,1)});
            for projidx = 1: nprojk
                % further specify noise level
                if iscell(noiseLevel_k)
                    noiseLevel_kj = noiseLevel_k{1}(projidx);
                else
                    noiseLevel_kj = noiseLevel_k;
                end;
                %further specify logPDF function
                if iscell( logPDFfun_k )
                    logPDFfun_kj = logPDFfun_k{projidx};
                else
                    logPDFfun_kj = logPDFfun_k;
                end;
                % get project:
                project_j = cell(1,size(opts.project,2));
                for idx = 1 : numel(project_j)
                    if iscell(opts.project{k,idx})
                        project_j{idx} = opts.project{k,idx}{projidx};
                    else
                        project_j{idx} = opts.project{k,idx};
                    end;
                end;
                if iscell(projectParameterPrior_k{1})
                    projectParameterPrior_kj = projectParameterPrior_k{1}(projidx,:);
                else
                    projectParameterPrior_kj = projectParameterPrior_k;
                end;
                [tmp1,tmp2,tmp3] = optimizeProjectParameters( opts.data_in{k}{projidx}, predictedImagesps, project_j  , currentPoint.projectParameters{k}{projidx} , logPDFfun_kj, noiseLevel_kj, opts.projectOptimization, projectParameterPrior_kj);
                currentPoint.projectParameters{k}{projidx} = tmp1;
                fval{k}(projidx) = tmp2;
                projectGrad{k}.projectGrad{projidx} = tmp3.projectGrad;
                projectGrad{k}.logPDFgrad{projidx}  = tmp3.logPDFgrad;
                projectGrad{k}.logPDFhess{projidx}  = tmp3.logPDFhess;
%                         projected{k}{projidx}               = tmp3.projected;
%                         par_upd = opts.projectParameters{k}{projidx} - par_init; % DEBUG
%                         disp(par_upd) % DEBUG
                done = done+1;progressbar(done);
            end;
        else
            [currentPoint.projectParameters{k}, fval(k), projectGrad{k}] = optimizeProjectParameters( opts.data_in{k}, predictedImagesps, opts.project(k,:) , currentPoint.projectParameters{k} , logPDFfun_k, noiseLevel_k, opts.projectOptimization, projectParameterPrior_k);
            done = done+1;progressbar(done);
        end;
    end;
    progressbar('ready');
    % DEBUG compactify projectGrad:
    % for k1=1:numel(projectGrad);for k2 = 1:numel(projectGrad{k1}.projectGrad);projectGrad{k1}.projectGrad{k2}.img_out = [];projectGrad{k1}.projectGrad{k2}.pargradient = [];end;end;
    % save dummy projectGrad
    currentPoint.projectGrad = projectGrad;
    clear predictedImagesp
end;