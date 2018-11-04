function [LL, gradifo, grad, gradprojectpar] = projectAllAndLL( data, opts , prepgradhess, projectParameters)
% [LL, gradifo, grad] = projectAllAndLL( data, opts, prepgradhess, projectParameters )
% From all images (output of predictAllImagesfun)
% compute - log likelihood for each image.
% prepgradhess : switch that determines which info is put into gradifo.
%     0  : dont prepare gradient and/or hessian (into gradifo)
%     1  : just prepare gradient
%     2  : also prepare hessian.
% 
% INPUTS:
%  data : images predicted by fun 
%  opts : fit_MRI option structure; (uses data_in field as measured data)
%  prepgradhess : scalar integer (see above)
%  projectParameters : override for opts.projectParameters (to allow
%                      optimizing these)
% OUTPUTS:
%  LL : minus log likelihood of current point
%  gradifo : all gradient (and hessian) related information; a cell array
%            with user info (obtained from project and logPDF funcs).
%  grad : gradient of LL with respect to data
%  gradprojectpar : gradient with respect to projectPar.
%
% Created by Dirk Poot, Erasmus MC, 27-10-2011
%          12-3-2014: added projectParameters and gradprojectpar output 
if nargin<3 || isempty(prepgradhess)
    prepgradhess = 1;
end;
if nargin<4
    projectParameters = opts.projectParameters;
end;
if nargout>=2
    prepgradhess = max(1,prepgradhess);
end;
szData = size(data);
if isempty(opts.project)
    if prepgradhess>0
        error('Apparently this case is (unexpectedly) used; but its not implemented yet. (I can only use global optimization in fit_MRI when project is used).');
    end;
    % make code correct:
    [lgpdf] = opts.logPDFfun( opts.data_in, data , opts.noiseLevel, [false true false]);
    LL = - sum( lgpdf(:) );
else
    datap = permute( data, [2:numel(szData) 1]);
    sel = repmat({':'},1,numel(szData));
    LL = zeros(size(opts.project,1),1);
    projgradifo   = cell( size(opts.project,1) , min(1,prepgradhess) );
    lgpdfgradhess = cell( size(opts.project,1) , prepgradhess );
    if nargout>=3
        grad = cell(1,size(opts.project,1)); % gradient of cost function with respect to data
        if nargout>=4
            gradprojectpar = cell(1,size(opts.project,1));
        end;
    end;
    for k=1:size(opts.project,1)
        sel{end} = opts.projectSelectSource{k};
        if size(opts.noiseLevel,1)==1
            noiseLevel = opts.noiseLevel;
        else
            noiseLevel = opts.noiseLevel(k,:);
        end;
        if numel(noiseLevel)==1 && iscell(noiseLevel)
            noiseLevel = noiseLevel{1};
        end;
        [proj, projgradifo{k,:}] = project1(   datap(sel{:}) , projectParameters{k} , opts.project(k,:), 1);
        if iscell( opts.logPDFfun ) || iscell(opts.project{k,1})
            % feed part of data and project to each logPDFfun.
            lgpdf = zeros(size(proj));
            lgpdfgradhessk = cell( numel(proj) , size(lgpdfgradhess,2) );
            for pdfidx = 1:numel(lgpdf)
                if iscell( opts.logPDFfun )
                    pdffun = opts.logPDFfun{pdfidx};
                else
                    pdffun = opts.logPDFfun;
                end;
                if iscell(noiseLevel)
                    noiseLevelk = noiseLevel{1}(pdfidx);
                else
                    noiseLevelk = noiseLevel;
                end;
                [lgpdfk, lgpdfgradhessk{pdfidx,:}] = pdffun( opts.data_in{k}{ pdfidx }, proj{ pdfidx } , noiseLevelk, [false true false]);
                lgpdf(pdfidx) = sum(lgpdfk(:));
            end;
            for nargs = 1:size(lgpdfgradhess,2)
                lgpdfgradhess{k,nargs} = {lgpdfgradhessk{:,nargs}};
            end;
        else
            [lgpdf, lgpdfgradhess{k, :}] = opts.logPDFfun( opts.data_in{k}, proj , noiseLevel, [false true false]);
        end;
        LL(k) = - sum( lgpdf(:) );
        if nargout>=3
            % compute grad =  dLL/ d data_sel
            %          = dLL/d lgpdf * d lgpdf / dproj * dproj / d datap_sel
            % dLL/d lgpdf         = -1
            % d lgpdf / dproj     = lgpdfgradhess{1}
            % dproj / d datap_sel = project{k,3}( .., projgrad{1} )
            [grad{k} , projgradifo{k,1}] = project1( lgpdfgradhess{k,1} , projgradifo{k, 1} , opts.project(k,:), 3 );
            if numel(sel{end})==1
                % special case when we can use (faster) reshape instead of
                % permute
                grad{k} = reshape(-grad{k},[1 szData(2:end)]);
            else
                % general; permute
                grad{k} = permute( reshape(-grad{k},[szData(2:end) numel(sel{end})]), [numel(szData) 1:numel(szData)-1]);
            end;
            if nargout>=4
                [gradprojectpar{k} , projgradifo{k,1}] = project1( lgpdfgradhess{k,1} , projgradifo{k, 1} , opts.project(k,:), 5 );
                gradprojectpar{k} =-gradprojectpar{k};
            end;
        end;
    end;
    if nargout>=2
        gradifo.projectGradinfo = projgradifo;
        %gradifo.lgpdfgrad = {lgpdfgradhess{:,1}};
        if prepgradhess>=2
            gradifo.lgpdfhess = { lgpdfgradhess{:,2} };
        end;
        if nargout>=3
            if isequal([opts.projectSelectSource{:}],1:szData(1))
                % if each project takes 1 predicted image, we can construct
                % the grad output efficiently:
                grad = vertcat(grad{:});
            else
                % general case with arbitrary projectSelectSource:
                tmp = zeros(szData);
                for k=1:numel(grad)
                    tmp(opts.projectSelectSource{k},:) = tmp(opts.projectSelectSource{k},:) + grad{k}(:,:);
                end;
                grad = tmp;
            end;
        end;
    end;
end;