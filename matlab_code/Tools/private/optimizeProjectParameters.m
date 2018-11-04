function [alignpar , fval, projout] = optimizeProjectParameters( data, predicted, project , alignpar_in, logPDFfun, noiseLevel, opts, paramprior)
%[alignpar , fval,  data] = adjustProjectParameters( predicted , data,  project , alignpar_in, opts, paramprior);
% Adjusts the alignment (project adjust) parameters to improve the fit. 
% INPUTS
%   predicted : ND array with predicted magnitudes
%   data   : (image) data for project function
%   project   : [data , dData_dAlignpar ] = project( data, alignpar )
%               function that warps the data to ND MR images, dependent on alignpar
%   alignpar_in: initial vector for alignpar
%   opts  : option structure with (at least) the fields 
%        maxIter_align : maximum number of iterations.
%        tolFun_align  : stopping tolerance.
%        logPDFfun  : log likelihood function. Measures distance (at each voxel) between predicted and data.
%                        [logpdf, dLogPDF_dA, d2Logpdf_dA2] = logPDFfun( data, predicted , opts.noiseLevel , [false true false])
%        noiseLevel    : current noise level estimate.
%   paramprior : 1x1 or 1x2 cell array with functions specifying the log likelihood of the prior for the paramaters
%               { logpriorfun [, hessmul] }
%
% OUTPUT
%   alignpar : optimized aligment parameters
%   fval     : final - log likelihood of the fit. 
%   data     = project( data, alignpar )
%
% Created by Dirk Poot, Erasmus MC, 1-4-2011
logpdffun = @(A) logPDFfun(data, A, noiseLevel , [false true false]); 

optfun = @(alignpar) adjustProjectParametersCritFun( predicted , logpdffun, project, alignpar, paramprior) ;
if isempty(alignpar_in) || opts.skip
    % dont align when there are no parameters to adjust the alignment.
    alignpar = alignpar_in;
    if nargout>=2
        if nargout>=3
            [fval, dum, projout] = optfun( alignpar );
            clear optfun logpdffun logPDFfun
%             projout = rmfield(projout, {'hessMulFun','h'});
        else
            fval = optfun( alignpar );
        end;
    end;
else
    opt = optimset('fminunc');
    hessmulfun = @(hessinfo, Y) hessinfo.hessMulFun(hessinfo, full(Y)); % The input Y might be sparse if optimization includes just 1 parameter.
%     makePreconditioner = @(hessinfo , upperbandw, DM, DG) hessinfo.makePreconditioner( hessinfo , upperbandw, DM, DG );

%, 'Preconditioner',@(A,upbw, DM, DG) aprecon_full(A.h, upbw, DM, DG)
    opt = optimset(opt, 'gradObj', 'on' ,'Hessian','on','Display','iter', 'MaxIter', opts.maxIter, 'TolFun',opts.tolFun, 'largescale','on','PrecondBandWidth',inf,'HessMult',hessmulfun);
    [alignpar, fval, exflag,outp, grad, projout] = fminunc( optfun , alignpar_in, opt);
end;

function [f, g, hinfo] = adjustProjectParametersCritFun( predicted , logpdffun, project , alignpar, paramprior)
% [f, g, h] = adjustProjectParametersCritFun( predicted , logpdffun, project , alignpar)
% Criterium function for optimization of projection adjustment paramters (Alginment..)
% Created by Dirk Poot, Erasmus MC, 18-10-2011

projout = cell(1,min(max(1,nargout),2)); % do not request hessian from project.
[projout{:}] = project{1}(  predicted, alignpar );

lgpdfout = cell(1, nargout);
[lgpdfout{:}] = logpdffun( projout{1} );
f = - sum( lgpdfout{1}(:) );
if ~isempty(paramprior{1})
    lgprior = cell(1, nargout);
    [lgprior{:}] = paramprior{1}( alignpar );
    f = f - lgprior{1};
end;
if nargout>=2
    % also compute gradient:
    % f = - sum_i LL_i
    % d f/ dpar_k = -sum_i  dLL_i/dpar_k = -sum_i dLL_i/dA dA/dpar_k
    %   = -sum_i dLL_i/dA_j dA_j/dpar_k == -sum_i (delta_ij dLL_i/dA_i)  dA_j/dpar_k
    %   = -sum_i dLL_i/dA_i dA_i/dpar_k == - lgpdfout{2}(:)' * projout{2}(:,k);
    gradientParAdjMul = project{5}; 
    if isempty(gradientParAdjMul) % hack function to allow calling without parameter gradients function, just to return hinfo.
        g=[];
    else
    [g , projout{2}] = gradientParAdjMul( lgpdfout{2} , projout{2} ); 
    g = -g;
    if ~isempty(paramprior{1})
        g = g - lgprior{2}(:);
    end;
    if nargout>=3
        % additionally compute full hessian:

        % - sum_i d2 LL_i/d par_j d par_k = -sum_i d2 LL_i/dA_i d A_i * dA_i/dpar_j dA_i/dpar_k   (+ dLL_i / dA_i  d2A_i/dpar_j dpar_k )
        explicithessian = false;
        lgpdfout{3} = lgpdfout{3}(:);
        gradientParMul =  project{4};
        if explicithessian
            savememory = true;
            if savememory
                h = zeros(numel(alignpar));
                for k=1:numel(alignpar)
                    dpar = alignpar; dpar(:)=0; dpar(k)=1; % preserve type of alignpar, set all to zero, except element k
                    [tmp , projout{2}] = gradientParMul( dpar , projout{2});
                    [h(:,k) , projout{2}]  = gradientParAdjMul( tmp(:).*lgpdfout{3}, projout{2} );
                end;
            else
                dpar = eye(numel(alignpar), numel(alignpar), class(alignpar));
                [tmp , projout{2}] = gradientParMul( dpar , projout{2});
                h  = gradientParAdjMul( bsxfun(@times, tmp , lgpdfout{3}), projout{2} );
            end;
            hinfo.h = -h;
            hinfo.hessMulFun  = @adjustProjectParametersCritFun_hessmul;
        else
            hinfo.hessMulFun  = @adjustProjectParametersCritFun_hessmul2;
            hinfo.gradientParMul =gradientParMul;
            hinfo.gradientParAdjMul = gradientParAdjMul;
        end;
        if ~isempty(paramprior{1})
            hinfo.hprior = lgprior{3};
            if size(paramprior,2)==1
                hinfo.hpriormul = @(Hinfo, x) Hinfo*x;
            else
                hinfo.hpriormul = paramprior{2};
            end;
        else
            hinfo.hpriormul = [];    
        end;
    end;
    end;
    if nargout>=3
        hinfo.logPDFhess  = -lgpdfout{3};
        hinfo.projectGrad = projout{2};
        hinfo.logPDFgrad  = -lgpdfout{2};
        hinfo.projected   = projout{1};
    end;
end;

function [Hx] = adjustProjectParametersCritFun_hessmul(hessinfo, x)
Hx = hessinfo.h * x;
if ~isempty(hinfo.hpriormul)
    Hx = Hx - hessinfo.hpriormul(hessinfo.hprior,  x );
end;

function [Hx] = adjustProjectParametersCritFun_hessmul2(hessinfo, x)
if size(x,2)~=1
    Hx = x;
    for k=1:size(x,2)
        Hx(:,k) =  adjustProjectParametersCritFun_hessmul2(hessinfo, x(:,k));
    end;
    return;
end;
[tmp , hessinfo.projectGrad ] = hessinfo.gradientParMul( x , hessinfo.projectGrad );
[Hx                         ] = hessinfo.gradientParAdjMul( tmp(:).*hessinfo.logPDFhess, hessinfo.projectGrad );
if ~isempty(hessinfo.hpriormul)
    Hx = Hx - hessinfo.hpriormul(hessinfo.hprior, x );
end;