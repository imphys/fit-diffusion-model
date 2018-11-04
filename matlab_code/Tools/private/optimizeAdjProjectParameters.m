function [alignpar , fval, data] = optimizeAdjProjectParameters( predicted , data, project , alignpar_in, opts, noiseLevel, paramprior)
%[alignpar , fval,  data] = adjustProjectParameters( predicted , data,  project , alignpar_in, opts);
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
%                        [logpdf, dLogPDF_dData, d2Logpdf_dData2] = logPDFfun( data, predicted , opts.noiseLevel, [true false false])
%        noiseLevel    : current noise level estimate.
%
%
% OUTPUT
%   alignpar : optimized aligment parameters
%   fval     : final - log likelihood of the fit. 
%   data     = project( data, alignpar )
%
% Created by Dirk Poot, Erasmus MC, 1-4-2011
opt = optimset('fminunc');
hessmulfun = @(hessinfo, Y) hessinfo.hessMulFun(hessinfo, full(Y)); % The input Y might be sparse if optimization includes just 1 parameter.
opt = optimset(opt, 'gradObj', 'on' ,'Hessian','on','Display','iter', 'MaxIter', opts.maxIter, 'TolFun',opts.tolFun, 'largescale','on','PrecondBandWidth',inf,'HessMult',hessmulfun);
logpdffun = @(data) opts.logPDFfun(data, predicted, noiseLevel, [true false false]); 
optfun = @(alignpar) adjustAdjProjectParametersCritFun(  data, logpdffun, project , alignpar, paramprior) ;
if isempty(alignpar_in)
    % dont align when there are no parameters to adjust the alignment.
    alignpar = alignpar_in;
    fval = optfun( alignpar );
else
    [alignpar, fval,exflag,outp] = fminunc( optfun , alignpar_in, opt);
end;
if nargout>=3
    data = project{1}(  data, alignpar );
end;

function [f, g, hinfo] = adjustAdjProjectParametersCritFun( data, logpdffun, project , alignpar, paramprior)
projout = cell(1,min(nargout,2)); % do not request hessian from project.
[projout{:}] = project{1}(  data, alignpar );
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
    % dLL/dpar = dLL/dA dA/dpar
    % dLL_i/dA_j == delta_ij dLL_i/dA_i
    % dA_i/dpar_j == projout{2}(i,j)
    gradientParAdjMul = project{3}; 
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
            lgpdfout{3} = lgpdfout{3}(:);
            gradientParMul =  project{2};

            hinfo.hessMulFun  = @adjustAdjProjectParametersCritFun_hessmul;
            hinfo.gradientParMul =gradientParMul;
            hinfo.gradientParAdjMul = gradientParAdjMul;
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
        % additionally compute full hessian:

        % D2_LL_i/dpar_j dpar_k = d2LL_i/dA_idA_i  dA_i/dpar_j dA_i/dpar_k   (+ dLL_i/dA_i  d2A_i/dpar_j dpar_k )
%         lgpdfout{3} = lgpdfout{3}(:);
%        
%         h = - (projout{2}'* bsxfun(@times , lgpdfout{3} , projout{2}));
%         
        hinfo.logPDFhess  = -lgpdfout{3};
        hinfo.projectGrad = projout{2};
        hinfo.logPDFgrad  = -lgpdfout{2};
        hinfo.projected   = projout{1};
    end;
    
end;



function [Hx] = adjustAdjProjectParametersCritFun_hessmul(hessinfo, x)
if size(x,2)~=1
    Hx = x;
    for k=1:size(x,2)
        Hx(:,k) =  adjustAdjProjectParametersCritFun_hessmul(hessinfo, x(:,k));
    end;
    return;
end;
[tmp , hessinfo.projectGrad ] = hessinfo.gradientParMul( x , hessinfo.projectGrad );
[Hx                         ] = hessinfo.gradientParAdjMul( tmp(:).*hessinfo.logPDFhess, hessinfo.projectGrad );
if ~isempty(hessinfo.hpriormul)
    Hx = Hx - hessinfo.hpriormul(hessinfo.hprior, x );
end;