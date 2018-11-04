function [ll, grad, hessinfo] = voxelLLfun_sigma_m( par, fun, logpdffun, data, maxfunArgsOut, explicitHessian, logprior, npar_sigma, funHessian_I, funHessian_J)
% [ll, grad, hess] = voxelLLfun_sigma_m( par, fun, logpdffun, data , maxfunArgsOut, explicitHessian, logprior, npar_sigma)
% computes the likelihood of the MR data with likelihood function logpdffun.
% INPUTS:
% par           : a vector/matrix with the parameters for fun (each column of par specifies a voxel)
% fun(par)      : a function predicting the magnitudes of the MR images for each column of par
%                 (Each column specifies the parameters of one voxel)
%                 [A, dA, d2A] = fun( par )
%                 par : (npar  x numtr) input parameters with sigma part removed.
%                 A  : ( numim  x  numtr  ) , where nelpp is the number of components per datapoint (typically 1, but might be 2 for complex data)
%                 dA : ( numim  x  numtr x npar ) , derivative of A wrt to each element of par (of the corresponding voxel/column in par)
%                 d2A: ( numim  x  numtr x (npar*(npar+1)/2) ) , d2 A /d par(I) d par(J) with [I,J] = find(triu(ones(npar)));
% logpdffun(data, A, sigma) : log(pdf(data, A, sigma)), log likelihood function, e.g. @(data, A, sigma) logricepdf(data, A, sigma, [false, true, true])
% data          : Measurements; (magnitude) MR values.
% maxfunArgsOut : maximum number of output arguments of fun
% explicitHessian : boolean: should explicit or compacted hessian be returned?
% logprior.fun   : if non empty: function that computes the log likelihood of the prior distribution of par.
%                               as well as the gradient and compactified hessian.
%  
% npar_sigma    : number of parameters describing the pdf, passed as third element to logpdffun
%
% OUTPUTS:
% ll : - log likelihood value
% grad: gradient of ll 
% hess: hessian of ll
%
% handles multiple data vectors simultanuously, each column of par should correspond to a column in data (and par after reshape).
%
% NOTE: update code in conjunction with voxelLLfun_m
%
% Created by Dirk Poot, Erasmus MC, 22-3-2011
szdta = size(data);
numtr = prod(szdta(2:end));

par = reshape(par,[],numtr);
sigma = par(end-npar_sigma+1:end,:);
parr = par(1:end-npar_sigma,:);

fargs = cell(1, min( max(1,nargout) , maxfunArgsOut ));
[fargs{:} ] = fun( parr );
lrargs = cell(1, max(1,nargout));
[lrargs{:}] = logpdffun(data(:,:), fargs{1} , sigma, [false true true]); % derivative w.r.t. A and sigma.
ll = - sum( lrargs{1}(:) );
if ~isempty(logprior.fun)
    % Evaluate prior for the parameters:
    lpriorargs = cell(1, nargout);
    [lpriorargs{:}] = logprior.fun( par );
    ll = ll - sum(lpriorargs{1});
end;
if nargout>=2
    % also compute gradient:
    % d ll /d par_i =
    %  i<=nparrs     = sum_j  d ll / d lgpdf_j * d lgpdf_j/d f_j * d f_j / dpar_i
    %                  + d ll / d lgprior * d lgprior / d par_i
    %                =  -1 * sum( lrargs{2}(:,:,1) .* fargs{2}(:,:,i) ,1)'
    %                   -1 * lpriorargs{2}(:,i);
    %  i>nparrs      = sum_j  d ll / d lgpdf_j * d lgpdf_j/d dpar_i
    %                  + d ll / d lgprior * d lgprior / d par_i
    %                =  -1 * sum( lrargs{2}(:,:,i-nparrs) ,1)
    %                   -1 * lpriorargs{2}(:,i);
    numim = size(data,1);
    nparrs = size(parr,1);
    if size(fargs{2},1)==numim*numtr
        fargs{2} = reshape(fargs{2},[numim, numtr, nparrs]);
    end;
    lrargs{2} = reshape( lrargs{2}, numim, numtr, 1+npar_sigma);
    % compute gradient, see DifusionTensorLL_gradient.nb
    % evaluate dLLdA_i * dA_idpar
    % first part: derivative wrt to parr, where continued derivative is used
    % second part: derivative wrt to sigma. 
    grad = [reshape(-sum( bsxfun(@times, lrargs{2}(:,:,1), fargs{2}), 1), numtr, nparrs)'; reshape(-sum(lrargs{2}(:,:,2:end),1),numtr, npar_sigma)'];
    
    if ~isempty(logprior.fun)
        grad(logprior.gradient_linindex,:) = grad(logprior.gradient_linindex,:) - lpriorargs{2};
    end;
        
    if nargout>=3
        % additionally compute hessian:
        % D2_ll /d par_i d par_j = 
        % i,j<=nparrs   = sum_k d ll / d lgpdf_k * d2 lgpdf_k/d f_k d f_k  d f_k / dpar_i  d f_k / dpar_j
        %                      +d ll / d lgpdf_k * d lgpdf_k/d f_k * d2 f_k / d par_i d par_j
        %                 + d ll / d lgprior * d2 lgprior / d par_i d par_j    
        %               = -1 * sum(lrargs{3}(:,:,1) .* fargs{2}(:,:,i) .* fargs{2}(:,:,j)) 
        %                 -1 * sum( lrargs{2}(:,:,1) .* fargs{3}(:,:,ij) )
        %                 -1 * lpriorargs{3}(:,ij)  
        % i<=nparrs, j>nparrs =  sum_k  d ll / d lgpdf_k * d2 lgpdf_k/d f_k d dpar_j * d f_k / d par_i
        %                        + d ll / d lgprior * d2 lgprior / d par_i d par_j   
        %               =  -1 * sum( lrargs{3}(:,:,1(j-nparrs)) .* fargs{2}(:,:,i) ,1)
        %                  -1 * lpriorargs{3}(:,ij)  
        % i,j > nparrs  =  sum_k  d ll / d lgpdf_k * d2 lgpdf_k/d par_i d dpar_j 
        %                        + d ll / d lgprior * d2 lgprior / d par_i d par_j   
        %               =  -1 * sum( lrargs{3}(:,:,(i-nparrs)(j-nparrs)) ,1)
        %                  -1 * lpriorargs{3}(:,ij)  
        npars = size(par,1);
%         [I,J] = find(triu(ones(npars)));
        red1 = (funHessian_I<=nparrs) & (funHessian_J<=nparrs);
        Ir = funHessian_I(red1);
        Jr = funHessian_J(red1); 
        nIr = numel(Ir);
        lrargs{3} = reshape(lrargs{3}, numim,  numtr, (npar_sigma+2).*(npar_sigma+1)/2);
        hessLng = lrargs{3}(:,:,ones(1,nIr)) .* (fargs{2}(:,:,Ir).*fargs{2}(:,:,Jr));
        if maxfunArgsOut>=3 && ~isempty(fargs{3})
            % Main curvature typically comes from rice distribution, so allow ignoring of function curvature (by fargs{3}==[])
            % However, when provided, include hessian:
            if size(fargs{3},1)==numim*numtr
                fargs{3} = reshape(fargs{3},[numim, numtr, nparrs*(nparrs+1)/2]);
            end;
            hessLng = hessLng + lrargs{2}(:,:,ones(1,nIr)) .* fargs{3} ;
        end;
%         hesscoeff = reshape( -sum( reshape( hessLng, numim, numtr*nI) ,1) , numtr, nI)';
        hesscoeff1 = reshape( -sum( hessLng ,1) , numtr, nIr)';
        
        [Io, Jo] = find(triu(ones(1+npar_sigma)));
        red2 = (funHessian_I<=nparrs) & (funHessian_J>nparrs); % dL /d par d sigma
        lrind = find( (Io==1) & (Jo>1))';
        hessLng2 = lrargs{3}(:,:,ones(nparrs,1)*lrind) .* fargs{2}(:,:,(1:nparrs)'*ones(1,npar_sigma));
        hesscoeff2 = reshape( -sum( hessLng2,1) , numtr, nnz(red2))';
        
        red3 = (funHessian_I>nparrs) & (funHessian_J>nparrs); % dL /d sigma d sigma;    red3 could be simplified to I>nparrs, since I<=J
        hessLng3 = lrargs{3}(: , : , Io>1 ) ;
        hesscoeff3 = reshape( -sum( hessLng3,1) , numtr, nnz(red3))';
        
        hesscoeff = zeros(numel(funHessian_I),numtr);
        hesscoeff(red1,:) = hesscoeff1;
        hesscoeff(red2,:) = hesscoeff2;
        hesscoeff(red3,:) = hesscoeff3;
        
        if ~isempty(logprior.fun)
            hesscoeff(logprior.hessian_write_funhessindex,:) = hesscoeff(logprior.hessian_write_funhessindex,:) - lpriorargs{3};
        end;
        
%         hess = zeros(npars*numtr);
%         
%         indadj = ones(numel(I),1)*(npars*(0:numtr-1));
%         Ia = I(:,ones(1,numtr)) + indadj;
%         Ja = J(:,ones(1,numtr)) + indadj;
%         hess(Ia + (Ja-1)*npars*numtr) = hesscoeff;
%         hess(Ja + (Ia-1)*npars*numtr) = hesscoeff;
        if explicitHessian
            rev = funHessian_I~=funHessian_J;
            hessblk = [hesscoeff(:,:);hesscoeff(rev,:)];
            Ir = bsxfun(@plus, [funHessian_I;funHessian_J(rev)] , 0:npars:numel(par)-1);
            Jr = bsxfun(@plus, [funHessian_J;funHessian_I(rev)] , 0:npars:numel(par)-1);
            hessinfo = sparse(Ir,Jr, hessblk, numel(par), numel(par));
        else
            hessinfo.hesscoeff = hesscoeff;
            hessinfo.hessMulFun = @compactedHessianMul;
            hessinfo.makePreconditioner = @compactedHessianMakePreconditioner;            
        end;
    end;
end;
if ~isfinite(ll)
    ll = inf;
end;
