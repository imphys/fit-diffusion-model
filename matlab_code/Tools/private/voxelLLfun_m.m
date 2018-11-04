function [ll, grad, hessinfo] = voxelLLfun_m( par, fun, logpdffun, data, maxfunArgsOut, explicitHessian, logprior, funHessian_I, funHessian_J)
% [ll, grad, hess] = voxelLLfun_m( par, fun, logpdffun, data , maxfunArgsOut, explicitHessian, logpriorfun)
% computes the likelihood of the MR data with likelihood function logpdffun.
% INPUTS:
% par           : a vector/matrix with the parameters for fun (each column of par specifies a voxel)
% fun(par)      : a function predicting the magnitudes of the MR images for each column of par
%                 (Each column specifies the parameters of one voxel)
% logpdffun(data, A) : log(pdf(data, A)), log likelihood function, e.g. @(data, A) logricepdf(data, A, sigma, [false, true, false])
% data          : Measurements; (magnitude) MR values.
% maxfunArgsOut : maximum number of output arguments of fun
% explicitHessian : boolean: should explicit or compacted hessian be returned?
% logpriorfun   : if non empty: function that computes the log likelihood of the prior distribution of par.
%                               as well as the gradient and compactified hessian.
%
% OUTPUTS:
% ll : - log likelihood value
% grad: gradient of ll 
% hess: hessian of ll
%
% handles multiple data vectors simultanuously, each column of par should correspond to a column in data (and par after reshape).
%
% NOTE: update code in conjunction with voxelLLfun_sigma_m
% 
% Created by Dirk Poot, Erasmus MC, 22-3-2011
szdta = size(data);
numtr = prod(szdta(2:end));
par = reshape(par,[],numtr);
fargs = cell(1, min( max(1,nargout) , maxfunArgsOut ));
[fargs{:} ] = fun( par );
lrargs = cell(1, max(1,nargout));
[lrargs{:}] = logpdffun(data(:,:), fargs{1} );
ll = - sum( lrargs{1}(:) );
if ~isempty(logprior.fun)
    lpriorargs = cell(1, nargout);
    [lpriorargs{:}] = logprior.fun( par );
    ll = ll - sum(lpriorargs{1});
end;
if nargout>=2
    % also compute gradient:
    numim = size(data,1);
    npars = size(par,1);
    if size(fargs{2},1)==numim*numtr
        fargs{2} = reshape(fargs{2},[numim, numtr, npars]);
    end;
    lrargs{2} = reshape( lrargs{2}, numim, numtr);
    % compute gradient, see DifusionTensorLL_gradient.nb
    % evaluate dLLdA_i * dA_idpar
    
    grad = reshape(-sum( bsxfun(@times, lrargs{2}, fargs{2}), 1), numtr, npars)';
    
    if ~isempty(logprior.fun)
        grad(logprior.gradient_linindex,:) = grad(logprior.gradient_linindex,:) - lpriorargs{2};
    end;
        
    if nargout>=3
        % additionally compute hessian:
%         error('cant yet compute hessian');
        % D2_ll_term = D2_ricedist_dAdA * (dAdpar1*dAdpar2) + D_riceDist_dA * d2AdPar1dPar2
        % D2_llvoxel = sum_i D2_ricedist_dAidAi * (dAidpar1*dAidpar2) + D_riceDist_dAi * d2AidPar1dPar2
%         [I,J] = find(triu(ones(npars)));
        nI = numel(funHessian_I);
        lrargs{3} = reshape(lrargs{3}, numim, numtr);
        
        hesscoeff = hess_prep_fun_fitMRI( lrargs{3} , fargs{2}, funHessian_I, funHessian_J)';
        
%         hessLng = lrargs{3}(:,:,ones(1,nI)) .* (fargs{2}(:,:,funHessian_I).*fargs{2}(:,:,funHessian_J));
        if maxfunArgsOut>=3 && ~isempty(fargs{3})
            % Main curvature typically comes from rice distribution, so allow ignoring of function curvature (by fargs{3}==[])
            % However, when provided, include hessian:
            if size(fargs{3},1)==numim*numtr
                fargs{3} = reshape(fargs{3},[numim, numtr, npars*(npars+1)/2]);
            end;
%             hessLng = hessLng + lrargs{2}(:,:,ones(1,nI)) .* fargs{3} ;
            hesscoeff = hesscoeff + reshape( sum( lrargs{2}(:,:,ones(1,nI)) .* fargs{3} ,1), size(hesscoeff));
        end;
%         hesscoeff_orig = reshape( -sum( hessLng ,1) , numtr, nI)';
%         if any(any( abs(hesscoeff-hesscoeff_orig)> 100*eps(abs(hesscoeff_orig)+.01*max(abs(hesscoeff_orig(:))))))
%             warning('new code differs');
%         end;
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
