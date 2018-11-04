function [ll, grad, hessinfo] = voxelLLfun_proj_m( par, fun, fpred_in, dLLsel, psfmul, psfdiag, maxfunArgsOut, explicitHessian, logprior, project_PSFscaleBlock, funHessian_I, funHessian_J)
% [ll, grad, hess] = voxelLLfun_m( par, fun, logpdffun, data , maxfunArgsOut, explicitHessian, logpriorfun, project_PSFscaleBlock)
% computes the likelihood of the MR data with likelihood function logpdffun.
% INPUTS:
% par           : a vector/matrix with the parameters for fun (each column of par specifies a voxel)
% fun(par)      : a function predicting the magnitudes of the MR images for each column of par
%                 (Each column specifies the parameters of one voxel)
% fpred_in      : fun( par_init ), where dLLsel and psfmul are constructed with par_init.
% dLLsel        : numimg x numtraces gradient of the log likelihood
%                  dLLsel(i,j) = dLL / d f_ij, where f_ij is i,j element of fun( par ) 
% psfmul        : function that computes  [Hx] = psfmul( x ) 
%                 H(i,j) = d2 LL / d f_i d f_j
%                 where f_i is element i of first output of fun( par )
% maxfunArgsOut : maximum number of output arguments of fun
% explicitHessian : boolean: should explicit or compacted hessian be returned?
% logpriorfun   : if non empty: function that computes the log likelihood of the prior distribution of par.
%                               as well as the gradient and compactified hessian.
% project_PSFscaleBlock : scale factor for ouput of psfmul. 
%
% OUTPUTS:
% ll : - log likelihood value
% grad: gradient of ll 
% hess: hessian of ll
%
% handles multiple data vectors simultanuously, each column of par should correspond to a column in data (and par after reshape).
%
% NOTE: update code in conjunction with voxelLLfun_sigma_m and voxelLLfun_m
% 
% Created by Dirk Poot, Erasmus MC, 22-3-2011
numtr = size(dLLsel,2);
par = reshape(par,[],numtr);
fargs = cell(1, min( max(1,nargout) , maxfunArgsOut ));
[fargs{:} ] = fun( par );
deltaf = fargs{1}-fpred_in; % size = numim x numtr
Hdeltaf = psfmul( deltaf );
if ~isequal(project_PSFscaleBlock,1)
    Hdeltaf  = Hdeltaf  .* project_PSFscaleBlock;
end;
ll =  reshape(deltaf,1,[])*dLLsel(:) + 0.5 * (deltaf(:)'*Hdeltaf(:)); % approximate -log likelihood function when projecting.
if ~isempty(logprior.fun)
    lpriorargs = cell(1, nargout);
    [lpriorargs{:}] = logprior.fun( par );
    ll = ll - sum(lpriorargs{1});
end;
if nargout>=2
    % also compute gradient:
    % d ll / d par_ij = sum_m d ll / d f_mj * d f_mj / d par_ij
    numim = size(dLLsel,1);
    npars = size(par,1);
    if size(fargs{2},1)==numim*numtr
        fargs{2} = reshape(fargs{2},[numim, numtr, npars]);
    end;
    
    % compute gradient.
    % evaluate dLLdA_i * dA_idpar
    dLLseldf = dLLsel +  reshape(Hdeltaf,size(dLLsel));
    grad = reshape( sum( bsxfun(@times, dLLseldf, fargs{2}), 1), numtr, npars)';
    
    if ~isempty(logprior.fun)
        grad(logprior.gradient_linindex,:) = grad(logprior.gradient_linindex,:) - lpriorargs{2};
    end;
        
    if nargout>=3
        % additionally compute hessian:
        % d2 ll / d par_ij d par_kl = sum_mn  d2 ll / d f_mj d f_nl * (d f_mj/ d par_ij) * (d f_nl / d par_kl)
        %                             + sum_m  d ll / d f_mj * (d2 f_mj / d par_ij d par_kl )            
        %                           =   sum_mn delta_mn delta_jl d2 ll / d f_mj d f_ml * (d f_mj/ d par_ij) * (d f_ml / d par_kl)
        %                             + delta_jl  sum_m  d ll / d f_mj        * (d2 f_mj / d par_ij d par_kj )
        %  => always with delta_jl if no psf =>
        % d2 ll / d par_ij d par_kj =   sum_m  d2 ll / d f_mj d f_mj * (d f_mj/ d par_ij) * (d f_mj / d par_kj)
        %                             + sum_m   d ll / d f_mj        * (d2 f_mj / d par_ij d par_kj )
%         [I,J] = find(triu(ones(npars)));
        nI = numel(funHessian_I);
        hesscoeff = 0;
        if maxfunArgsOut>=3 && ~isempty(fargs{3})
            % Main curvature typically comes from distribution, so allow ignoring of function curvature (by fargs{3}==[])
            % However, when provided, include hessian:
            if size(fargs{3},1)==numim*numtr
                fargs{3} = reshape(fargs{3},[numim, numtr, npars*(npars+1)/2]);
            end;
            hessLng = bsxfun(@times, dLLsel , fargs{3} );
            hesscoeff = reshape( sum( hessLng ,1) , numtr, nI)';  % sum_m d ll / d f_mj * (d2 f_mj / d par_ij d par_kj )          
        end;
        if ~isempty(logprior.fun)
            hesscoeff(logprior.hessian_write_funhessindex,:) = hesscoeff(logprior.hessian_write_funhessindex,:) - lpriorargs{3};
        end;
        hessinfo.hesscoeff = hesscoeff;
        hessinfo.psfmul = psfmul;
        hessinfo.dfdpar = permute( fargs{2} ,[1 3 2]);
        hessinfo.hessMulFun = @compactedPSFHessianMul;
        hessinfo.makePreconditioner = @compactedHessianMakePreconditioner;      
        hessinfo.project_PSFscaleBlock =  project_PSFscaleBlock;
        hessinfo.psfdiag = psfdiag;
        if explicitHessian
            error('explicit hessian should not be computed in voxelLLfun_proj_m');
        end;
    end;
end;
if ~isfinite(ll)
    ll = inf;
end;

