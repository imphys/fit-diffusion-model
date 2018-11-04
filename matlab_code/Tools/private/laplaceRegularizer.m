function [f, g, h] = laplaceRegularizer(x , spacing , weights, hessianmul, cyclic)
% [f, g, h] = laplaceRegularizer(x , spacing , weights, hessianmul)
% Evaluates the log likelihood regularization value for the Laplace 
% regularizer in the spatial dimensions.
% INPUTS:
% x         : N+1 dimensional matrix for N dimensional image. First dimension = different parameters
% spacing   : spacing(i): spatial distance between samples in dimension i
% weights   : weights(i): weight of each parameter i  ( = weigth of x(i, ..) )
% hessianmul: empty/not provided: use laplaceRegularizerHessMul to multiply with the hessian.
%             0 : (default) h argument is explicit hessian. 
%             1 : h is implicit hessian, use laplaceRegularizerHessMul to multiply with hessian 
%
%             (2 : indicates that hessian * x is requested. 
%             	  should only be specified by laplaceRegularizerHessMul)
%
% Created by Dirk Poot, Erasmus MC, 22-3-2011
ismatweight = all(size(weights)==size(x,1)*[1 1]);
if nargin<4
    hessianmul = 0;
end;
if nargout>=3 && isequal(hessianmul,0)
    if nargin<5 
        [Lx L] = laplaceMulND(x,1, spacing);
    else
        [Lx L] = laplaceMulND(x,1, spacing, cyclic);
    end;
    h = L'*L;
    nx = numel(x)/size(x,1);
    if ismatweight
        repw = repmat({weights}, 1, nx );
        W = sparse(blkdiag(repw{:}));
        h = W * h;   % h is symmetric since h is scaled identity at each block.
    else
        h = spdiags(repmat(weights(:),nx,1),0,numel(x),numel(x)) * h;
    end;
else
    if nargin<5
        Lx = laplaceMulND(x,1, spacing);
    else
        Lx = laplaceMulND(x,1, spacing, cyclic);
    end;
    h=size(x);
end;
if ismatweight
    g = weights * Lx(:,:) ;
else
    g = bsxfun( @times, Lx , weights(:)) ;
end;
if nargin>=4
    if hessianmul==2
        f = g(:);
        return;
    elseif hessianmul==0 || hessianmul==1 
        % compute the hessian explicitly or implicitly.
        % but this is done above.
    elseif ~isempty(hessianmul)
        error('invalid hessianmul option');
    end;
end;
f = 0.5*(x(:)'*g(:));

