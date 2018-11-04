function [predicted, grad, hess] = predictAllImagesfun( theta, opts , gradmuldir)
% [predicted, grad, hess] = predictAllImagesfun( theta, opts [, gradmuldir])
% Function to predict all images with fun.
% INPUTS:
%   theta   : ND matrix with all parameters
%   opts    : options structure of fit_MRI
%   gradmuldir : optional 'vector' with same number of elements as theta, with which the gradient is multiplied. 
%                (mainly for line search in which you are only interested in the gradient along the line)
% OUTPUTS:
%   predicted: nimg x size(theta,1) x .... x size(theta,end) with result of fun(theta) for each voxel.
%   grad     : if 'gradmuldir' : size(predicted) with the gradient along gradmuldir in each predicted value.
%                   else       : nimg x num_voxels x size(theta,1), full gradient in each voxel
%   hess     : if 'gradmuldir' : size(predicted) with the hessian along gradmuldir in each predicted value.
%                   else       : nimg x num_voxels x (size(theta,1)*(size(theta,1)+1)/2), full hessian in each voxel.
%             NOTE: only request hess when opts.maxfunArgsOut>=3  (error
%                   otherwise)
%
% Created by Dirk Poot, Erasmus MC, 18-10-2011
fieldarg ={};
if opts.numPDFoptpar>0
    thetar = theta(1:end-opts.numPDFoptpar,:);
else
    thetar = theta;
end;
ntraces = prod(opts.spatialSize);

grhessout = cell(1,nargout-1);
maxninblock = opts.maxTracesInBlock(nargout);
% Test if we can do one call, or need a loop:
if  maxninblock >= ntraces 
    % one call suffices; dont waste memory by first initializing variables.
    if ~isempty(opts.fields)
        fieldarg = {opts.fields(:,:)};
    end;
    [predicted , grhessout{:}]= opts.function(thetar(:,:), fieldarg{:} );
    predicted = reshape(predicted, [opts.numImages opts.spatialSize]);
    if nargout>=2
        if isequal(size(grhessout{1}),[opts.numImages*ntraces size(thetar,1)])
            grhessout{1} = reshape(grhessout{1}, [opts.numImages ntraces size(thetar,1)]);
            if (numel(grhessout)>=2) && isequal(size(grhessout{2}),[opts.numImages*ntraces numel(opts.funHessian_I)])
                grhessout{2} = reshape(grhessout{2}, [opts.numImages ntraces numel(opts.funHessian_I)]);
            end;
        end;
        if nargin>=3
            grhessout{1} = reshape(grhessout{1}, [opts.numImages prod(opts.spatialSize) size(thetar,1)]);
            grad = reshape( sum( bsxfun(@times, grhessout{1}, permute(reshape(gradmuldir, size(theta,1),[]),[3 2 1])) ,3) , size(predicted));
            if nargout>=3
                tmp = permute(reshape(gradmuldir, size(theta,1),[]),[3 2 1]);
                cnt = 1+(opts.funHessian_I(:)~=opts.funHessian_J(:));
                tmp2 = bsxfun(@times, grhessout{1}, tmp(:,:,opts.funHessian_I).*tmp(:,:,opts.funHessian_J));
                hess = reshape( reshape( tmp2 ,[], size(tmp2,3)) * cnt  , size(predicted));
            end;
        else
            grad = grhessout{1};
            if nargout>=3
                hess = grhessout{2};
            end;
        end;
    end;
else
    predicted = zeros([opts.numImages opts.spatialSize]);
    if nargout>=2
        if nargin>=3
            gradmuldir = reshape(gradmuldir, size(theta,1),[]);
            grad = predicted;
            if nargout>=3
                cnt = 1+(opts.funHessian_I(:)~=opts.funHessian_J(:));
                hess = predicted;
            end;
        else
            grad = zeros([opts.numImages prod(opts.spatialSize) size(theta,1)]);
            if nargout>=3
                hess = zeros([opts.numImages prod(opts.spatialSize) size(theta,1)*(size(theta,1)+1)/2]);
            end;
        end;
    end;    
    for k = 1 : maxninblock : ntraces
        ed = min(k+ maxninblock -1,ntraces);
        if ~isempty(opts.fields)
            fieldarg = {opts.fields(:,k:ed)};
        end;
        [predicted(:,k:ed), grhessout{:}] = opts.function(thetar(:,k:ed), fieldarg{:});
        if nargout>=2
            if isequal(size(grhessout{1}),[opts.numImages*(ed-k+1) size(thetar,1)])
                % allow concatenation of predicted values in gradient and hessian (so they are 2D)
                grhessout{1} = reshape(grhessout{1}, [opts.numImages (ed-k+1) size(thetar,1)]);
                if (nargout>=3) && isequal(size(grhessout{2}),[opts.numImages*(ed-k+1) numel(opts.funHessian_I)])
                    grhessout{2} = reshape(grhessout{2}, [opts.numImages (ed-k+1) numel(opts.funHessian_I)]);
                end;
            end;
            if nargin>=3
                tmp = permute(gradmuldir(:,k:ed),[3 2 1]);
                grad(:,k:ed) = sum( bsxfun(@times, grhessout{1}, tmp) ,3) ;
                if nargout>=3
                    hess(:,k:ed) = reshape( reshape( bsxfun(@times, grhessout{1}, tmp(:,:,opts.funHessian_I).*tmp(:,:,opts.funHessian_J)) , [], numel(cnt))*cnt, [size(hess,1), ed-k+1]) ;
                end;
            else
                grad(:,k:ed,:) = grhessout{1};
                if nargout>=3
                    hess(:,k:ed,:) = grhessout{2};
                end;
            end;
        end;
    end;
end;

