function [val, dval, hessinfo] = regularizedLLFun(par, voxLLfun, regularizer, outerpar, outerparsel)
% [val, dval, hess] = regularizedLLFun(par, voxelLLfun, regularizerinfo, outerpar, outerparsel)
% Computes the combined log likelihood of voxelLLfun_m and the regularization function
% specified in opts.
%
% Created by Dirk Poot, Erasmus MC, 22-3-2011


args1 = cell(1,nargout);
[args1{:}] = voxLLfun(par(:,:));
if nargin<5
    outerpar = par;
else
    outerpar(:,outerparsel) = par(:,:);
end;
args2 = cell(1,nargout);
[args2{:}] = regularizer.fun( outerpar );
val = args1{1} + args2{1};
if nargout>=2
    if nargin>=5
        args2{2} = reshape(args2{2}(:, outerparsel), size(par));
    else
        args2{2} = reshape(args2{2}, size(par));
    end;
    dval = reshape(args1{2},size(par)) + args2{2};
    if nargout>=3
        if nargin>=5 
            selall = bsxfun(@plus, (1:size(par,1))', (outerparsel(:)'-1)*size(par,1));
            if isempty(regularizer.hessMulFun)
                args2{3} = args2{3}(selall,selall); 
            end;
        end;
        if regularizer.explicitHessian
            if ~isempty(regularizer.hessMulFun)
                if nargin<5
                    selall = 1:numel(outerpar);
                end;
                tst = full(sparse(selall(:), (1:numel(selall))',ones(numel(selall),1), numel(outerpar),numel(selall)));
                tst = regularizer.hessMulFun( args2{3}, tst );
                args2{3} = tst(selall,:);
            end;
            hessinfo = args1{3} + args2{3};
        else
            hessinfo.Hfun = args1{3};
            hessinfo.Hregularizer = args2{3};
            if nargin>=5 && ~isempty(regularizer.hessMulFun)
                hessinfo.parmask = selall;
                hessinfo.outerparsz = size(outerpar);
            else
                hessinfo.parmask = [];
            end;
            hessinfo.hessMulFun = make2arg_anonfun( @regularizedHessMulFun, regularizer.hessMulFun);
            hessinfo.makePreconditioner = @compactedHessianMakePreconditioner;              
        end;        
    end;
end;
