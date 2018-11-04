function [LL, gL, H] = fullThetaCostfun( theta, opts )
% [LL, gL, H] = fullThetaCostfun( theta, opts )
% Evaluates the full cost function with the voxel-wise model parameters as parameters
% to which the derivative is taken.
%
% theta : [nparameters  spatialsize]
% opts  : fit_MRI option structure.
%
% LL : scalar log likelihood at current theta.
% gL : gradient of LL w.r.t theta
% H  : hessian multiplication info.
%
% Created by Dirk Poot, Erasmus MC, 2-2-2012

data = cell(1,min(nargout, opts.maxfunArgsOut));
[data{:}] = predictAllImagesfun( theta, opts );
if nargout==1
    LLf = projectAllAndLL(data{1}, opts);
elseif nargout==2
    [LLf ,projgradifo, LLgrads] = projectAllAndLL(data{1}, opts , 1);
else
    [LLf ,projgradifo, LLgrads] = projectAllAndLL(data{1}, opts , 2);
end;
LL = sum(LLf(:));
if opts.doRegularize
    regul = cell(1,nargout);
    [regul{:}] = opts.spatialRegularizer{1}(theta);
    LL = LL + regul{1};
end;
if nargout>=2
    data{1} = []; % free memory
    gL = opts.function_jacMulAdj( data{2} , LLgrads(:,:) ) ;
    if opts.doRegularize
        gL = gL + regul{2}(:,:);
    end;
    if nargout>=3
        % store partial function values, for output.
        if opts.doComputeDerivativeRegularizationScale
            H.dRegularizationdTheta = regul{2}(:);
        end;    
        regul{2}=[]; % free memory
        H.projGradHess = projgradifo;
        H.gfun = data{2};
        if numel(data)>=3
            if numel(data{3})>1e8
                % construct tmp in parts (to use less memory):
                szd3 = size(data{3});
                tmp = zeros(1,szd3(2),szd3(3) );
                LLgrads = LLgrads(:,:);
                step = max(10,floor(1e5/(szd3(1)*szd3(3)))); % treat about 800kB of data{3} per iteration.
                for k = 1 : step : szd3(2)
                    ed = min( k+step-1, szd3(2) );
                    tmp(:,k:ed,:) = sum(bsxfun(@times, data{3}(:,k:ed,:), LLgrads(:,k:ed)),1);
                end;
            else
                tmp = sum(bsxfun(@times, data{3}, LLgrads(:,:)),1);
            end;
            mat = zeros(size(theta,1));
            mat(opts.funHessian_I + (opts.funHessian_J-1) * size(theta,1)) = 1:numel(opts.funHessian_I);
            mat(opts.funHessian_J + (opts.funHessian_I-1) * size(theta,1)) = 1:numel(opts.funHessian_I);
            tmp = tmp(:,:,max(1,mat(:)));
            if any(mat(:)==0)
                tmp(:,:,mat(:)==0)=0;
            end;
            H.HPfun = permute( reshape(tmp, [size(data{3},2), size(mat)]),[2 1 3]); 
        else
            H.HPfun = [];
        end;
        
        H.Hregul = [];
        if opts.doRegularize    
            H.Hregul = regul{3};
            H.sizeTheta = size(theta);
            H.LLregularization = regul{1};
        else
            H.LLregularization = 0;
        end;
        H.LLfun            = LLf;
        
%         hessinfo.hessMulFun = @fullThetaCostfun_HessMul;
    end;
end;