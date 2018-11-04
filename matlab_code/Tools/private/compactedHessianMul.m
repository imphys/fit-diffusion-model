function [HY] = compactedHessianMul(hessinfo, Y)
% [HY] = compactedHessianMul(hessinfo, Y)
% Mulitiplies a vector (or matrix with few columns) with a 'compactified' hessian 
% as returned by voxelLLfun_m
%
% Created by Dirk Poot (d.poot@erasmusmc.nl), Erasmus medical center, Rotterdam
if size(Y,2)>1
    % sometimes multiple right hand sides are given, but we cant handle that efficiently.
    % So split and return combined results:
    for k=1:size(Y,2)
        Y(:,k) = compactedHessianMul(hessinfo, Y(:,k));
    end;
    HY = Y;
    return;
end;
hessinfo = hessinfo.hesscoeff;
szHinf = size(hessinfo);
Y = reshape(Y,[],szHinf(2:end));
if 0 
    % fastest method with pure MATLAB code:
    HY = zeros(size(Y));
    rdindx = 1;
    for k1 = 1:size(Y,1)
        for k2 = 1:k1-1
            HY(k1,:) = HY(k1,:) + hessinfo(rdindx,:).*Y(k2,:);
            HY(k2,:) = HY(k2,:) + hessinfo(rdindx,:).*Y(k1,:);
            rdindx = rdindx+1;
        end;
        HY(k1,:) = HY(k1,:) + hessinfo(rdindx,:).*Y(k1,:);
        rdindx = rdindx+1;
    end;
    if rdindx~=size(hessinfo,1)+1
        error('wrong size of hessinfo or Y');
    end;
elseif 1
    % mex routine for extra speed
    HY = compactedHessianMultiply_c( hessinfo, Y );
elseif 0
    %% alternative method.
    % Strangely enough this appears to be slower than the first method above.
    m = triu(true(size(Y,1)));
    mi = double(m);mi( m ) = 1:size(hessinfo,1);mit = mi'; mi(~m) = mit(~m);
    HYc = cell(size(Y,1),1);
    for k=1:size(Y,1)
        HYc{k} = sum( hessinfo(mi(:,k),:).*Y(:,:),1);
    end;
    HY2 = vertcat(HYc{:});
%     disp(['max abs difference: ' num2str( max(abs(HY2(:)-HY(:))) ) ])
end;
HY = HY(:);

