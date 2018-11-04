function [W] = compactedJacobianMul(Jinfo, Y, flag)
% [W] = compactedJacobianMul(Jinfo, Y, flag)
% Mulitiplies a vector or matrix with the jacobian of the least squares fit function.
%
% Created by Dirk Poot (d.poot@erasmusmc.nl), Erasmus medical center, Rotterdam
if size(Y,2)>1
    % split multiple right hand sides (multiple columns in Y)
    W = cell(1,size(Y,2));
    for k=1:size(Y,2)
        W{k} = compactedJacobianMul(Jinfo, Y(:,k), flag);
    end;
    W = [W{:}]; % concatenate resulting columns.
    return;
end;
% NOTE: jacobian multiplication functions should behave as:
% flag ==0 : W = J'*J*Y
% flag >0  : W = J*Y
% flag <0  : W = J'*Y
if flag==0
    W = compactedJacobianMul(Jinfo, compactedJacobianMul(Jinfo, Y, 1), -1);
elseif flag >0
    % Jinfo = (nMRI x numel(index_range) x size(par,1))
    % J(i,j) = d( A_i ) /d( par_j )
    % Y = (size(par,1) x numel(index_range))
    % W = ( nMRI x numel(index_range)) 
    Y = reshape(Y,size(Jinfo,3),[])';
    Y = reshape(Y, 1, [], size(Jinfo,3));
    W = sum(bsxfun(@times, Jinfo, Y ), 3);
    W = W(:);
else
    % Y = ( nMRI x numel(index_range)) 
    % W = (size(par,1) x numel(index_range))
    Y = reshape(Y,size(Jinfo,1) ,[]);
    W = permute(sum(bsxfun(@times, Y, Jinfo),1) , [3 2 1]);
    W = W(:);
end;