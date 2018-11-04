function [HY] = compactedPSFHessianMul(hessinfo, Y)
% [HY] = compactedHessianMul(hessinfo, Y)
% Mulitiplies a vector (or matrix with few columns) with the PSF part of the hessian 
% as returned by  voxelLLfun_proj_m
%
% Created by Dirk Poot (d.poot@erasmusmc.nl), Erasmus medical center, RotterdamHY = 
if isequal(hessinfo.hesscoeff,0)
    HY = zeros(size(Y));
else
    HY = compactedHessianMul(hessinfo, Y); % handles hesscoeff part.
end;
sz = [size(hessinfo.dfdpar) 1];
ncol = size(Y,2);
tmp = reshape( sum( bsxfun(@times, hessinfo.dfdpar , reshape( Y , [1, sz(2), sz(3) , ncol])) , 2) , [sz(1) sz(3) ncol]);  % sum_i (d f_mj/ d par_ij) y_ijn
tmp2 = hessinfo.psfmul( tmp );  %  sum_j  d2 ll / d f_mj d f_ml * sum_i (d f_mj/ d par_ij) y_ijn
if ~isequal(hessinfo.project_PSFscaleBlock, 1)
    tmp2 = hessinfo.project_PSFscaleBlock .* tmp2;
end;
HY = HY + reshape( sum(bsxfun( @times, reshape(tmp2 , [sz(1) 1 sz(3) ncol]) , hessinfo.dfdpar),1) , [sz(2)*sz(3) ncol]); % tmp2 * (d f_ml / d par_kl)

