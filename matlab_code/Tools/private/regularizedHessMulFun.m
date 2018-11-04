function [HY] = regularizedHessMulFun(hessinfo, Y, spatialRegularizerHmul)
% [HY] = regularizedHessMulFun(hessinfo, Y, spatialRegularizerHmul)
% Mulitiplies a vector or matrix with the hessian of the fit function as well as with 
% the hessian of the regularization function. 
%
% Created by Dirk Poot (d.poot@erasmusmc.nl), Erasmus medical center, Rotterdam
if isempty(spatialRegularizerHmul) 
    HY2 = hessinfo.Hregularizer * Y;
else
    if ~isempty(hessinfo.parmask)
        Y2 = zeros(prod(hessinfo.outerparsz),size(Y,2));
        Y2(hessinfo.parmask,:) = Y;
    else
        Y2 = Y;
    end;
    HY2 = spatialRegularizerHmul( hessinfo.Hregularizer , Y2 );
    if ~isempty(hessinfo.parmask)
        HY2 = HY2(hessinfo.parmask,:);
    end;
end;
HY = hessinfo.Hfun.hessMulFun(hessinfo.Hfun, Y);
HY = HY + HY2;