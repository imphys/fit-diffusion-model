function [Hx] = laplaceRegularizerHessMul(hessinfo, x  , spacing , weights, varargin)
% [Hx] = laplaceRegularizerHessMul(hessinfo, x  , spacing , weights)
% Multiplies a vector with the hessian of the the log likelihood 
% regularization value for the Laplace regularizer in the spatial dimensions.
% INPUTS:
% hessinfo  : third output of laplaceRegularizer
% x         : vector or multicolumn matrix with parameter vectors that should be multiplied by 
%             the hessian
% spacing   : spacing(i): spatial distance between samples in dimension i
% weights   : weights(i): weight of each parameter i  ( = weigth of x(i, ..) )
%
% Created by Dirk Poot, Erasmus MC, 22-3-2011
if numel(x)>prod(hessinfo)
    x = reshape(x,prod(hessinfo),[]); % fails if x has invalid number of elements.
    Hx = zeros(size(x));
    for k=1:size(x,2)
        Hx(:,k) = laplaceRegularizer(reshape(x(:,k),hessinfo) , spacing , weights, 2, varargin{:});
    end;
else
    Hx = laplaceRegularizer(reshape(x,hessinfo) , spacing , weights, 2, varargin{:});
end;