function [varargout] = clampargCallFun( fun, par, minval, maxval )
% [f, grad, hess] = clampargCallFun( fun, par, minval, maxval )
% clamps the parameters par between minval and maxval and calls fun:
%   [f, grad, hess] = fun( max(minval, min(maxval, par)) ) ;
%   grad(clamped_values) = 0;
% INPUTS:
%  fun : a function taking 1 argument
%  par : a numeric array passed to fun
%  minval: empty, scalar or size(par) numeric array with minimum value for par
%  maxval: empty, scalar or size(par) numeric array with maximum value for par
% OUTPUTS:
%  f    : function value evaluated at the clamped par
%  grad : gradient as returned by the fun; except for the clamped value at which it is set to zero.
%  hess : hessian as returned by the fun.
%
% Created by Dirk Poot, Erasmus MC.
% Extracted from fit_MRI at 27-1-2012

clamped = [];
if nargin>=3 && ~isempty(minval)
    clamped = par<minval;
    if numel(minval)>1
        % if non scalar minval: assume same size as par. Subselect the clamped values.
        minval = minval(clamped);
    end;
    par(clamped) = minval;
end;
if nargin>=4 && ~isempty(maxval)
    clampedmax = par>maxval;
    if numel(maxval)>1
        % if non scalar minval: assume same size as par. Subselect the clamped values.
        maxval = maxval(clampedmax);
    end;
    par(clampedmax) = maxval;
    if ~isempty(clamped)
        clamped = clamped | clampedmax;
    else
        clamped = clampedmax;
    end;
end;
varargout = cell(1,nargout);
[varargout{:}] = fun( par ); % Call function.

if nargout>=2 && any(clamped(:))
    varargout{2}(clamped) = 0; % gradient of clamped values is zero.
    % dont treat hessian, as the format is different
end;