function [f , g , hessinfo] = LSfun2normalfun( x, LSfun , jacobmulfun)
% [f , g , hessinfo] = LSfun2normalfun( x, LSfun , jacobmulfun)
% Conversion function that wraps a function that computes the 'least squares' evaluation 
% (as requested by lsqnonlin) to a 'normal' function value, as requested by 
% fminunc (etc.). 
% INPUTS:
% LSfun(x) : function that computes the non-squared function values. (final criterium is sum(f.^2)) 
%            as well as the jacobian.
% x        : current value (passed through to LSfun)
% jacobmulfun : if empty: LSfun computes jacobian explicitly, 
%               otherwize a jacobian multiplication function with arguments: (Jinfo, Y, flag)
%               where Jinfo is the second output of LSfun.
% 
% NOTE: use LSfun2normalhessmulfun as hessian multiply function.
%
%                          fun           LSfun
% function value            f     =     sum_i f_i^2
% gradient                  grad  =     sum_i 2 * f_i * df_i/dx_j = 2 f * J
% hessian                   hess  =    2 sum_i df_i/dx_j1 df_i/dx_j2 + f_i * d2f_i/dx_j1 dx_j2
%
% We ignore the second derivative of F (since that is not computed in LS functions)
%
% Created by Dirk Poot, Erasmus MC, 22-3-2011
if isequal(x,'hessmulfun')
    f = @(jacobinfo, x) LSfun2normalhessmulfun( jacobinfo , x, jacobmulfun);
    return;
end;
args = cell(1,nargout);
[args{:}] = LSfun(x);
f = args{1}(:)'*args{1}(:);
if nargout>1
    if isempty(jacobmulfun)
        g = 2 * (f(:)' * args{2});
    else
        g = 2*jacobmulfun( args{2}, f, -1);
    end;
    hessinfo = args{2};
end;

function [Hx] = LSfun2normalhessmulfun( jacobinfo , x, jacobmulfun)
% [Hx] = LSfun2normalhessmulfun( jacobinfo , x, jacobmulfun)
% Multiplies a vector (or matrix) x with the approximate hessian of what originally is 
% a least squares function that has been wrapped with LSfun2normalfun.
% 
% 
% Created by Dirk Poot, Erasmus MC, 22-3-2011
if isempty(jacobmulfun)
    Hx = 2*((jacobinfo * x)' * jacobinfo)';
else
    Hx = 2*jacobmulfun( jacobinfo , x, 0);
end;