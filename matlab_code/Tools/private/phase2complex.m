function [A, dA, hessA] = phase2complex(fun, x , fields, funHessian_I, funHessian_J)
% [A, dA, hessA] = phase2complex(fun, x, fields, funHessian_I, funHessian_J)
% Conversion function that changes the output of a function that returns a 
% phase image into a complex image in which the fitting is easier 
% (doesn't suffer from phase wraps)
%
% Created by Dirk Poot (d.poot@erasmusmc.nl), Erasmus medical center, Rotterdam
args = cell(1,nargout);
if isempty(fields)
    [args{:}] = fun(x); % call function 
else
    [args{:}] = fun(x, fields); % call function 
end;
A = exp(1i*args{1}); % convert function value
if nargout>1
    % f(A(x)) = e^(i*A(x))
    % df/dx = df/dA dA/dx = f(A(x)) * 1i* dA/dx
    dA = bsxfun(@times, 1i*A, args{2}); % convert derivative
    if nargout>2
        % d2f/dx1dx2 = d2f/d2A dA/dx1 dA/dx2 + df/dA d2A/dx1dx2
        %            = f(A(x)) (-1* dA/dx1 * dA/dx2 + 1i * d2A/dx1dx2)
        if isempty(args{3})
            hessA = bsxfun(@times, A, -args{2}(:,funHessian_I).*args{2}(:,funHessian_J));
        else
            hessA = bsxfun(@times, A, -args{2}(:,funHessian_I).*args{2}(:,funHessian_J) + 1i*args{3});
        end;
    end;
end;