function [A, dA, hessA] = splitsComplex(fun, x, fields)
% [A, dA, hessA] = splitsComplex(fun, x, fields)
% Conversion function that changes the complex output of a function to a longer 
% real output, as the optimization functions require real values.
%
% Created by Dirk Poot (d.poot@erasmusmc.nl), Erasmus medical center, Rotterdam
args = cell(1,nargout);
if nargin<3 %isempty(fields)
    [args{:}] = fun(x);
else
    [args{:}] = fun(x, fields);
end;
A = [real(args{1});imag(args{1})]; % concatenate real and imaginary parts.
if nargout>1
    if size(args{1},1) ~= size(args{2},1)
        error('incorrect size of gradient, cant wrap complex -> real');
    end;
    dA = [real(args{2});imag(args{2})]; % concatenate real and imaginary parts of derivative. 
    if nargout>2
        if size(args{1},1) ~= size(args{3},1)
            error('incorrect size of hessian, cant wrap complex -> real');
        end;
        hessA = [real(args{3});imag(args{3})];
    end;
end;