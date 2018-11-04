function [xout , state] = HaltonSequence( npnt, base , x0 )
% [xout , final_state] = HaltonSequence( npnt, base , x0 )
% OR: [xout , final_state] = HaltonSequence( npnt, initial_state )
% Creates a halton sequence starting from an arbitrary point 0<=x0<1.
% This sequence is a pseudorandom sequence. 
% 
% INPUTS:
%   npnt : scalar integer that specifies the number of points to generate.
%   base : A postive integer base for the Halton sequence. To 'uniformly' fill a
%          multi-dimensional space use a different base for the coordinate
%          in each dimension. The bases should be relative prime and small.
%   x0   : starting point (default =0). 
%   initial_state : second output of a previous call to HaltonSequence.
%                   Continue with that sequence.
%
% OUTPUT: 
%   xout : pseudo random sequence; numel(x0)  x npnt
%   final_state : final state which can be used as input to a subsequent
%                 call to HaltonSequence to proceed with the sequence.
%
% Created by Dirk Poot, TUDelft, 22-11-2013

if nargin<3
    x0 = 0;
else
    x0 = mod(x0,1);
end;
if isstruct(base)
    state = base;
    base = state.base;
    digits = state.digits;
else
    maxdigits = max( -log(eps)./log(base) );
    mindigits = ceil( log(npnt)./min(log(base)) );
    % decompose x0 into digits
    digits = zeros(numel(x0), mindigits);
    k = 1;
    while any(x0)>0 && k<maxdigits;
        x1 = x0.*base;
        digits(:,k) = floor(x1);
        x0 = x1-digits(:,k);
        k = k+1;
    end;
end;
resid = (base-1-digits)*base.^(0:size(digits,2)-1)';
if resid<npnt
    lastdigit = ceil( log(npnt+base.^size(digits,2)-resid)/log(base) );
    digits(1,end+1:lastdigit)=0;
end;
toValue = base.^(size(digits,2)-1:-1:0)'; % should be all integers.
scale = base.^size(digits,2);
xout = zeros(numel(x0),npnt);
for k=1:npnt
    xout(:,k) = (digits * toValue)/scale; % roundoff only in division by (the integer) scale
    admsk = ones(size(xout,1),1);
    i = 1;
    while any(admsk) 
        digits(:,i) = digits(:,i)+admsk;
        admsk = digits(:,i)>=base;
        digits(admsk,i)=0;
        i=i+1;
    end;
end;
if nargout>1
    state.base = base;
    state.digits = digits;
end;
    