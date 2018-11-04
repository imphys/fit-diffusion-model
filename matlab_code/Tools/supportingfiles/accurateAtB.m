function [dotprod_hi, dotprod_lo] = accurateAtB(A, B, maxabsA, maxabsB)
% [dotprod_hi, dotprod_lo] = accurateAtB(A, B)
% computes:
%   dotprod = A'*B 
% in a accurate way (more accurate than just typing A'*B in MATLAB).
% Typically the roundoff error due to the standard MATLAB matrix multiplication
% will be of order:  eps(dotprod)*sqrt(size(A,1))
% but in worst case it might be  eps(abs(A)'*abs(B)) * size(A,1) /16
% Especially for single precision arguments with large size(A,1), this
% might be significant.
%
% For double precision arguments:
% The first output argument is the most significant and will almost
% always be the correctly rounded result of A'*B (when the product would have been
% computed with infinite precision). The second output is an approximation
% of the remainder that is typically only accurate for some 10-20 bits.
% (Almost) always:  abs(dotprod_lo) < eps(dotprod_hi)
%
% For (GPU)single arguments:
% (Implicitly) the arguments are converted to double and then they are
% multiplied, returning (GPU)double precision results. The second output is
% always zero.
% 
% Created by Dirk Poot, University of Antwerp,
% 4-11-2010.

szA = size(A);
szB = size(A);
if numel(szA)>2 || numel(szB)>2 || szA(1)~=szB(1)
    error('A and B should be matrices with matching size of first dimension.');
end;
% first handle special classes:
if isa(A,'single') || isa(B,'single');
    dotprod_hi = double(A')*double(B);
    dotprod_lo = 0;
    return;
end;
if isa(A,'GPUsingle') && isa(B,'GPUsingle')
    dotprod_hi = zeros([szA(2),szB(2),60],GPUdouble);
    dotprod_lo = 0;
    GPUdotproddouble(A,B,dotprod_hi);
    setSize(dotprod_hi,[szA(2),szB(2)]);
    return;
end;
if ~isa(A,'double') || ~isa(B,'double')
    warning('AccurateAtB:ConvertClass',['A and/or B are converted from "' class(A) '", "' class(B) '" to "double"; I don''t know the precision of the original type, some precision might be lost in this conversion.']);
    A = double(A);
    B = double(B);
end;
if isreal(A) && isreal(B)
    % more accurate and faster version of the general code below:
    [dotprod_hi,dotprod_lo] = accurate_dotprod_c(A,B);
    return;
end;
%% normal double precision accurate multiplication:
if nargin<3
    maxabsA = sum(A.*conj(A),1);
else
    maxabsA = maxabsA.^2 * szA(1); 
    if ~isequal(size(maxabsA),[1 szA(2)]) || any(maxabsA<0)
        error('incorrect maxAbsA specified');
    end;
end;
if nargin<4
    maxabsB = sum(B.*conj(B));
else
    maxabsB = maxabsB.^2 * szB(1);
    if ~isequal(size(maxabsB),[1 szB(2)]) || any(maxabsB<0)
        error('incorrect maxAbsB specified');
    end;
end;
% if each column of A & B is multiplied by 1/sqrt(maxabsA) and 1/sqrt(maxabsB) 
% then for all elements abs(dotprod)<=1
% in order to guarantee that the scaling does not introduce roundoff errors
% we scale with the first power of 2 larger than sqrt(maxabsA).
% Additionally we scale such that abs(dotprod)<=2^53 as this is the
% maximum integer exactly representable (with double precision numbers) 
% 
Ascpwr2 = ceil( 0.5*log2(  maxabsA +realmin ) ) - 26; %26
Bscpwr2 = ceil( 0.5*log2(  maxabsB +realmin ) ) - 27; %27

A_hi = bsxfun(@times, round( bsxfun(@times, A, pow2(-Ascpwr2) ) ), pow2(Ascpwr2));
A_lo = A- A_hi;
clear A

B_hi = bsxfun(@times, round( bsxfun(@times, B, pow2(-Bscpwr2) ) ), pow2(Bscpwr2));
B_lo = B- B_hi;
clear B

dotprod_hi1 = A_hi' * B_hi; % this one is exact (i.e. no roundoff errors anywhere)
dotprod_lo1 = A_lo' * B_hi; 
dotprod_lo2 = A_hi' * B_lo; 
dotprod_lo3 = A_lo' * B_lo;
% next lines: try to make sure that dotprod_hi is correctly rounded
%             and dotprod_lo preserves as much precision as possible.
dotprod_hi = dotprod_hi1 + (dotprod_lo1 + dotprod_lo2 + dotprod_lo3);
if nargout>1
    dotprod_lo = (dotprod_hi1 - dotprod_hi) + dotprod_lo1 + dotprod_lo2 + dotprod_lo3;
end;