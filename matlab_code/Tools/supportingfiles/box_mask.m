function range = box_mask( mask, varargin )
% box = box_mask( mask , option_1, value_1, ...)
% Find a box around the nonzero elements of mask.
% I.e.: all nonzero elements of mask are contained in 
%                mask( box(1,1):box(2,1), ... , box(1,end):box(2,end) ) 
% 
% INPUTS:
%  mask : N-D logical array 
% Several options that allow processing of mask, prior to detection:
%  elementThreshold : equivalent to box_mask( mask>=elementThreshold )
%                     (but slightly faster)
%  sumThreshold     : for range in dimension i: sums mask in all dimensions ~=i
%                     and then finds elements>=sumThreshold
%  filter           : specifies a filter with which the reduced mask is filtered prior to detecting.
%                     the range. Passed as argument to separableConvN_c and thus equivalent to:
%                     box_mask( separableConvN_c( x(:) , {filter} ) )
% 
% OUTPUT:
%   box : 2 x n integer matrix, 
%        box(1,i):box(2:i) specifies the extend of the mask in dimension i.
%
% Created by Dirk Poot, Erasmus MC, 
% 31-8-2011

% default options.
opts.elementThreshold = [];
opts.sumThreshold = [];
opts.filter = []; 

if nargin <=1 
    reduce = @any;
    binarize = @logical;
else
    opts = parse_defaults_optionvaluepairs( opts, varargin{:});
    if ~isempty(opts.elementThreshold)
        reduce = @(val, dim) max(val, [],dim);
        binarize = @(x) x>=opts.elementThreshold;
    elseif ~isempty(opts.sumThreshold)
        reduce = @sum;
        binarize = @(x) x>=opts.sumThreshold;
    else
        reduce = @any;
        binarize = @logical;
    end;
    if ~isempty(opts.filter)
        binarize = @(x) binarize( separableConvN_c( x(:) , {opts.filter} ));
    end;
end;

nspatialdims = ndims(mask);
range = [1;0]*ones(1,nspatialdims);
for dim=1:nspatialdims
    % all dimensions prior to dim are reduced, so make first dimension equal to size
    % mask in dimension dim and marginalize over all other dimensions:
    dimrng = binarize( reduce( reshape(mask,size(mask,dim),[]), 2 ) );
    if nnz(dimrng)==0
        % return early if empty mask. Have to check as find cant find any if mask is empty.
        return;
    end;
    range(:,dim) = [find(dimrng,1,'first');find(dimrng,1,'last')];
    mask = reduce( mask, dim); % reduce mask to reduce work in subsequent iterations.
end;
