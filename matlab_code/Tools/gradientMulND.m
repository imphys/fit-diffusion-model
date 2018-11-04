function [image, L] = gradientMulND( image, lambda, voxelspacing)
% [ f ] = gradientMulND( image, lambda [, voxelspacing])
% Let
%   G = cell(1,ndims(image));
%   [G{:}] = gradient(image , voxelspacing(2), voxelspacing(1), voxelspacing(3),...,voxelspacing(ndims(image)) );
% This function computes f such that 
%   f(:)' * image(:) = lambda * sum_k sum( G{k}(:).^2 )
% 
% Thus this function multiplies with the gradient operator and it's adjoint.
% It is much faster and more memory efficient than an equivalent MATLAB implementation.
%
% INPUTS:
%  image  : N-dimensional image
%  lambda : scalar with which the result is scaled.
%  voxelspacing :  default = ones(1,ndims(image))
%                  Specifies the voxel spacing with which the gradient is computed.
% OUTPUTS:
%  f      : the image multiplied by the gradient and its adjoint.
%           Note the gradient is not computed explicitly.
% 
% Created by Dirk Poot
% Erasmus MC 28-6-2011

if nargin>=3 && ~isempty( voxelspacing )
    rvoxelspacing = {1./voxelspacing};
else
    rvoxelspacing = {};
end;

% Call c++ subroutine that does the hard work:
image = gradientMulND_c(image, lambda, rvoxelspacing{:});

if 0
%% Compute gradient operator in 1D & 2D:
    for l1 = 2:7
        for l2 = 1:7
            G = zeros(l1*l2*(1+(l2>1)),l1*l2);
            Gc= zeros(l1*l2);
            for k=1:l1*l2
                v = zeros(l1,l2);v(k)=1;
                if l2>1
                    [g1,g2] = gradient(v);
                else
                    [g1] = gradient(v);
                    g2 = zeros(0,1);
                end;
                G(:,k) = [g1(:);g2(:)];
                g3 = gradientMulND_c(v,1);
                Gc(:,k) = g3(:);
            end;
            [l1 l2]
            err = Gc-G'*G;
            if any(err(:))
                err
            end;
        end;
    end;
    
end;