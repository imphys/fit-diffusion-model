function [ psf , grad ] = estPSFpoints( d2L, grad, project , voxind1, voxind2, numvox , basevoxidx)
% Estimates the PSF of several voxels. This is a general, but very inefficient method (although in general we cant do any better).
% A custom replacement should not project grad, but use project specific specializations. 
%
% Created by Dirk Poot, Erasmus MC, 28-10-2011
swap = (numel(voxind1)<numel(voxind2) ) || strcmp(voxind2,':');
if swap
    [voxind1, voxind2] = deal(voxind2, voxind1);
end;
if strcmp(voxind1,':')
    nvox1 = numvox;
else
    nvox1 = numel(voxind1);
end;
if nargin>=7
    % use basevoxidx to derive all psf's from 
    x = zeros(numvox,1);x(basevoxidx)=1;
    [Gx , grad] = project1( x, grad , project, 2);
    [GHGx , grad] = project1( Gx(:) .* d2L(:), grad , project, 3);
    if strcmp(voxind1,':')
        voxind1 = 1:numvox;
    end;
    psf = GHGx( max(1,min(numvox, bsxfun(@minus, voxind1(:), voxind2(:)'-basevoxidx))) ); % index into GHGx is difference betweeen voxind1- voxind2 + basevoxidx
    % don't care about internal wrapping around, as the PSF should typically not larger than the image size - blocksize,
    % so no actual wrapping around happens.
else
    psf = zeros( nvox1, numel(voxind2) );
    for k = 1:numel(voxind2)
        x = zeros(numvox,1);x(voxind2(k))=1;
        [Gx , grad] = project1( x, grad ,  project, 2);
        if iscell(Gx)
            Gxd2L = cell(size(Gx));
            for i = 1:numel(Gx)
                Gxd2L{i} = Gx{i}(:) .* d2L{i}(:);
            end;
        else
            Gxd2L = Gx(:) .* d2L(:);
        end;
        [GHGx , grad] = project1( Gxd2L, grad , project, 3);
        psf(:,k) = GHGx(voxind1);
    end;
end;
if swap 
    psf = psf.';
end;