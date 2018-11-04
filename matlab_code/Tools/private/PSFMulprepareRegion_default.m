function [ selLrgPSF , hLLifoBlk ] = PSFMulprepareRegion_default( sellin, hLLifo )
% [ selLrgPSF ] = PSFMulprepareRegion_default( sellin, hLLifo )
%    (default for: opts.projectScaledPartPSF.prepareRegion)
% Prepares a list of voxels in which the PSF is(/might be) nonzero.
% (This list might be reduced by the masked and/or image limits)
% INPUTS:
% sellin : a row vector with linear indices (sub2ind) of the reconstruction volume
%          of the voxels that are optimized within the current block.
% hLLifo : output of PSFMulprepareGlobal_default (opts.projectScaledPartPSF.prepareGlobal)
% OUTPUTS:
% selLrgPSF : n x ndims   matrix with full indices of voxels that (migth) have non-zero PSF.
%             These points should be unique, but may be outside the reconstruction 
%             volume (they are clipped afterwards).
% 
% Created by Dirk Poot, Erasmus MC, 19-10-2011
ndim = numel(hLLifo.spatialSize);
idx = zeros(numel(sellin),ndim);
sellin = sellin.'-1;
for k=1:ndim;
    idxk = mod(sellin,hLLifo.spatialSize(k));
    sellin = (sellin - idxk)./hLLifo.spatialSize(k);
    idx(:,k) = idxk+1;
end;
selLrgPSF = reshape(bsxfun(@plus, hLLifo.projectPSFnonzeroOffsets , reshape(idx,size(idx,1),1,size(idx,2)) ) , [] , ndim);
% selLrgPSF = unique(selLrgPSF,'rows'); % make sure we dont include voxels through multiple routes.
[dum sel]= unique(selLrgPSF*cumprod([1 hLLifo.spatialSize(1:end-1)])'); % unique 'rows' is slow, so do unique on linear index (+C)
selLrgPSF = selLrgPSF(sel,:); 
hLLifoBlk = struct;
if hLLifo.projectPSFspatiallyInvariant==.5
    dst = max(idx,[],1)-min(idx,[],1);
    hLLifoBlk.basevoxidx = max( dst+1, min(hLLifo.spatialSize-dst, round(median(idx,1)))); % use median position except when the current processing  box is (too) close to the edge of the volume.
    hLLifoBlk.basevoxidx = (hLLifoBlk.basevoxidx-1)*cumprod([1 hLLifo.spatialSize(1:end-1)])' + 1;
end;
