function [ hLLifoBlk , psfdiag] = PSFMulprepareBlock_default( sellin, selLrgPSF, hLLifo, hLLifoBlk)
% [ hLLifoBlk ] = PSFMulprepareBlock_default( sellin, selLrgPSF, hLLifo)
%    (default for: opts.projectScaledPartPSF.prepareBlock)
% For each processing block, prepares hLLifoBlk, which is a structure with which PSFMulSelectedBlock_default
% and PSFMulLargeBlock_default can perform their multiplications with the psf 
%
% INPUTS:
% sellin : a row vector with linear indices (sub2ind) of the reconstruction volume
%          of the voxels that are optimized within the current block.
% selLrgPSF : a row vector with linear indices (sub2ind) of the reconstruction volume
%          of the voxels of which the 'dLLdf_full' needs to be updated. 
% hLLifo : output of PSFMulprepareGlobal_default (opts.projectScaledPartPSF.prepareGlobal)
%
% OUTPUTS:
% hLLifoBlk : Structure with which PSFMulSelectedBlock_default (opts.projectScaledPartPSF.mulsel)
%             and PSFMulLargeBlock_default (opts.projectScaledPartPSF.mulfull) can do their 
%             multiplications with the psf.
% 
% Created by Dirk Poot, Erasmus MC, 19-10-2011
if hLLifo.projectPSFspatiallyInvariant==1
    offsets = bsxfun(@minus, [sellin selLrgPSF], sellin(:));
    [dummy, loc] = ismember( offsets, hLLifo.offsetsPSF);
    tmp = hLLifo.psf(:,max(loc,1));
    tmp(:,loc==0)=0;
    hLLifoBlk.psfWithinBlock = reshape(tmp(:,1:numel(sellin)*numel(sellin)),size(tmp,1),numel(sellin),numel(sellin));
    hLLifoBlk.psfOutsideBlock = reshape(tmp(:,numel(sellin)*numel(sellin)+1:end),size(tmp,1),numel(sellin),numel(selLrgPSF));
    psfdiag = reshape(hLLifoBlk.psfWithinBlock(:, 1:numel(sellin)+1:end), size(tmp,1), 1, numel(sellin));
else
    hLLifoBlk.psfWithinBlock = cell(numel(hLLifo.projectEqualPSFgroups), 1);
    hLLifoBlk.psfOutsideBlock = cell(numel(hLLifo.projectEqualPSFgroups), 1);

    for k = 1:numel(hLLifo.projectEqualPSFgroups)
        idx = hLLifo.projectEqualPSFgroups{k}(1);
        if hLLifo.projectPSFspatiallyInvariant==.5
            tmp = estPSFpoints( hLLifo.projectGrad{idx}.logPDFhess , hLLifo.projectGrad{idx}.projectGrad, hLLifo.project(idx,:) , sellin, [sellin selLrgPSF], hLLifo.numVoxels , hLLifoBlk.basevoxidx);
            tmp(:,1:numel(sellin)) = 0.5* ( tmp(:,1:numel(sellin)) + tmp(:,1:numel(sellin))' ); % force to be symmetric.
        else
            tmp = estPSFpoints( hLLifo.projectGrad{idx}.logPDFhess , hLLifo.projectGrad{idx}.projectGrad, hLLifo.project(idx,:) , sellin, [sellin selLrgPSF], hLLifo.numVoxels );
        end;
        hLLifoBlk.psfWithinBlock{k} =  tmp(:,1:numel(sellin));
        hLLifoBlk.psfOutsideBlock{k} = tmp(:,numel(sellin)+1:end);
    end;
    tmp = cat(3,hLLifoBlk.psfWithinBlock{:});
    hLLifoBlk.psfWithinBlock  = permute(tmp(:,:,hLLifo.projectEqualPSFgroups_v),[3 1 2]);
    tmp = cat(3,hLLifoBlk.psfOutsideBlock{:});
    hLLifoBlk.psfOutsideBlock = permute(tmp(:,:,hLLifo.projectEqualPSFgroups_v), [3 1 2]);
    psfdiag = reshape(hLLifoBlk.psfWithinBlock(:, 1:numel(sellin)+1:end), size(tmp,1), 1, numel(sellin));
end;