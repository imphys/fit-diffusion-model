function [ hLLifo ] = PSFMulprepareGlobal_default( dLLdf_full, opts )
% [ hLLifo ] = PSFMulprepareGlobal_default( dLLdf_full, opts )
%    (default for: opts.projectScaledPartPSF.prepareGlobal)
% Prepares the hLLifo structure. With this structure, for each processing block,
% PSFMulprepareRegion_default constructs a list of voxels in which the PSF is nonzero 
% and PSFMulprepareBlock_default constructs the hLLifoBlk structure. 
% hLLifo is not used for anything else.
% 
% Created by Dirk Poot, Erasmus MC, 19-10-2011
hLLifo.spatialSize = opts.spatialSize;
hLLifo.projectPSFspatiallyInvariant = opts.projectPSFspatiallyInvariant;
rep = zeros(opts.numImages,1);
for k = 1:numel(opts.projectEqualPSFgroups)
    rep(opts.projectEqualPSFgroups{k}) = k;
end;
hLLifo.projectEqualPSFgroups_v =rep;
if opts.projectPSFspatiallyInvariant==1
    hLLifo.psf = cell(1,numel(opts.projectEqualPSFgroups));
    numvox = prod(opts.spatialSize);
    steps = cumprod([1 opts.spatialSize(1:end-1)])';
    signalvoxel = floor(opts.spatialSize/2)+1;
    voxind2 = (signalvoxel-1)*steps + 1; % just the central voxel.
    flip1 = cell(1, numel(opts.spatialSize));
    flip2 = cell(1, numel(opts.spatialSize));
    for k=1 : numel(opts.spatialSize)
        baseind = max( -signalvoxel(k)+1 , -opts.spatialSize(k) + signalvoxel(k)) : min( opts.spatialSize(k) - signalvoxel(k) , signalvoxel(k)-1);
        flip1{k} = signalvoxel(k) + baseind ;
        flip2{k} = signalvoxel(k) - baseind ;
    end;
    for k = 1:numel(opts.projectEqualPSFgroups)
        idx = opts.projectEqualPSFgroups{k}(1);
        hLLifo.psf{k} = estPSFpoints( opts.projectGrad{idx}.logPDFhess , opts.projectGrad{idx}.projectGrad, opts.project(idx,:) , ':' , voxind2, numvox );
        % force symmetric (for the part in which we can do that:
        hLLifo.psf{k} = reshape(hLLifo.psf{k},opts.spatialSize);
        hLLifo.psf{k}(flip1{:}) = 0.5 * (hLLifo.psf{k}(flip1{:}) + hLLifo.psf{k}(flip2{:}));
        hLLifo.psf{k} = reshape(hLLifo.psf{k},numvox, []);
        if isempty(opts.projectPSFnonzeroOffsets{k})
            % determine PSF extend automatically:
            I = find(abs(hLLifo.psf{k})> opts.projectPSFrelativecutoff * max(abs(hLLifo.psf{k}(:))));
            Ic = cell(1,numel(opts.spatialSize));
            [Ic{:}]= ind2sub(opts.spatialSize, I);
            opts.projectPSFnonzeroOffsets{k} = bsxfun(@minus, [Ic{:}], signalvoxel);
        end;
    end;
    hLLifo.projectPSFnonzeroOffsets = unique(vertcat(opts.projectPSFnonzeroOffsets{:}),'rows');    
    hLLifo.psf = [hLLifo.psf{:}];
    hLLifo.offsetsPSF = hLLifo.projectPSFnonzeroOffsets*steps;
%     % Force psf to be symmetric (for the voxels for which we can force that):
%     inflip = all( (bsxfun(@plus, signalvoxel, -hLLifo.projectPSFnonzeroOffsets)>=1) & (bsxfun(@plus, signalvoxel-opts.spatialSize, -hLLifo.projectPSFnonzeroOffsets)<=0),2);
%     hLLifo.psf( hLLifo.offsetsPSF(inflip)+ voxind2 ,:) = 0.5*( hLLifo.psf( hLLifo.offsetsPSF(inflip)+ voxind2 , : ) + hLLifo.psf( -hLLifo.offsetsPSF(inflip)+ voxind2 , : ));
    % select psf part:
    hLLifo.psf = hLLifo.psf( hLLifo.offsetsPSF + voxind2 , hLLifo.projectEqualPSFgroups_v ).';
    
%     hLLifo.offsetsPSF =
%     size(hLLifo.psf) == numImages x numel(hLLifo.offsetsPSF)
else
    hLLifo.projectEqualPSFgroups = opts.projectEqualPSFgroups;

    hLLifo.projectGrad = opts.projectGrad;
    hLLifo.project = opts.project;
    hLLifo.numVoxels = prod(opts.spatialSize);
    hLLifo.projectPSFnonzeroOffsets = unique(vertcat(opts.projectPSFnonzeroOffsets{:}),'rows');
end;
hLLifo.projectPSFnonzeroOffsets = reshape( hLLifo.projectPSFnonzeroOffsets , [1 size(hLLifo.projectPSFnonzeroOffsets)]);
