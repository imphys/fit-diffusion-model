function [nii] = make_niiT(img, T, varargin)
% [nii] = make_niiT(img, T, [datatype], [description] )
%  Make NIfTI structure specified by an N-D matrix. This routine extends 
%  make_nii by allowing an [4 x 4] affine transformation matrix T to specify the 
%  voxel spacing, direction cosinuses and origin. 
%  Other arguments are passed to make_nii.
%  Subsequently call
%    save_nii( nii, [filename '.nii']);
%  to store the image to disk
%
%  Definition of T:
%    [world-coordinate , 1 ] = [image_index , 1] * T;
%  
%
% Created by Dirk Poot, Erasmsus MC
% 23-2-2011

if size(T,1)==3
    T(4,1:4) = [T(3,[1 2]) 0 T(3,3)];
    T(3,:) = [0 0 1 0];
end;
T(end,:) = sum(T,1); % compensate for nifti specifying data starts at [0 0 0] and we in matlab start at [1 1 1]
T = T([2 1 3:end],[2 1 3:end])';
T(1:2,:) = -T(1:2,:);

voxsp = sqrt(sum(T(:,1:end-1).^2,1));
img = permute(img,[2 1 3:ndims(img)]);
nii = make_nii( img, voxsp, zeros(size(voxsp)) , varargin{:});
nii.hdr.hist.srow_x = T(1,:);%./[voxsp 1];
nii.hdr.hist.srow_y = T(2,:);%./[voxsp 1];
nii.hdr.hist.srow_z = T(3,:);%./[voxsp 1];
nii.hdr.hist.sform_code = 1;
