function [img, imgT] = load_niiT( filename ) 
% [nii, T] = load_niiT( filename )
%  kind of inverts the store preparing of make_niiT.
%  This routine extends load_(untouch_)nii by allowing an [4 x 4] affine 
%  transformation matrix T to specify the 
%  voxel spacing, direction cosinuses and origin. 
%
%  [world-coordinate , 1 ] = [image_index , 1] * T;
%  
% OUTPUT:
%  nii : structure with nii information and image
%  T   : transformation of nii.img.
% 
% Created by Dirk Poot, Erasmsus MC
% 23-2-2011

img = load_untouch_nii( filename );
img.img = permute(img.img,[2 1 3 4:ndims(img.img)]);
if img.filetype==0
    % analize format.
    warning('Analize image loaded, origin and direction cosinusses are unknown and set to identity');
    imgT = diag([img.hdr.dime.pixdim(2:4) 1]);
    return;
end;
imgT = [img.hdr.hist.srow_x;img.hdr.hist.srow_y;img.hdr.hist.srow_z]' ;
imgT = imgT([2 1 3 4],[2 1 3]); 
imgT(:,[1 2]) = -imgT(:,[1 2]);
imgT(4,4)=1;


imgT(1:3,1:3) = diag(img.hdr.dime.pixdim([3 2 4])./sqrt( sum(imgT(1:3,1:3).^2,2) )') * imgT(1:3,1:3);
imgT(end,:) = [-1 -1 -1 1]*imgT; % compensate for nifti specifying first voxel at [0 0 0 ] and we in matlab use [1 1 1]
%% Quaternion representation:
if img.hdr.hist.qform_code>0
b = img.hdr.hist.quatern_b;
c = img.hdr.hist.quatern_c;
d = img.hdr.hist.quatern_d;
a = sqrt(1 - (b*b+c*c+d*d));
R = [ a*a+b*b-c*c-d*d   2*b*c-2*a*d       2*b*d+2*a*c     ;
      2*b*c+2*a*d       a*a+c*c-b*b-d*d   2*c*d-2*a*b     ;
      2*b*d-2*a*c       2*c*d+2*a*b       a*a+d*d-c*c-b*b  ]';
qfac = img.hdr.dime.pixdim(1);
if qfac==0
    warning('qfac should be 1 or -1');
    qfac=1;
end;
imgTq = [diag(  (img.hdr.dime.pixdim([2 3 4]).*[1 1 qfac]) * R ) [0; 0 ; 0]; 
        img.hdr.hist.qoffset_x  img.hdr.hist.qoffset_y img.hdr.hist.qoffset_z 1];
imgTq = imgTq([2 1 3 4],[2 1 3 4]);    
imgTq(:,[1 2]) = -imgTq(:,[1 2]);
imgTq(end,:) = [-1 -1 -1 1]*imgTq; % compensate for nifti specifying first voxel at [0 0 0 ] and we in matlab use [1 1 1]
end;
