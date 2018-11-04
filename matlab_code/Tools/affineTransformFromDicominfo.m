function [transform] = affineTransformFromDicominfo(info)
% [transform] = affineTransformFromDicominfo(info)
% [transform] = affineTransformFromDicominfo(filename)
% Extracts the affine transformation parameters from the dicom info
% structure. 
%
% INPUTS:
% info  : dicominfo structure of the dicom file from which you want the
%         affine transformation
% filename: string with the name of a dicom file. 
% 
% OUTPUT:
% transform: 4 x 4 affine transformation matrix.
%            When img is an image of which the dicom information is provided:
%            img(x1,x2,x3) has physical location [x1 x2 x3 1]*transform
%               
% Created by Dirk Poot, Erasmus MC
% 28-5-2010

if ischar(info)
    info = dicominfo(info);
end;
if isfield(info, 'PixelSpacing')
    pxsp = info.PixelSpacing;
else
    warning('PixelSpacing not defined in info structure. Is it a dicom header?\nUsing PixelSpacing = [1 1]');
    pxsp = [1; 1];
end;
slicesp = 0;
if isfield(info, 'SliceLocations') && numel(info.SliceLocations)>1
    slicesp = diff(info.SliceLocations(1:2));
end;
isMOSIAC = isfield(info,'ImageType') && ~isempty(strfind(info.ImageType,'MOSAIC'));
% if (isMOSIAC || slicesp == 0) && isfield(info,'SpacingBetweenSlices')
if isMOSIAC || isfield(info,'SpacingBetweenSlices')    
    slicesp = info.SpacingBetweenSlices;
end
if ndims(info.ImageOrientationPatient)>2
    warning('ImageOrientationPatient not equal for all slices/volumes; using the info of the first slice for all slices/images');
end;
if size(info.ImageOrientationPatient(:,1),1)==9
    orient = reshape(info.ImageOrientationPatient(:,1),3,3)';
    slicedir =orient(3,:);
    orient = orient(1:2,:);
else
    orient = reshape(info.ImageOrientationPatient(:,1),3,2)';
    if isfield(info,'ImagePositionsPatient') && size(info.ImagePositionsPatient,2)>1;
        slicedir = median(diff(info.ImagePositionsPatient,[],2),2)';
        if abs(norm(slicedir)-slicesp)<.05*slicesp || ~isfield(info,'SpacingBetweenSlices')
            slicesp = 1;
        else
            % ImagePositionsPatient does not give slice locations, so better trust SpacingBetweenSlices
            % happens for 'mosaic' images (in DTI/fMRI)
            if isfield(info,'Private_0029_1010') && isfield( info.Private_0029_1010,'SliceNormalVector')
                slicedir = zeros(1,3);
                for k=1:3;slicedir(k) = str2double(info.Private_0029_1010.SliceNormalVector.value{k});end;
            else
                slicedir = -cross(orient(2,:),orient(1,:)); % TODO: validate sign of slicedir.
            end;
        end;
    else
        slicedir = -cross(orient(2,:),orient(1,:));
    end;
end;

transform = zeros(4,3);
if ndims(pxsp)>2
    if max( std(pxsp(:,:,:),[],3)./mean(pxsp(:,:,:),3) )<.001
        pxsp = mean( pxsp(:,:,:) ,3);
    else
        error('pixel spacing of the different slices not equal; cannot combine into a single volume');
    end;
end;
transform([2 1],:) = diag(pxsp) * orient;
transform(3,:) = slicesp * slicedir;
transform(4,:) = info.ImagePositionPatient(:,1)';
transform(:,[1 2]) = transform(:,[2 1]);
transform(4,4) = 1;

% Compensate for dicom having center of first voxel at location [0 0 0] and MATLAB at [1 1 1]
transform(4,:) = [-1 -1 -1 1]*transform;