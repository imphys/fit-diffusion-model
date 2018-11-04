function [transform] = affineTransformFromDicominfo(info)
% [transform] = affineTransformFromDicominfo(info)
% [transform] = affineTransformFromDicominfo(filename)
% Extracts the affine transformation parameters from the dicom info
% structure. 
% INPUTS:
% info  : dicominfo structure of the dicom file from which you want the
%         affine transformation
% filename: string with the name of a dicom file. 
%
% Created by Dirk Poot, Erasmus MC
% 28-5-2010

if ischar(info)
    info = dicominfo(info);
end;

pxsp = info.PixelSpacing;
slicesp = 0;
if isfield(info, 'SliceLocations')
    slicesp = diff(info.SliceLocations(1:2));
end;
if slicesp == 0;
    slicesp = info.SpacingBetweenSlices;
end
orient = reshape(info.ImageOrientationPatient,3,2)';
if isfield(info,'ImagePositionsPatient');
    slicedir = median(diff(info.ImagePositionsPatient,[],2),2)';
    slicesp = 1;
else
    slicedir = -cross(orient(2,:),orient(1,:));
end;

transform = zeros(4,3);

transform([2 1],:) = diag(pxsp) * orient;
transform(3,:) = slicesp * slicedir;
transform(4,:) = info.ImagePositionPatient';
transform(:,[1 2]) = transform(:,[2 1]);
transform(4,4) = 1;