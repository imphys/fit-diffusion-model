function write_dicom3D(filename, img, info, objT, seriesnr, frameOfReferenceUID)
% write_dicom3D(filename, img, info, Tf, SeriesNumber, frameOfReferenceUID)
% writes a 3D dicom image
% 
% filename : path and base of file name, the slice number and extension
%           '.dcm' are automatically added.
% img      : 3D matrix with the image
% info     : a dicom info structure that is used to provide all dicom
%            information. Typically from the source image (e.g. loaded by my loaddicom.m function)
%            I modify some fields:
%            - A random series number is created.
%            - orientations and slice locations are set from Tf
% Tf       : affine (4 x 4) transformation matrix that specifies the
%            transformation from image index to object space.
%            img(x1,x2,x3) has physical location [x1 x2 x3 1]*Tf
% SeriesNumber : optional: new series number; default: random value between 0 and 1000.
% frameOfReferenceUID: Optional argument that specifies the frame of reference unique identifier.
%            Note that typically this will be the same within a single session, 
%            but might differ between sessions.
%
% Created by Dirk Poot, 1-2-2011
% Erasmus MC

% Check inputs:
if any(size(objT)~=[4 4]) || any(objT(1:end-1,end)~=0) || objT(end,end)~=1
    warning('write_dicom3D:badObjT','Input objT should be a 4 x 4 affine transform matrix with last column [0; 0; 0; 1]');
end;

if nargin>=5 && ~isempty(seriesnr)
    if numel(seriesnr)~=1
        warning('write_dicom3D:badSeriesNumber','New SeriesNumber is not a scalar, this will probably create an error when writing the file; did you confuse it with "frameOfReferenceUID"');
    end;
else
    seriesnr=round(rand(1)*1000); % random series number
end;

% Compensate for MATLAB having first voxel at [1 1 1] and dicom at [0 0 0]
objT(4,:) = [1 1 1 1]*objT;

% remove some fields that might have been added:
rmfields = {'datafiles','SliceLocations','ImagePositionsPatient','Stimulus','Dataset'};
foundfields = isfield(info,rmfields);
infres = rmfield(info, rmfields(foundfields));

if nargin>=6 && ~isempty(frameOfReferenceUID)
    infres.frameOfReferenceUID = frameOfReferenceUID;
end;
infres.SeriesInstanceUID = dicomuid; % create new series instance UID
infres.SeriesNumber=seriesnr;
infres.AcquisitionNumber=seriesnr;
infres.SpacingBetweenSlices = norm(objT(3,:));
infres.SliceThickness = norm(objT(3,:));
infres.PixelSpacing = [norm(objT(2,:)),norm(objT(1,:))];
if ~isfield(infres,'Manufacturer')
    infres.Manufacturer=[];
end;
infres.Manufacturer = ['Matlab Convert, from ' infres.Manufacturer];
infres.ImageOrientationPatient = reshape((diag(1./infres.PixelSpacing)*objT([2 1],[2 1 3]))',6,1);

% undo loaddicom slice & volume expansion:
fieldnm = fieldnames(infres);
selfields  = false(1,numel(fieldnm));
for j2 = 1:numel(fieldnm)
    selfields(j2) = ndims( infres.( fieldnm{j2} ) ) > 2;
end;
selfields = find(selfields);
if ~isempty(selfields)
    inf_orig = infres;
end;
           
            
progressbar('start',size(img,3),'saving dicom','EstTimeLeft','on');
for k=1:size(img,3)
    dcmname = [filename sprintf('%03d',k) '.dcm'];
    for j2 = selfields;
        if iscell(inf_orig.( fieldnm{j2} ))
            infres.( fieldnm{j2} ) = inf_orig.( fieldnm{j2} ){:,:,min(k,end)};
        else
            infres.( fieldnm{j2} ) = inf_orig.( fieldnm{j2} )(:,:,min(k,end));
        end;
    end;
    infres.SliceLocation = k*infres.SpacingBetweenSlices ;
    infres.ImagePositionPatient = [0 0 k-1 1]*objT(:,[2 1 3]);
    infres.InstanceNumber = k;
    % k, infres.SliceLocation, infres.ImagePositionPatient, infres.SpacingBetweenSlices 
    % note: also store private fields. However note that data type is cast to unknown (typically interpreted as uint8 or char).
    % use createmode 'copy' since otherwise rescaleSlope & rescaleIntercept
    % are not stored.
%     dicomwrite( uint16( img(:,:,k) ), dcmname, infres, 'CreateMode','Copy','WritePrivate', true );
    dicomwrite( uint16( img(:,:,k) ), dcmname, infres);
    progressbar(k);
end;
progressbar('ready')
