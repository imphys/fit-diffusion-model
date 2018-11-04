function [M, Mx, My, sliceno, masksource] = drawROI( images, copy2slices, destfolder )
% [M, Mx, My, sliceno, masksource] = drawROI( images, copy2slices)
% Asks the user to draw an ROI or load from disk.
% images : a set of images as can be input for fit_MRI
%         so first dimension is different contrasts
%         and ROI drawn on 2nd and 3rd dimension.
%         This ROI is copied to all slices (4th dimension) if copy2slices
%         == true;
%
% Created by Dirk Poot, TUDelft, 27-3-2014
figno = figure; 
descr = 'All images on which the fit will be performed. Please draw the ROI in one slice (in dimension 2 - 3).';
if copy2slices
    descr = [descr ' This ROI will be copied to all other slices.'];
else
    descr = [descr ' The slice on which you drawed the ROI will be the slice on which the fitting is performed.'];
end;
newline = sprintf('\n');
descr  = [descr newline 'Alternatively you can close this figure and a file browser allows you to select a mat file from which the variable mask is loaded as mask.'];
figno = imagebrowse(images,[],'description',descr);
ud = get(figno,'userdata');
axid = ud.subfigs( ( ud.subfigdims(:,1)==2 & ud.subfigdims(:,2)==3 ) );
set(figno,'CurrentAxes',axid,'name','IR - Images');

[Mslice, Mx, My] =roipoly;
if ~ishandle(figno)
    % load from disk:
    if nargin<3
        destfolder = '';
    end
    [FileName,PathName,FilterIndex] = uigetfile([destfolder '*.mat']);
    masksource = [PathName FileName];
    loadmask = load(masksource,'mask','mask_x','mask_y');
    szimg = size(images);
    szmask = size(loadmask.mask);
    if ~isequal( szimg(2:end) , szmask)
        error('The mask loaded from the file should have the same size as the current image');
    end;
    M = loadmask.mask;
    if isfield(loadmask,'mask_x') && isfield(loadmask,'mask_y')
        Mx = loadmask.mask_x;
        My = loadmask.mask_y;
    end;
    szmask(end+1)=1; % make sure a 3rd element exists
    sliceno = find(any(reshape(M,[prod(szmask(1:2)) szmask(3)]),1));
else
    M = repmat( Mslice,[1 1 size(images,4) ]);
    if copy2slices ||size(images,4)==1
        sliceno = 1:size(images,4);
    else
        sliceno = round( get(ud.sliders(4),'value')) ;
        M(:)=false;M(:,:,sliceno)=Mslice;
    end;
    close(figno);
    masksource = 'user drawn';
end;
