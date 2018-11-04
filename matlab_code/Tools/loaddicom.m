function [M_slice, info] = loaddicom( dataset, doRecurseDir, useUint16) 
% [M_slice, info] = loaddata( Path [, searchSubDirectories [, return_uint16]] )  
%
% This function loads (f)MRI/DWI dicom datasets.
%   
% Inputs:
% Path : 
%   Loads the dicom data in the directory given by Path. 
%   Leave empty to interactively locate a directory.
%   IMPORTANT: include a final fileseparator if just a directory name.
%   When only 1 series is in the directory this series is loaded, otherwize
%   you get a question which series to load. This question can be avoided 
%   by specifiing the number or description of the series as
%   extension. Example: 
%       Data = loaddata('c:\temp\.3'); 
%       Loads dicom series number 3 from directory 'c:\temp' into the variable Data.
%
%   With the extension '.series' a structure with the available series is returned (instead of an image).
%   Note that the available files and some of the the dicom header information 
%   is cached into 'DicomInfoCache*.mat'. If files are added/removed from the 
%   directory or this header information changes, removing the cache file forces 
%   re-scanning of all dicom files in the directory.
% 
% searchSubDirectories :  default = false
%   If true all subdirectories are recursively searched for dicom files. 
% return_uint16 :  default = false
%   Boolean to specify whether the data should be loaded into an uint16
%   variable to save memory. By default a double precision array is retured.
%
% Outputs: 
%  M_slice :
%   4-D Data matrix, dimensions: [x y z time]
%   For multi slice images, the z-direction might be the 4th dimension.
%  info :
%   dicominfo structure of the first file in the series.
%   The fields I do add:
%   'Dataset'  : Call loaddata with this parameter to get the same dataset
%                again, also correct for interactively selected datasets.
%   For images that might be DTI images:
%   'dti_bMatrix', 'dti_bValue', 'dti_grad_dir', 
%   'dti_bValue_alt', 'dti_grad_dir_alt'
%                Specify the diffusion weighting applied. The 'alt*' variables
%                are alternative interpretations/ read from different dicom fields.
%                They are not always equal to the probably slightly more 
%                reliable non 'alt*' versions.
%
% Copyright D.H.J. Poot, University of Antwerp & Erasmus MC

% last edit: 
%    10-05-2006  : Added Dataset & Stimulus to info structure.
%     9-05-2006  : Added support for loading dicom directories.
%                  Added support for returning uint16 data.
%    18-06-2008  : Added support for loading (not previously known) files
%                  with manual method specification 
%    26-05-2011  : Remove old features ('known datasets') and update help text.
%    27-11-2012  : Use ImagePositionPatient for main ordering & detection
%                  of 3D volumes. Also now the non stationary dicom fields are
%                  stored in array's in the 3rd or 4th dimension.
 
persistent cache_infofilename cache_files; % preserve last dicominfo, to avoid excessive loading and excessive 'dicominfo' calls when the cache cannot be saved.

tryToSave = 1;      % try to save dicom info from files. 
displayNotDicomFiles = true;

if nargin<1 || isempty(dataset)
    % If called without argument, interactively select dataset directory:
    dataset = uigetdir;
    if isequal(dataset, 0)
        disp('Dataset loading cancelled per user request');
        return;
    end;
    dataset= [dataset filesep];
end;
% Set default option for recurseDir when not provided.
if nargin<2 || isempty(doRecurseDir)
    doRecurseDir = false;
end;
% set default to not use uint16
if nargin<3 || isempty(useUint16);
    useUint16 = false;
end;

% Split pathname:
[dsPath, dsName, dsExt]= fileparts( dataset ); 

%  Load all dicom files in current directory and get info, get different
%  series
if doRecurseDir
    dicomCacheName = [dsPath filesep 'DicomInfoCache_r.mat'];
else
    dicomCacheName = [dsPath filesep 'DicomInfoCache.mat'];
end;
% Get dicominfo from internal cache, cache file, or re-create:
if isequal( cache_infofilename, dicomCacheName)
    files = cache_files;
elseif exist( dicomCacheName,'file');
    load(dicomCacheName,'files')
    cache_infofilename = dicomCacheName;
    cache_files = files;
else
    % Read all filenames in the specified directory: 
    progressbar('start',[],'Loading dicom datset.','EstTimeLeft','on');
    
    files = dir([dsPath filesep '*']);
    if doRecurseDir
        k=1;
        while k<=length(files)
            % recurse over all subdirectories
            if files(k).isdir && ~(isequal(files(k).name,'.') || isequal(files(k).name,'..')) 
                filesadd = dir([dsPath filesep files(k).name filesep '*']);
                l=1;
                % add current path to all new filenames:
                while l<=numel(filesadd)
                    if isequal(filesadd(l).name,'.') || isequal(filesadd(l).name,'..')
                        filesadd(l) = [];
                    else
                        filesadd(l).name = [files(k).name filesep filesadd(l).name];
                        l=l+1;
                    end;
                end;
                files = [files;filesadd]; %#ok<AGROW> : we don't know the final number of files.
            end;
            k=k+1;
        end;
    end;
    files([files.isdir])=[]; % remove directories from file list.
    if numel(files)>1
        % initialize fields.
        files(1).SeriesNumber = [];
        files(1).InstanceNumber = [];
        files(1).SeriesDescription = [];
    else
        progressbar('ready');
        error('No files found in directory : %s\nSet searchSubDirectories to also search subdirectories.', dsPath );
    end;
    progressbar('clearTimeEstimation');
    
    % For each file check if it is a valid dicom file:
    for k=1:length(files)
%                 waitbar(k/length(files), WB, 'Reading fileinfo.');
        progressbar(k/length(files),'Reading fileinfo.');
        if files(k).bytes>128 %& length(regexp(files(k).name,'\d{4,}','once'))>0% only use nonempty files with at least 4 adjacent digits in the name
           try
               info = dicominfo([dsPath filesep files(k).name]);
               if ~ischar(info.SeriesNumber)
                   files(k).SeriesNumber = info.SeriesNumber;
               else
                   files(k).SeriesNumber = -2;
               end;
               files(k).InstanceNumber = info.InstanceNumber;
               if isfield(info,'SeriesDescription')
                   files(k).SeriesDescription = info.SeriesDescription;
               else
                   files(k).SeriesDescription = 'empty';
               end;
               if isfield(info,'SeriesInstanceUID');
                   files(k).SeriesInstanceUID = info.SeriesInstanceUID;
               end;
               if isfield(info, 'ImagePositionPatient')
                   files(k).ImagePositionPatient = info.ImagePositionPatient;
                   files(k).ImageOrientationPatient = info.ImageOrientationPatient;
               end;
%                        files(k).AcquisitionNumber = info.AcquisitionNumber;
           catch
               files(k).SeriesNumber = -1; % no series.
           end;
        else
            files(k).SeriesNumber = -1; % no series.
        end;
    end;
    if isfield(files,'SeriesInstanceUID')
        for k=1:numel(files)
            if isempty(files(k).SeriesInstanceUID)
                % make sure that SeriesInstanceUID is filled for all files, when SeriesInstanceUID exists.
                files(k).SeriesInstanceUID = num2str( files(k).SeriesNumber );
            end;
        end;
    end;
    % fill memory cache :
    cache_infofilename = dicomCacheName;
    cache_files = files;
    if tryToSave
        %save cache file:
        dicomCacheVersion = '2';
        try 
            save( dicomCacheName, 'files','dicomCacheVersion');
        catch
            disp('Cannot write dicominfo cache file.')
        end;
    end;
    progressbar('ready');
end;

% Oleh:
% files = files(cellfun(@(x)isempty(strfind(x,'modvars')),{files.name},'UniformOutput',true));


% --  PROCEED WITH IDENIFYING, SORTING, AND LOADING THE REQUESTED DATASERIES ----


% processing fileinfo:
seriesNrs = [files.SeriesNumber];
if isfield(files,'SeriesInstanceUID')
    seriesIds = {files.SeriesInstanceUID};
else
    seriesIds = seriesNrs;
end;
    
% q=([files.InstanceNumber])-32;seriesNrs = mod(q,66); % hack for Philips DTI data of Bianca.

if any(seriesNrs==-2)
    % for valid dicom files which have no SeriesNumber specified, set to value higher than highest found:
    % different number for each series.
    seriesNrs(seriesNrs==-2) = max(0,max(seriesNrs))+(1:nnz(seriesNrs==-2));
end;

% find the unique series Id's, find the corresponding series numbers and
% sort them:
[uniqueSeriesIds , idIdxToImg, imgToIdIdx] = unique(seriesIds);
% uniqueSeriesNrs = sort(seriesNrs);
% uniqueSeriesNrs = uniqueSeriesNrs([1 diff(uniqueSeriesNrs)] & uniqueSeriesNrs~=-1);
uniqueSeriesNrs = seriesNrs(idIdxToImg);
[uniqueSeriesNrs , remapm] = sort(uniqueSeriesNrs);
uniqueSeriesIds = uniqueSeriesIds(remapm);
idIdxToImg = idIdxToImg(remapm);
iremapm = remapm;
iremapm(remapm) = 1:numel(remapm);
imgToIdIdx = iremapm(imgToIdIdx);
tmp = (uniqueSeriesNrs ==-1);
if any(tmp)
    % remove files with id ==-1, which is created when a file is not a
    % (valid) dicom file.
    uniqueSeriesNrs = uniqueSeriesNrs(~tmp);
    uniqueSeriesIds = uniqueSeriesIds(~tmp);
    idIdxToImg = idIdxToImg(~tmp);
    newidx = cumsum(~tmp);newidx(tmp)=0;
    imgToIdIdx = newidx(imgToIdIdx);
end;
if displayNotDicomFiles && any(seriesNrs==-1)
    disp('files found that appear not to be dicom files:');
    for k=find(seriesNrs==-1)
        tst = [dsPath filesep files(k).name];
        [dum1,dum2,tstext] = fileparts(tst);
        if isequal(tstext,'.dcm')
            [dum4, dum5] = fileparts(dsPath);
            disp([dsPath filesep files(k).name]);
%                     disp(['/mnt/edra_ergo/Rotterdam_Study/' dum5 '/' strrep(files(k).name,'\','/')]);
        end;
    end;
end;

LoadSeries = -2; % Initial value, should never be used as series index. 
returnAvailableSeries = false;
if length(dsExt)>1
    if strcmpi(dsExt,'.series')
        returnAvailableSeries = true;
    else
        LoadSeries = str2num(dsExt(2:end));
    end
end;
if any( uniqueSeriesNrs==LoadSeries(1) )
    % Select serie that is requested (when it exists):
    tmp = find( uniqueSeriesNrs==LoadSeries(1) );
    if numel(LoadSeries)>1
        tmp = tmp( LoadSeries(2) ); % select from duplicate series numbers. 
        % these have been identified by their (different) seriesUID
    end;
    uniqueSeriesIds = uniqueSeriesIds( tmp );
    uniqueSeriesNrs = uniqueSeriesNrs( tmp );
    idIdxToImg = idIdxToImg( tmp );
    [dummy, imgToIdIdx] = ismember(imgToIdIdx,tmp);
    disp('Loading series:');
else
    if ~returnAvailableSeries
        disp('Available series:');
    end;
end;

% Create structure containing info about each series:
DicomSeries = struct('SeriesNumber',cell(numel(uniqueSeriesNrs),1),'SeriesDescription',[],'files',[],'InstanceNumbers',[],'instanceStep',[]);

numduplicates = 0;
for k = 1:numel(uniqueSeriesIds)
    seriesNr = uniqueSeriesNrs(k);
    if nnz(uniqueSeriesNrs==seriesNr)>1
        seriesNrStr = [num2str(seriesNr) ',' num2str(find(find(uniqueSeriesNrs==seriesNr)==k))];
    else
        seriesNrStr = num2str(seriesNr);
    end;
    DicomSeries(k).SeriesNumber = seriesNr;
    % Find all files belonging to this series:
    seriesSel = imgToIdIdx==k;
    instanceNrs = sort([files(seriesSel).InstanceNumber]);
    if length(instanceNrs)==1
        instanceStep = 1;
    else
        % Fill in instanceNrs when not provided in the files:
        if ~any(diff(instanceNrs)) && isfield(files, 'ImagePositionPatient')
            % some datasets don't have different instanceNrs to
            % differentiate slices. Try to guess from imageposition.
            slicepos = [files(seriesSel).ImagePositionPatient];
            orient = [files(seriesSel).ImageOrientationPatient];
            if any(any(diff(orient,[],2)))
                error( 'image orientation not constant in 3D image');
            end;
            orient = reshape(orient(:,1),3,2);
            slicedir = -cross(orient(:,2),orient(:,1));
            [sliceip, instanceNrs] = sort( slicedir' * slicepos );
            if any(abs( diff(sliceip)/median(diff(sliceip)) -1 ) > .1 )
                warning('The image slices ordered by slice position do not seem to be regular spaced.' );
            end;
            fsersel = find(seriesSel);
            for i_nr = 1 : numel(instanceNrs)
                files( fsersel( instanceNrs( i_nr ) ) ).InstanceNumber =  i_nr ;
            end;
            instanceNrs = 1:numel(instanceNrs);
        end;
        instanceStep = min(diff(instanceNrs));
        if instanceStep==0
            numduplicates = numduplicates + sum(diff(instanceNrs(:))==0);
            instanceNrs = instanceNrs([true ; diff(instanceNrs(:))>0]);
            instanceStep = min(diff(instanceNrs));
            if isempty(instanceStep)
                instanceStep = 0;
            end;
        end;
    end;
    DicomSeries(k).instanceStep = instanceStep;
    fsersel = find(seriesSel);
    DicomSeries(k).filesidx = fsersel;
    for l_ind= 1:numel(fsersel);
        l = fsersel(l_ind);
        if isempty(files(l).InstanceNumber)
            fileNr = l_ind; % Instance number should be specified.
        else
            fileNr = round((files(l).InstanceNumber - instanceNrs(1))/max(1,instanceStep))+1;
        end;
        DicomSeries(k).files{fileNr} = [dsPath filesep files(l).name];
        DicomSeries(k).SeriesDescription = files(l).SeriesDescription;
    end;

    if length(DicomSeries(k).files)>0 && ~returnAvailableSeries
        % Print info for interactive selection of series
        fprintf('Series nr. %s; %3d images; Description: %s',seriesNrStr,length(DicomSeries(k).files),DicomSeries(k).SeriesDescription);
        missing = [];
        for l=1:length(DicomSeries(k).files)
            if isempty(DicomSeries(k).files{l})
                missing = [missing l];
            end;
        end;
        if length(missing)==0
            fprintf('\n');
        else
            missing = [missing missing(end)+10];
            m = 2;st = 1;s = '';
            while m<=length(missing)
                if missing(m) - missing(st) > m-st
                    if m==st+1
                        s = [s sprintf('  %d',missing(st))];
                    else
                        s = [s sprintf('  %d-%d',missing([st m-1]))];
                    end;
                    st=m;
                end;
                m=m+1;
            end;
            fprintf(', there are %d files missing in this series. (numbers:%s)\n',length(missing)-1,s);%sprintf(' %d',missing));
        end;
    end;
end;
% clear loop variables so we dont accidently use them afterwards:
clear seriesNr seriesNrStr seriesSel fsersel instanceNrs instanceStep fileNr
if numduplicates >0
    warning('Loaddata:Duplicates',['There appear to be ' num2str( numduplicates ) ' duplicate images, using only the first instance found.']);
end;
if returnAvailableSeries
    M_slice = DicomSeries;
    return;
end;

%         if isempty(LoadSeries) | round(LoadSeries)~=LoadSeries | LoadSeries<=0 | LoadSeries>length(DicomSeries) 
% If multiple series are still present, ask user to select series that
% he/she wants to load:
if length(uniqueSeriesNrs)==0
    error(['No dicom files found in directory :' dsPath ]);
    return;
elseif length(uniqueSeriesNrs)>1
    LoadSeries = input('Which series do you want to load (Nr)? : ','s');
    LoadSeries = str2num(LoadSeries);
    if isempty(LoadSeries) || numel(LoadSeries)>2
        error('invalid series selected');
    end;
    dataset = [dsPath filesep '.' num2str(LoadSeries(1))];
    tmp = find([DicomSeries.SeriesNumber]==LoadSeries(1));
    if numel(LoadSeries)==2
        dataset = [dataset ',' num2str(LoadSeries(2))];
        tmp = tmp(LoadSeries(2));
    end;
    DicomSeries = DicomSeries(tmp);
else
    LoadSeries = uniqueSeriesNrs;
end;

LoadSeries = 1;%find([DicomSeries.SeriesNumber]==LoadSeries); % map seriesNumber to index in DicomSeries.
if numel(DicomSeries)~=1
    error('Requested series not found.');
end;
disp(['Loading series ' num2str(DicomSeries.SeriesNumber)]);

numfiles = length(DicomSeries.files);

% get info from first file of the series (most info will be constant; not fully checked).
info = dicominfo( DicomSeries.files{1} ); 
info_sl1 = info;
info_fields = fieldnames(info);
iterate_first_over_z_dim = true;
didwarn_improperlocations = false;
% Extract/Compute size of return image 
if isfield(info,'ImageType') && ~isempty(strfind(info.ImageType,'MOSAIC'))
    % data is probably in mosaic formL
    info = parseSiemensCSAheader(info);
    zSize = str2double( info.Private_0029_1010.NumberOfImagesInMosaic.value{1} );
    nTiles= ceil( sqrt( zSize )); % n_row_blocks == n_col_blocks;
    
    xSize = double(info.Rows)/ nTiles;
    ySize = double(info.Columns)/ nTiles;

    % Correct ImagePositionPatient to first voxel of first slice
    % (from: http://nipy.sourceforge.net/nibabel/dicom/dicom_mosaic.html#dicom-mosaic)
    % Since apparently it is stored incorrectly.
    F = reshape(info.ImageOrientationPatient,[3 2]);
    slicedir = zeros(3,1);
    for k=1:3;slicedir(k) = str2double( info.Private_0029_1010.SliceNormalVector.value{k} );end;
    RS = [F(:,[2 1]) slicedir] * diag([info.PixelSpacing; info.SpacingBetweenSlices]);
    info.ImagePositionPatient = info.ImagePositionPatient + RS(:,1:2)*[(double(info.Rows) - xSize)/2;(double(info.Columns) - ySize)/2];
    nb_vol = numfiles;
else
    zSize = 1;
    nTiles = 1;
    ySize = info.Rows;
    xSize = info.Columns;
    nb_vol = numfiles;
    if isfield(files, 'ImagePositionPatient')
        slicepos = [files(DicomSeries.filesidx).ImagePositionPatient];
        [uniqueslicepos, slicepos_m, slicepos_n] = unique(slicepos','rows');
        expectednumslices = size( uniqueslicepos ,1);
        if mod( numfiles, expectednumslices)==0 
            zSize = expectednumslices;
            nb_vol = numfiles/expectednumslices;
            if zSize>1 && diff( slicepos_n(1:2) )==0
                iterate_first_over_z_dim = false;
            end;
        else
            warning('LOADDICOM:IMPROPER3D','The number of images in this series divided by the number of slice positions is not integer. Thus I cannot form a proper 3D image and will return each slice as different ''volume''.');
            didwarn_improperlocations = true;
        end;
    elseif isfield(info,'NumberOfFrames')
        if mod(numfiles, info.NumberOfFrames)==0
            zSize = info.NumberOfFrames;
            nb_vol = numfiles/zSize;
        end;
    elseif isfield(info,'ImagesInAcquisition');
        if mod(numfiles, info.ImagesInAcquisition)==0
            zSize = info.ImagesInAcquisition;
            nb_vol = numfiles/zSize;
        end;
    end;
end;
fprintf('DICOM image parameters: y size = %d, x size = %d, z-size <= %d, number of tiles = %d x %d\n',ySize, xSize,zSize,nTiles,nTiles);

ySize  = double(ySize);
xSize  = double(xSize);

zSize = double(zSize);
nTiles = double(nTiles);
info.datafiles = cell(numfiles,1);

% Should check if   Private_0019_10xx_Creator == 'SIEMENS MR HEADER'
if isfield(info,'Manufacturer') && ~isempty(info.Manufacturer)
vendor = info.Manufacturer;
if isequal(vendor(1:2),'GE')
    hasgradient = ((numel(info.SoftwareVersion)==2) && all(info.SoftwareVersion>='11') || numel(info.SoftwareVersion)~=2 ) &&  isfield(info,'Private_0019_10e0') && numel(info.Private_0019_10e0)==1 && info.Private_0019_10e0>0;
elseif isfield(info,'Private_2001_10xx_Creator') && strcmpi(info.Private_2001_10xx_Creator(1:min(7,end)),'Philips')
    hasgradient = isfield(info,'Private_2001_1004');
    vendor = 'Philips';
else
    hasgradient = isfield(info,'Private_0019_100c');
end;
else
    hasgradient =false;
end;
if isfield(info,'SliceLocation')
    info.SliceLocations = nan(1, zSize);
end;
if isfield(info,'ImagePositionPatient')
    info.ImagePositionsPatient = nan( 3, zSize);
end;
if hasgradient
    fprintf(['reading ' vendor ' diffusion gradients, for DWI''s']);
    info.dti_grad_dir = nan(nb_vol,3);
    info.dti_grad_dir_alt = nan(nb_vol,3);
    info.dti_grad_dir_alt2 = nan(nb_vol,3);
    info.dti_bValue = nan(nb_vol,1);
    info.dti_bValue_alt = nan(nb_vol,1);
    info.dti_bMatrix = nan(nb_vol,6);
    dti_bValue_alt = nan;
end;    

% be memory efficient, use uint16 (later changed to double when needed)
M_slice = repmat( uint16(0), double([ySize xSize zSize nb_vol])); 
progressbar('start',numfiles,'Reading images','EstTimeLeft','on');
for file_idx = 1:numfiles
    str_read = DicomSeries.files{file_idx};
    progressbar(file_idx, sprintf('Busy with file %d',file_idx));
    if ~exist(str_read,'file')
        disp(['Slice/volume ' num2str(file_idx) ' not present in this series, replacing by zeros in output.']);
        continue;
    end;
    M_vol_temp = dicomread(str_read);
    di = dicominfo(str_read);
    di_fields = fieldnames(di);

    %dicom_slice    
    vol_idx = file_idx;
    slice_idx = 1;
    if (nTiles==1) %&& ndims(M_vol_temp)>2
        if size(M_vol_temp,3)==1 && zSize>1
            if iterate_first_over_z_dim
                M_slice(:, :, file_idx) = M_vol_temp(:,:);
                vol_idx = ceil( file_idx / zSize );
                slice_idx = mod( file_idx-1 , zSize )+1;
            else
                vol_idx = mod( file_idx-1, nb_vol )+1;
                slice_idx = ceil( file_idx / nb_vol );
                M_slice(:, :, slice_idx + zSize * (vol_idx-1) ) = M_vol_temp(:,:);
            end;
        else
            M_slice(:, :, :, file_idx) = M_vol_temp(:,:,:);
        end;
    else
        % MOSAIC
        for co_sl = 1:zSize
            row_start=floor((co_sl-1)/nTiles)*ySize+1;
            row_end=floor((co_sl-1)/nTiles)*ySize+ySize;
            col_start=mod(co_sl-1,nTiles)*xSize+1;
            col_end=mod(co_sl-1,nTiles)*xSize+xSize;
            M_slice_temp=M_vol_temp(row_start:row_end,col_start:col_end);
            M_slice(:, :, co_sl, file_idx) = reshape(M_slice_temp,ySize,xSize,1);      
        end;
    end;        
    
    [field_is]= ismember(info_fields, di_fields );
    if ~all( field_is ) || (numel(di_fields)~=numel(info_fields)) % ~isequal(di_fields,info_fields)
        warning(['File ' num2str(file_idx) ' does not have '  num2str(nnz(~field_is)) ' fields that are present in the first image and has ' num2str(numel(di_fields)-nnz(field_is)) ' that are not present in the first scan. New field not added and not present fields keep value of first scan.']);
    end;
    for fieldidx  = find( field_is )' %1: numel(info_fields)
        if ~isequal( info_sl1.(info_fields{fieldidx}), di.(info_fields{fieldidx}) );
            % value in new scan is not equal to first scan, or there have
            % been differences in previous files so this field should be set for each image.
            
            fn = info_fields{fieldidx} ; % get field name in convienient variable name.
            if ischar( info.(fn) )  
                % convert character arrays to scalar cell array:
                info.(fn) = {info.( fn )};
            end;
            sl_p = min( size( info.(fn) ,3) , slice_idx ); % get slice position
            if vol_idx>1 && size( info.(fn) ,4)==1 && ~isequal( info.(fn)(:,:, sl_p,1), di.(fn) ) 
                % Expand field in volumes dimension:
                info.( fn ) = cat(4, info.(fn), repmat( info.( fn )(:,:,1,1) , [1 1 size( info.(fn) ,3) nb_vol-1] ) );
            end;
            v_p = min( size( info.(fn) ,4) , vol_idx);
            if slice_idx>1 && size(info.(fn),3)==1 && ~isequal( info.(fn)(:,:,1,v_p), di.(fn) )
                % Expand in slice (/z) dimension:
                info.(fn) = cat(3, info.(fn), repmat( info.(fn)(:,:,1,1) , [1 1 zSize-1 size( info.(fn) ,4)] ) );
                sl_p = slice_idx ;
            end;
            if ischar( di.(fn) ) || iscell(info.(fn))
                info.(fn)(:,:, sl_p, v_p ) = {di.(fn)};
            else
                if size(info.(fn),1)~=size(di.(fn),1) || size(info.(fn),2)~=size(di.(fn),2) 
                    warning(['Field "' fn '" has inconsistent size, returning as cell array']);
                    info.(fn) = mat2cell(info.(fn), size(info.(fn),1), size(info.(fn),2), ones(1,size(info.(fn),3)), ones(1,size(info.(fn),4)) );
                    info.(fn)(:,:, sl_p, v_p ) = {di.(fn)};
                else
                    info.(fn)(:,:, sl_p, v_p ) = di.(fn);
                end;
            end;
        end;
    end;

    info.datafiles{file_idx} = str_read;
    
    if isfield(di,'SliceLocation');
        if ~isequal( info.SliceLocations(:,slice_idx), di.SliceLocation)
            if isnan( info.SliceLocations(:,slice_idx) )
                info.SliceLocations(:, slice_idx) = di.SliceLocation;
            else
                if ~didwarn_improperlocations
                    warning('LOADDICOM:ImproperSliceLocations','Slice locations in each volume should be equal');
                    didwarn_improperlocations = true;
                end;
            end;
        end;
    end;
    if isfield(di,'ImagePositionPatient');
        if ~isequal( info.ImagePositionsPatient(:,slice_idx) , di.ImagePositionPatient)
            if all(isnan( info.ImagePositionsPatient(:,slice_idx) ) )
                info.ImagePositionsPatient(:,slice_idx )= di.ImagePositionPatient;
            else
                if ~didwarn_improperlocations
                    warning('LOADDICOM:NonConstantVolumeLocations','ImagePositionPatient should be equal in each volume');
                    didwarn_improperlocations = true;
                end;
            end;
        end;
    end;
    
    if hasgradient
        % TODO: check if rotation of field of view is treated consistently in gradient direction
        if isequal(vendor(1:2),'GE')
            % GE:
            dti_bValue = di.Private_0043_1039(1);
            dti_grad_dir = [di.Private_0019_10bb di.Private_0019_10bc di.Private_0019_10bd];
        elseif isequal(vendor(1:2),'Ph') 
            % Philips:
            dti_bValue = di.Private_2001_1003;
            if strcmp( di.Private_2001_1004 ,'O')
                dti_grad_dir = [di.Private_2005_10b0 di.Private_2005_10b1 di.Private_2005_10b2];
            else
                % no diffusion direction (b-value ==0, or other problem?)
                dti_grad_dir=0;
                if info.dti_bValue~=0
                    warning('no diffusion direction specified, but b-value non zero');
                end;
            end;
        else
            % Siemens
            dti_bValue_alt = di.Private_0019_100c;
            if isfield(di, 'Private_0019_1027')
                info.dti_bMatrix(vol_idx,:) = di.Private_0019_1027;
                dti_bValue = norm(sqrt(info.dti_bMatrix(vol_idx,[1 4 6])))^2;
                dti_grad_dir = sign(sign(info.dti_bMatrix(vol_idx,(1:3))) + 0.5).*sqrt(info.dti_bMatrix(vol_idx,[1 4 6])/dti_bValue);
                dti_grad_dir = dti_grad_dir([2 1 3]);
            else
                info.dti_bMatrix(vol_idx,:) = nan;
                dti_grad_dir = [1 0 0];
                dti_bValue = 0;
            end;
            if (di.Private_0019_100c ~= 0) && isfield(di,'Private_0019_100e');
                dti_grad_dir_alt =  di.Private_0019_100e([2 1 3])';
            else
                dti_grad_dir_alt = [1 0 0];
            end;
            if isfield(di,'Private_0029_1010')
                info.dti_grad_dir_alt2(vol_idx,[2 1 3]) = GetDiffusionGradientDirection(char(di.Private_0029_1010)');
            end;
            if ~isequal(info.dti_bValue_alt(vol_idx), dti_bValue_alt)
                if ~isnan(info.dti_bValue_alt(vol_idx))
                    error('Each slice in a volume should have the same b-value (alternative)')
                else
                    info.dti_bValue_alt(vol_idx) = dti_bValue_alt;
                end;
            end;            
            if ~isequal(info.dti_grad_dir_alt(vol_idx,:), dti_grad_dir_alt)
                if ~isnan(info.dti_grad_dir_alt(vol_idx,:))
                    error('Each slice in a volume should have the same gradient direction (alternative)')
                else
                    info.dti_grad_dir_alt(vol_idx,:) = dti_grad_dir_alt;
                end;
            end;  
        end;
        if ~isequal( info.dti_bValue(vol_idx) , dti_bValue)
            if ~isnan( info.dti_bValue(vol_idx) )
                error('Each slice in a volume should have the same b-value')
            else
                info.dti_bValue(vol_idx) = dti_bValue;
            end;
        end;
        if ~isequal( info.dti_grad_dir(vol_idx,:) , dti_grad_dir)
            if any(~isnan( info.dti_grad_dir(vol_idx,:) ))
                error('Each slice in a volume should have the same gradient direction')
            else
                info.dti_grad_dir(vol_idx,:) = dti_grad_dir;
            end;
        end;        
    end;
end
if hasgradient
    noteq = find(max(abs(info.dti_grad_dir-info.dti_grad_dir_alt),[],2)>.02);
    if ~isempty(noteq)
        warning('Loaddata:gradientsNotEq',['The directions of ' num2str(numel(noteq)) ' gradients wrong in info.dti_grad_dir or info.dti_grad_dir_alt. (' sprintf(' %d,',noteq) ')']);
    end;
end;
if nTiles^2==zSize
    for k=zSize:-1:1
%             waitbar(1-k/zSize, WB, 'removing empty tiles.');
        progressbar(1-k/zSize, 'removing empty tiles.');
        if ~any(any(any(M_slice(:,:,k,:)))) % remove empty tiles.
            M_slice(:,:,k,:) = [];
        end;
    end;
end;
%     close(WB);
progressbar('ready');


if useUint16 && ~isa(M_slice,'uint16')
    M_slice = uint16(M_slice); % compact.
end;
if ~useUint16 && isa(M_slice,'uint16')
    M_slice = double(M_slice);
end;

info.Dataset = dataset;



function [ifo] = parseSiemensCSAheader(ifo)
%%
if ~isfield(ifo,'Private_0029_1008') || ~isequal(ifo.Private_0029_1008,'IMAGE NUM 4')
    return;
end;
for stridx = 1:2
    if stridx==1
        procfield = 'Private_0029_1010';
    else
        procfield = 'Private_0029_1020';
    end;
    str = ifo.(procfield)';
    reslt = struct;
isCSA2type = isequal(str(1:4),'SV10');

m = [1;2^8;2^16;2^24];    
csa_position = 8;
csa_max_pos = numel(str);

ntags = double(str(csa_position+(1:4)))*m; csa_position= csa_position+4;
if ntags<1 || ntags>128
    continue;
end;
unused = double(str(csa_position+(1:4)))*m; csa_position= csa_position+4;

for k=1:ntags
    name = char(str(csa_position+(1:64))); csa_position= csa_position+64;
    firstzero = find(name==0,1,'first');
    if ~isempty(firstzero)
        name = name(1:firstzero-1);
    end;
    value = struct;
    value.vm = double(str(csa_position+(1:4)))*m; csa_position= csa_position+4;
%     vr = double(str(csa_position+(1:3)))*m(1:3); csa_position= csa_position+4;
    value.vr = char(str(csa_position+(1:4))); csa_position= csa_position+4;
    value.syngodt = double(str(csa_position+(1:4)))*m; csa_position= csa_position+4;
    nitems = double(str(csa_position+(1:4)))*m; csa_position= csa_position+4;
    value.nitems = nitems;
    value.xx = double(str(csa_position+(1:4)))*m; csa_position= csa_position+4;
    value.itemlen = zeros(1,nitems);
    value.itemxx2 = zeros(1,nitems);
    value.itemxx3 = zeros(1,nitems);
    value.itemxx4 = zeros(1,nitems);
    value.value = cell(1,nitems);

    for itidx = 1:nitems
        itemlen = double(str(csa_position+(1:4)))*m; csa_position= csa_position+4;
        value.itemlen(itidx) = itemlen;
        value.itemxx2(itidx) = double(str(csa_position+(1:4)))*m; csa_position= csa_position+4;
        value.itemxx3(itidx) = double(str(csa_position+(1:4)))*m; csa_position= csa_position+4;
        value.itemxx4(itidx) = double(str(csa_position+(1:4)))*m; csa_position= csa_position+4;
        value.value{itidx} = char(str(csa_position+(1:itemlen))); csa_position= csa_position+ceil(itemlen/4)*4;
    end;
    if isfield(reslt, name)
        disp(['Overwriting old value of field "' name '". Old value:']);
        disp(reslt.(name));
    end;
    reslt.(name) = value;
end;
ifo.(procfield) = reslt;
end;