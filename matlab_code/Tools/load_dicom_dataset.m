function [MRimgs, info, T , objT, inputs , relativeSliceWidth, szobj] = load_dicom_dataset( sourcefilesdir , seriesnrs , varargin)
% [MRimgs, info, T , objT, inputs , relativeSliceWidth] = load_dicom_dataset( sourcefilesdir , seriesnumbers/seriesnames , [option1, value1 [, ...]]  )
% OLD:
% function [MRimgs, info, T , objT, inputs , relativeSliceWidth] = load_dicom_dataset( sourcefilesdir , seriesnrs , resample_to_acquisitiongrid , szobj, objvoxel_size
% 
% Loads a number of dicom series, reads the transformations and resamples
% the images to the acquisition grid to remove the zero filling in k-space that often is present.
%
% INPUTS:
%  sourcefilesdir  : Base directory of the dicoms of the scans. The actual dicom files might be inside subdirectories.
%                    if empty, a dialog box opens that requests the user to
%                    specify the folder.
%  seriesnumbers   : vector with integers that specifies which series that should be loaded.
%  seriesnames     : cell array with the 'SeriesDescription' of the series that you want to load. 
%                  : if seriesnumbers is empty: all series are loaded.
% Options:
%  regexpMatchSeries : default:false; if true: use all  seriesnames  as
%                      regular expressions and load all series for which
%                      (at least) 1 matches. E.g: seriesnumbers = {'IR'} to
%                      match to all seriesDescriptions that contain 'IR'.
%                      If seriesnames{i} is 'Ask user' : show a box in
%                      which the user can adjust the current selection (up
%                      to all descriptions i-1).
%  resample_to_acquisitiongrid : default = false; Boolean that specifies whether the images should be resampled 
%                                back to the acquisition grid.
%                                If true, T and inputs is adjusted to make the outputs consistent. (info is not)
%  displayFFTbeforeDownsampling : default=false; when resampling, display
%                                   FFT of image before downsampling to
%                                   inspect if it's done correctly.
%  objectSize      : Size of the object space in voxels. The objectare centered in this 
%  objectVoxelSize : Voxel size of each object voxel, default: minimum voxel spacing in source images (after resampling)
%  objectTransform : The specification of the object space in world coordinates.
%  imagePermute    : default: empty = no permutation. 
%                    Cell array with a permutation vector per MR image (seriesnrs).
%  AcquisitionMatrix: default: empty = no override.
%                    Cell array with an override for each of the acquisition Matrices of the MR images.
%                    Mainly usefull when AcquisitionMatrix is not specified in the dicom header.
%
% OUTPUTS:
% MRimgs	: cell array with the MR images.
% info      : cell array with dicom info of each of the MR images.
% T         : cell array with transformation matrices of each of the images. 
%             If the images are resampled, T is adjusted to represent the resampled images.
% objT      : suggested transform for the object
% inputs    : inputs for 'TomographicMRI_Transf', based on objT
% relativeSliceWidth : relative slice width of each of the MR images.
%                      positive for 2D acquisitions (and when type of acquisition cannot be determined)
%                      negative for 3D acquisitions.
%
% Created by Dirk Poot, Erasmus MC
% 23-2-2011


% % objT            : use direction cosinuses and voxel spacing from this affine transformation matrix
% %                   as objT with which 'inputs' is created.


opts.resample_to_acquisitiongrid = false;
opts.objectSize = [];
opts.objectVoxelSize =[];
opts.objectTransform = [];
opts.imagePermute = [];
opts.regexpMatchSeries = false;
opts.AcquisitionMatrix = [];
opts.displayFFTbeforeDownsampling = false;
opts = parse_defaults_optionvaluepairs(opts, varargin{:});

if isempty(sourcefilesdir)
    sourcefilesdir = uigetdir;
    sourcefilesdir = [sourcefilesdir filesep];
    disp(['loading from : "'  sourcefilesdir '"'])
end;
ndim = 3; % number of dimensions, can be increased by first image.
ndatadim = 3;
progressbar('start',[-1 numel(seriesnrs)]);
if isempty(seriesnrs)
    disp('No series numbers provided, loading all dicom series in directory');
    [seriesinfo] = loaddicom([sourcefilesdir '.series'],1,1); 
    seriesnrs =  [seriesinfo.SeriesNumber];
end;
if iscell( seriesnrs )
    seriesdescription = seriesnrs;
    [seriesinfo] = loaddicom([sourcefilesdir '.series'],1,1); 
    if opts.regexpMatchSeries
        seriesnrs = [];
        sel = false(numel(seriesinfo),1);
        for k=1:numel(seriesdescription) 
            if isequal(seriesdescription{k},'Ask user')
                S = {seriesinfo.SeriesDescription}';
                for si = 1:numel(S)
                    S{si} = sprintf('%03d : %s',seriesinfo(si).SeriesNumber,S{si});
                end;
                [selupd, ok] =listdlg('ListString', S,'InitialValue', find(sel) , 'PromptString','Please update the selection of dicom series that will be loaded.','ListSize',[350 300]);
                if ok
                    % overwrite entire selection if updating
                    sel(:)=false;
                    sel(selupd)=true;
                    seriesnrs = [];
                else
                    % no additional series:
                    error('Cancel pressed. Not loading data');
                    selupd = [];
                end;
            else
                q = regexp({seriesinfo.SeriesDescription}',seriesdescription{k});
                selupd = ~sel & ~cellfun(@isempty,q);
                sel = sel | selupd;
            end;
            seriesnrs = [seriesnrs seriesinfo(selupd).SeriesNumber]; %#ok<AGROW> : not performance critical and difficult to avoid.
        end;
        if isempty(seriesnrs )
            warning('no series match the provided regular expression(s); no series are loaded');
            MRimgs=[];
            info ={};
            T= {}; 
            objT = eye(4);
            inputs= [];
            relativeSliceWidth = [];
            szobj = [];
            return;
        end;
    else
        seriesnrs = zeros(1,numel(seriesdescription));
        for k=1:numel(seriesdescription)
            lidx = find(strcmp({seriesinfo.SeriesDescription},seriesdescription{k}));
            if numel(lidx)~=1
                error(['cant uniquely find series ' seriesdescription{k} ' in dataset : ' sourcefilesdir ]);
            end;
            seriesnrs(k) = seriesinfo(lidx).SeriesNumber;
        end;
    end;
end;
MRimgs = cell( 1, numel(seriesnrs) );
info = cell( size(MRimgs) );
T = cell( size(MRimgs) );

szMRI = ones( numel(seriesnrs), ndatadim );
origszMRI = cell(numel(seriesnrs),1);
MRIdownsamplefact = ones( numel(seriesnrs), ndatadim );
voxelspacing = zeros( numel(seriesnrs), ndim );
imgCenterWorldSpace = zeros(numel(seriesnrs), ndim);
progressbar('adjustlimits',[-1 numel(seriesnrs)]);
progressbar(0);
for serienr = 1:numel(seriesnrs)
    tst = seriesnrs == seriesnrs(serienr);
    if nnz(tst)==1
        seriestr = num2str(seriesnrs(serienr));
    else
        seriestr = [num2str(seriesnrs(serienr)) ',' num2str(find(find(tst)==serienr))];
    end;
    [data, info{serienr}] = loaddicom([sourcefilesdir '.' seriestr],1,1);
    data = squeeze(data);
    origszMRI{serienr} = size(data);
    T{serienr} = affineTransformFromDicominfo(info{serienr});
%     T{serienr}(1:2,1:2) = T{serienr}([2 1],[2 1]);
    if ~isempty(opts.imagePermute) && ~isempty(opts.imagePermute{serienr})
        data = permute(data, opts.imagePermute{serienr});
        T{serienr}(1:numel(opts.imagePermute{serienr}),:) = T{serienr}(opts.imagePermute{serienr},:);
    end;
    MRimgs{serienr} = data;
    szCurImg = size(data);
    if (numel(szCurImg) > ndatadim) 
        ndatadim = numel(szCurImg);
        % expand matrices that structurally depend on number of dimensions:
        szMRI(:,end+1:ndatadim) = 1;
        MRIdownsamplefact(:,end+1:ndatadim)=1;
%         else
%             error('images differ in number of dimensions.');
%         end;
    end;
        
    if opts.resample_to_acquisitiongrid
        if ~isempty(opts.AcquisitionMatrix) && ~isempty(opts.AcquisitionMatrix{serienr})
            info{serienr}.AcquisitionMatrix = opts.AcquisitionMatrix{serienr};
        end;
        if isfield(info{serienr},'AcquisitionMatrix')
        % downsample images to original (=acquisitionmatrix) grid
        % Assume upsampling in slice direction when sliceThickness is twice the spacing between the slices.
        Acqmat1 = info{serienr}.AcquisitionMatrix([2 3]); % verified to be correct for: 'd:\Datasets\UZA\Volunteer\'
        Acqmat2 = info{serienr}.AcquisitionMatrix([4 1]); % verified to be correct for: 'd:\Datasets\UZA\Volunteer\'
        if all(Acqmat1==0)
            Acqmat = Acqmat2;
        elseif all(Acqmat2==0)
            Acqmat = Acqmat1;
        else
            error('I can''t (yet) interpret this acquisition matrix.');
        end;
        
        resamplesz = min( szCurImg , [ceil(double(Acqmat')*1.07+.9) szCurImg(3:end)]);
        if isfield(info{serienr},'SliceThickness')
            if isfield(info{serienr},'MRAcquisitionType') && isequal(info{serienr}.MRAcquisitionType,'3D')
                resamplesz(3) = min( szCurImg(3), ceil( szCurImg(3) * info{serienr}.SpacingBetweenSlices/ info{serienr}.SliceThickness *1.07+.9));
            else
                if abs(info{serienr}.SliceThickness/info{serienr}.SpacingBetweenSlices-2)<.01
                    % assume 2 fold zero-fill interpolation in slice direction:
                    resamplesz(3) = resamplesz(3)/2;
                end;
            end;
        end;
        MRIdownsamplefact(serienr, 1:numel(szCurImg) ) = szCurImg./resamplesz;
        if ~all(MRIdownsamplefact(serienr,:)==1)
            if opts.displayFFTbeforeDownsampling
                % display FFT of image, to inspect frequency content:
                wndw = cell(1,3);
                timg = double(MRimgs{serienr}(:,:,:,1));
                for k=1:3; 
                    wndw{k}=normpdf(linspace(-2.5,2.5,size(MRimgs{serienr},k)))';
                    timg = bsxfun(@times, permute(wndw{k},[2:k 1 k+1]), timg); 
                end;
                imagebrowse(abs(fftn(timg)));    
                resamplesz./szCurImg
            end;
            firstvoxloc = [1 1 1 1] * T{serienr};
            scale = szCurImg./resamplesz;
            T{serienr} = diag([scale(1:3) 1])* T{serienr};
            % compensate shift due to first voxel at [1 1 1], while fftshift preserves
            % same point at first voxel:
            T{serienr}(4,:) = firstvoxloc - [1 1 1 0] * T{serienr};
            MRimgs{serienr} =  fftresample(double(MRimgs{serienr}),resamplesz);
        end;
        else
            disp(['Series '  num2str(seriesnrs(serienr)) ' lacks field AcquisitionMatrix. Cannot resample image to acquisition grid']);
        end;
    end;
    szCurImg = size(MRimgs{serienr});
    szMRI(serienr,1:numel(szCurImg)) = szCurImg;
    imgCenterWorldSpace(serienr,:) = [(szMRI(serienr,1:ndim)+1)/2 1] * T{serienr}(:,1:ndim);
    voxelspacing(serienr,:) = sqrt(sum(T{serienr}(1:ndim,:).^2,2))';
    progressbar(serienr);
end;
progressbar('ready');
clear data
if isempty(opts.objectVoxelSize) || isequal(opts.objectVoxelSize,0)
    objvoxel_size = min(voxelspacing(:));
else
    objvoxel_size = opts.objectVoxelSize;
end;
if isempty(opts.objectTransform)
    objT = diag( [ ones(1,ndim) .* objvoxel_size 1] );
    % move mean image center to center of reconstructed object:
    szobj = opts.objectSize;
    if isempty(szobj)
        worldbox = [inf;-inf]*ones(1,ndim);
        for k=1:numel(T)
            minpo = zeros(1,ndim);
            for k2 = 1:ndim
                stvec = [szMRI(k,1:3) 1];
                stvec( T{k}(:,k2)>0 ) = 1;
                minpo(k2) = stvec*T{k}(:,k2);
            end;
            worldbox(1,:) = min(worldbox(1,:), minpo);

            maxpo = zeros(1,ndim);
            for k2 = 1:ndim
                stvec = [szMRI(k,1:3) 1];
                stvec( T{k}(:,k2)<0 ) = 1;
                maxpo(k2) = stvec*T{k}(:,k2);
            end;
            worldbox(2,:) = max(worldbox(2,:), maxpo);
        end;
        szobj = diff([worldbox [1;1]]/objT ,1,1);
        szobj = ceil(abs(szobj(1:ndim)));
    end;
    imgoffset = [(szobj+1)/2 1]*objT - [mean(imgCenterWorldSpace,1) 1];
    objT(end,:) = objT(end,:) - imgoffset;
else
    % objectTransform specified.
    objT = opts.objectTransform;
    szobj = opts.objectSize;
    if numel(opts.objectTransform)==1 
        % objectTransform specifies volume indec to use as object orientation.
        origobjT = affineTransformFromDicominfo(info{opts.objectTransform});
        if isempty(opts.objectVoxelSize)
            objT = origobjT;
            voxsc = ones(size(objT,1)-1,1);
        else
            voxsc = objvoxel_size(:)./sqrt(sum(origobjT(1:end-1,:).^2,2));
            objT = diag([voxsc; 1])*origobjT;
        end;
        % read size from originitating object
        origszobj = origszMRI{opts.objectTransform};
        origszobj = origszobj(1:min(3,numel(origszobj)));
        if isempty(szobj)
            szobj = ceil(origszobj./voxsc');
        end;
        % shift objT so it centers on the original image:
        origcenter = [(origszobj+1)/2 1]*origobjT;
        objT(end,:) = origcenter -[(szobj+1)/2 0]*objT;
    else
        % explicit object transform specified.
        objT = opts.objectTransform;
        if isempty(szobj)
            % no size specified. 
            % Take size so that everything is contained (without changing
            % object origin)
            szobj = ones(1,ndim);
            for k=1:numel(T)
                maxpo = zeros(1,ndim);
                Tk = T{k}/objT;
                for k2 = 1:ndim
                    stvec = [szMRI(k,1:3) 1];
                    stvec( Tk(:,k2)<0 ) = 1;
                    maxpo(k2) = stvec*Tk(:,k2);
                end;
                szobj = max(szobj, maxpo);
            end;
            szobj = ceil( szobj );
        end;
    end;
end;
inputs = zeros(ndim+1,ndim+1,numel(T));
for k=1:numel(T)
    inputs(:,:,k) = T{k}/objT;
end;
inputs = inputs(:,1:ndim,:);
relativeSliceWidth  = ones(1,numel(MRimgs));
for k=1:numel(MRimgs)
    if isfield(info{k},'SliceThickness')
        SpacingBetweenSlices = norm(T{k}(3,:));
        if isfield(info{k},'SpacingBetweenSlices')
            if abs(SpacingBetweenSlices/abs(info{k}.SpacingBetweenSlices*MRIdownsamplefact(k,3))-1)>1e-4
                error('SpacingBetweenSlices of transform differs to much from dicominfo.SpacingBetweenSlices');
            end;
        end;
        relativeSliceWidth(k) = info{k}.SliceThickness/ SpacingBetweenSlices;
    else
        disp(['Series '  num2str(seriesnrs(k)) ' lacks field SliceThickness. Cannot compute relative slice width, defaulting to 1']);
    end;
    if isfield(info{k},'MRAcquisitionType')
        if isequal(info{k}.MRAcquisitionType,'3D')
            relativeSliceWidth(k) = -relativeSliceWidth(k);
        elseif isequal(info{k}.MRAcquisitionType,'2D')
            % if known to be 2D then do nothing.
        elseif isfield(info{k},'ProtocolName') && ~isempty( strfind( info{k}.ProtocolName , '3D') ) 
            relativeSliceWidth(k) = -relativeSliceWidth(k);
            disp(['Series '  num2str(seriesnrs(k)) ' : unknown value for MRAcquisitionType. Protocolname contains 3D, so assuming 3D acquisition.']);
        elseif ~isequal(info{k}.MRAcquisitionType,'2D')
            disp(['Series '  num2str(seriesnrs(k)) ' : unknown value for MRAcquisitionType, assuming 2D']);
        end;
    else
        disp(['Series '  num2str(seriesnrs(k)) ' does not specify type of acquisition (2D or 3D), assuming 2D for relativeSliceWidth computation']);
    end;
end;