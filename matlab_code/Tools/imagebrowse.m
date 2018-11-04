function [fighnld] = imagebrowse(IMG, doscale, varargin)
% imagebrowse(IMG [, doscale/minmax , option_i, value_i,... ])
% depricated: imagebrowse(IMG [, doscale [, xrange]])
%
% Shows the N dimensional image IMG, in a browsable view of slices. Click
% on a point to focus there, or use the sliders to navigate through the
% image. The standard zoom buttons can also be used to view something in
% more detail.
% Complex data can be visualized as well. Use the radiobuttons to select
% between views of the magnitude, angle, real or imag part of the data.
% A log checkbox is present to easily allow a view of the log of the data.
%
% INPUTS:
%  IMG:  N dimensional image, any numeric datatype, real or complex.
%  doscale : default = 1; An optional scalar integer argument that
%            specifies the kind of automatic scaling that is performed:
%       0: no scaling, uint8 inputs are treated as is, and other inputs are
%          assumed to have a range of [0..1].
%       1: (default) scale so the maximum of the image is white and the minimum
%          is black.
%       2: scale each slice so that it contains black and white.
%  minmax : a 2 element vector that contains the custom black and white points (respectively).
%           A 4 x 2 minmax matrix can provide a limit for each of the 
%           magnitude, angle, real, imag part of img (respectively)
%
% Supported options:
%  colormap : default = gray(256), any colormap for the data (if not full
%             color)
%  labels   : default = {'Dimension ' i , } ; a N element cell array with string
%               values containing a (meaningfull) description of the
%               dimensions.
%  ROIs     : A cell array with ROI's. Each element is a m x N matrix that
%             specify the m edge points of the ROI (/curve). For correct
%             viewing, each ROI should be non-constant in (exactly) 2 dimensions 
%             nan values indicate all planes in that dimension.
%            (2 columns). 
%             Alternatively 'polygon' objects (by D.H.J. Poot)
%             can be drawn.   
%  ROIcolors: numel(ROIs) element cell array with a color for each ROI.
%  xrange   : The range of the axis can be specified by xrange. Specify the
%             location of the first and last element in IMG in each of the
%             respective dimensions. 
%           So xrange = [location_IMG(   1 ,:,..,:) location_IMG(:,   1 ,..) .. ;
%                        location_IMG( end ,:,..,:) location_IMG(:, end ,..) .. ];
%             Default: xrange = [ones(1, ndims);size(IMG) ];
%             Example: If you have a pixelspacing and the first voxel should be at 0
%             use xrange = [zeros(1,ndims(IMG));
%                           pixelspacing.*(size(IMG)-1)]
%  colordim : If non empty, a scalar integer that specifies a dimension in
%             which size(IMG, colordim)==3.
%             Default: the last dimension in which size(IMG)==3, if any such is present.
%  description : a string with a description of the content in this figure.
%  notes    : 1 x N cell array with in element 'i' is an   size(IMG,i)
%             element cell array with strings. Each string may contain a description for
%             the corresponding index in dimension i.
%
% Created by Dirk Poot, 6-12-2007, University of Antwerp, Erasmus MC, TUDelft.

% 27-1-2009: Added complex data support 
% 25-8-2009: Added center and range sliders. 
% 15-12-2011: Added scroll weel support. Unfortunately I can't enable it yet, as 
%             selecting 'zoom' overwrites the 'WindowScrollWheelFcn' option of the figure.
%             And when 'zoom' is activated, it cannot be set at all.
%    -----   : Enabled scroll weel support.
% 5-9-2012 : D. Poot: added 'option-value' support; added ROI and colormap support 
% 22-11-2013: D.Poot: made all subfigures zoom at the same time. (can be
%             disabled by setting option 'linkzoom' to false)
% 16-04-2014: D.Poot: added overlays and cleaned up some of the internal   
%             structure and data handling. 

% Implementation documentation:
% Internal structure and data storage:
% upon construction of the figure the info structure is expanded. 
% It contains all handles and input options/values.
% This info structure is stored in the 'userdata' of the figure
% Since this structure may be large (it contains IMG) it might be horribly
% slow to update it (=a bug in at least some versions of MATLAB) 
% Therefore all dynamically updated data is stored in a separate structure,
% dinfo, which is stored in the 'userdata' field of the panel that contains the
% subfigures. Since that contains just a few small fields that should stay
% fast....


info.colormap  = gray(256);
info.ROIs = {};
info.ROIcolors = [];
info.xrange = [];
info.labels = {};
info.colordim = 'lastpossible'; 
info.neutralColor = 0; % neutral color used only for creating the colorbar
info.description = [];
info.notes = [];
info.linkzoom = true; % link zooming across all subfigures. Before 22-11-2013 this was not supported; Old behaviour is linkzoom = false
info.overlays = [];  % overlays, structure with fields:
                     %  mask   : actually opacity map; 0 = fully transparant, 1 = fully opaque, same size as image.
                     %  image  : overlay image with size(image,i)==1 or size(IMG,i
                     %  colormap : colormap used by this overlay when image is not true color.
                     %  min    : value mapping to first entry in colormap
                     %  max    : value mapping to last entry in colormap.
                     %  name   : name of this overlay (shown in the
                     %           overlay select - checkboxes)
dinfo = struct;
dinfo.sliderpos = [];

sIMG = size(IMG);
ndims = length(sIMG);

if nargin<2 || isempty(doscale)
    doscale =true;
end;

% parse option value pairs (and handle deprecated calling order)
if nargin>=3 
    if nargin==3 && isnumeric(varargin{1}) && isequal( size(varargin{1}), [2 ndims])
        info.xrange = varargin{1};
        warning('IMAGEBROWSE:depricated', 'calling imagebrowse with third argument xrange is depricated, use option value pairs instead');
    else
        if ~exist('parse_defaults_optionvaluepairs', 'file')
            warning('IMAGEBROWSE:optionparsefunction_notfound', 'Function "parse_defaults_optionvaluepairs" not found. Option-value pair arguments are ignored, please get "parse_defaults_optionvaluepairs" from D.H.J. Poot.');
        else
            info = parse_defaults_optionvaluepairs(  info, varargin{:} );
        end;
    end;
end;


if issparse(IMG)
    if numel(IMG)-nnz(IMG)>1e7
        warning('imagebrowse:SparseIMG',['Many zero elements in sparse image (' num2str(numel(IMG)-nnz(IMG)) '). Therefore, the image is not converted to full and thus cannot be shown. Figure is not modified.']);
        return;
    end;
    IMG = full(IMG); % only case in which we 'copy' the input.
end;
if isreal(doscale) && (length(doscale)==2 || isequal(size(doscale),[4 2]))
    if size(doscale,1)==4
        mn = doscale(:,1);
        mx = doscale(:,2);
    else
        mn = doscale(1);
        mx = doscale(2);
    end;
    doscale = true;
elseif doscale
    if exist('advancedminmax','file')
        [mn,mx] = advancedminmax(IMG(:));
    else
        disp('Try to get advancedminmax (from D.H.J. Poot)!');
        mn = min(IMG(:));
        mx = max(IMG(:));
    end;
    if any(mn==mx)
        mx = mn+1;
    end;
elseif ~isa(IMG,'uint8')
    mn = 0;
    mx = 1;
%     doscale = 1;
else
    mn = 0;
    mx = 255;
end;

tileDims = 1:ndims;
rngIMG = [ones(1,ndims);sIMG];
if isequal( info.colordim , 'lastpossible');
    colordim = find(sIMG==3,1,'last');
else
    if numel(info.colordim)~=1 || sIMG( info.colordim )~=3
        if ~isempty(info.colordim)
            warning('Can only select a single "color" dimension and the size of that dimension should be 3 (RGB)');
        end;
        colordim = [];
    else
        colordim = info.colordim;
    end;
end;

if ~isempty(info.xrange)
    if ~isequal(size(info.xrange),[2 ndims])
        error('image range should be a (2 x ndims) matrix specifying coordinate of [first;last] element.');
    end;
    rngIMG = info.xrange;
end;
if ~isempty(colordim)
    tileDims(colordim) =[]; %remove color dim from tiling.
    ndims = ndims-1;
end;
squeezedims = sIMG( tileDims )==1;
if any(squeezedims)
    tileDims(squeezedims)=[]; % remove size =1 dimensions.
    ndims = ndims - nnz(squeezedims);
end;

% Test if we can preserve the current imagebrowse figure and just update
% the image (preserves settings such as log, position, range, center)
if isvalidImageBrowsefigure(get(0,'currentfigure'))
    fighnld = get(0,'currentfigure');
    ud = get(fighnld,'UserData');
    if isequal(ud.sIMG,size(IMG)) ... % size, range should be equal
         && isequal(ud.rngIMG, rngIMG) ... % limits of axis should be equal
         && isequal(ud.minmax(:,:,1), [double(mn) double(mx)])
        ud.IMG = IMG;
        set(fighnld,'UserData',ud);
        redrawsubplots(fighnld, 1:numel(ud.subfigs));
        return;
    end;
end;

%  Start constructing the figure:
clf
fighnld = gcf;
% rngIMG = [ones(1,ndims);sIMG];
info.fighnld = fighnld;
info.IMG = IMG;
info.sIMG = sIMG;
info.rngIMG = rngIMG;
info.sliders = zeros(ndims,1);
info.texts = zeros(ndims,1);
info.tileDims =  tileDims;
dinfo.drawmethod = 0; 
usedspace = 0;
% construct the position sliders (and accompagnying texts):
for k=1:ndims
    info.sliders(k) = uicontrol('Style','slider','Min',1,'Max',sIMG(tileDims(k)),...
                                'value'      ,  ceil(sIMG(tileDims(k))/2),       ...
                                'SliderStep' , [1/(sIMG(tileDims(k))-1) round(sqrt(sIMG(tileDims(k))))/(sIMG(tileDims(k))-1)],...
                                'callback',{@sliderupdate,[info.fighnld , k]},...
                                'units','pixels',...
                                'Position',[5 (ndims-k)*50+10 100 20]);

    msd = floor(log10(max(abs(rngIMG(:,tileDims(k))))));
    lsd = floor(log10(diff(rngIMG(:,tileDims(k)))/(sIMG(tileDims(k))-1)));
    if lsd<0
        posstr = ['[%' num2str(msd-lsd+2) '.' num2str(-lsd) 'f]'];
    elseif lsd<=2
        posstr = ['[%' num2str(msd+1) '.0f]'];
    end;
    info.texts(k) = uicontrol('Style','text',...
                                'units','pixels', ...
                                'Position',[5 (ndims-k)*50+30 100 18], ...
                                'UserData',['Dim ' num2str(info.tileDims(k)) ', ' posstr]);
end;
usedspace = usedspace+ ndims*50;
% Construct center and range sliders (and accompagnying texts):
if  doscale~=2
    usedspace = usedspace+10;

    info.sliders(ndims+2) = uicontrol('Style','slider','Min',-2,'Max',10,...
                                'value'      ,  0,       ...
                                'SliderStep' , [1/(120-1) 1/(12-1)],...
                                'callback',{@sliderupdate,[info.fighnld , k]},...
                                'units','pixels',...
                                'Position',[5 usedspace+(1)*50+10 100 20]);
    info.texts(ndims+2) = uicontrol('Style','text',...
                                'units','pixels', ...
                                'Position',[5 usedspace+(1)*50+30 100 18], ...
                                'UserData','Range: %1.4g');
    info.sliders(ndims+1) = uicontrol('Style','slider','Min',mn(1),'Max',mx(1),...
                                'value'      ,  (mx(1)+mn(1))/2,       ...
                                'SliderStep' , [1/999 1/30],...
                                'callback',{@sliderupdate,[info.fighnld , k]},...
                                'units','pixels',...
                                'Position',[5 usedspace+(0)*50+10 100 20]);
    info.texts(ndims+1) = uicontrol('Style','text',...
                                'units','pixels', ...
                                'Position',[5 usedspace+(0)*50+30 100 18], ...
                                'UserData','Center:  %1.4g');
    usedspace = usedspace+ 2*50;
end;


butsz = 20;
if ~isreal(IMG)
    % Create magnitude/real/imag/phase radio buttons:
    
    % Create four radio buttons in the button group.
    h = uibuttongroup('visible','off','units','pixels','Position',[5 usedspace+10 100 butsz*4+4]);
    usedspace = usedspace+ butsz*4+10+4; 
    u0 = uicontrol('Style','Radio','String','Magnitude',...
        'pos',[2 3*butsz+2 96 butsz],'parent',h,'Tag','M');
%    u1 = ...
         uicontrol('Style','Radio','String','Phase',...
        'pos',[2 2*butsz+2 96 butsz],'parent',h,'Tag','P');
%     u2 = ...
         uicontrol('Style','Radio','String','Real',...
        'pos',[2 1*butsz+2 96 butsz],'parent',h,'Tag','R');
%     u3 = ...
         uicontrol('Style','Radio','String','Imag',...
        'pos',[2 0*butsz+2 96 butsz],'parent',h,'Tag','I');
    % Initialize some button group properties. 
    set(h,'SelectionChangeFcn',{@selcbk,info.fighnld});
    set(h,'SelectedObject',u0);  % No selection
    set(h,'Visible','on');
    dinfo.drawmethod = 'M'; 
end;
% Create 'log' box:
dolog = 0;
checkloghndl = uicontrol('Style','checkbox','String','Log',...
        'pos',[7 usedspace+10 100 butsz],'Callback', {@togglelog,info.fighnld});
usedspace = usedspace+butsz+10;
buttonhndl = uicontrol('Style','pushbutton','String','100%',...
        'pos',[7 usedspace+10 50 butsz],'Callback', {@setzoom,info.fighnld});
usedspace = usedspace+butsz+10;

% add overlay checkboxes:
if ~isempty( info.overlays )
    if ~isfield(info.overlays,'name')
        info.overlays(1).name = '';
    end;
    if ~isfield(info.overlays,'colormap')
        info.overlays(1).colormap =[];
    end;
    h = uibuttongroup('visible','off','units','pixels','Position',[5 usedspace+10 100 butsz*numel(info.overlays)+4]);
    usedspace = usedspace+ butsz*numel(info.overlays)+10+4; 
    overlayboxes = zeros(1,numel(info.overlays));
    for k=1:numel(info.overlays)
        if isempty(info.overlays(k).name)
            info.overlays(k).name = ['Overlay ' num2str(k)];
        end;
        if isempty(info.overlays(k).colormap)
            info.overlays(k).colormap = jet;
        end;
        overlayboxes(k) = uicontrol('Style','checkbox','String',info.overlays(k).name,'pos',[2 (numel(info.overlays)-k)*butsz+2 96 butsz],'parent',h,'Value',1,'Callback',@(obj, dummy) checkboxchange( obj, dummy, k, h) );
        dinfo.overlays(k).min = info.overlays(k).min;% current plot range = default range.
        dinfo.overlays(k).max = info.overlays(k).max;% current plot range = default range.
    end;
    set(h,'Visible','on');
%      info.overlayboxes = overlayboxes;
    info.overlayselectors = h; 
end;
dinfo.activeOverlays = true(1,numel(info.overlays));

% Add colorbar:
colorbarheight = 130;
info.colorbar = axes('Units','pixels','Position',[10 usedspace+30 20 colorbarheight],'parent',info.fighnld);
usedspace = usedspace+30+colorbarheight;
info.colorbarname =  uicontrol('Style','text','units','pixels', 'Position',[5 usedspace+5 100 18]);
if numel(colordim)>0
    cmsk = [1 1 1; 0 0 1; 0 1 0; 1 0 0];
    colorbarimagemat = {reshape(((0:255)'/255-info.neutralColor)*kron(cmsk(:)',ones(1,10)) , [256,size(cmsk,1)*10, 3])+info.neutralColor};
else
    colorbarimagemat = {(0:size(info.colormap,1)-1)'*ones(1,10)};
end;
colorbarimagemat{2} = linspace(0,1,256)'*ones(1,10);
info.colorbarimage = image([0 1],[0 1], colorbarimagemat{1} );
info.colorbarimagemats = colorbarimagemat ;
set(info.colorbar,'xTick',[],'YDir','normal','YAxisLocation','right');
% update_colorbar( info.fighnld, 0); called after creating panel and
% setting info and dinfo
colormap(info.colormap);

% Create the panel containing the subfigures:
info.panel = uipanel('Position',[.2 0 .8 1]); % position is overwritten by call to resizefig
% Create all subfigures:
totalnumsubfigs = ndims*(ndims-1)/2;
info.subfigs = zeros(totalnumsubfigs,1);                 
doaddnotes_box = ~isempty( info.notes ) || ~isempty(info.description);
if doaddnotes_box
    totalnumsubfigs = totalnumsubfigs + 1;
end;
info.imgs = zeros(size(info.subfigs,1),1);                 
info.subfigdims = zeros(size(info.subfigs,1),2+numel(colordim));
info.doscale = doscale;
dinfo.dolog = dolog;
if numel(mn)>1
    midx = 2;
    if mx(midx)>4/5*pi && mn(midx)<-4/5*pi
        mn(midx) = -pi;
        mx(midx) = pi;
    end;
    [mn_lg,mx_lg]=deal([log(mn(1)) -pi log(mn(1)) -pi]',[log(mx(1)) pi log(mx(1)) pi]');
    mn_lg = max(mn_lg, mx_lg-30);
else
    if mx<0
        mn_lg = log(-mx);
        mx_lg = log(-mn);
    else
        mn_lg = log(max(eps*abs(mx),mn));
        mx_lg = log(abs(mx));
    end;
    mx_lg = real(mx_lg);
    mn_lg = max(mx_lg-30,real(mn_lg));
end;
info.minmax = cat(3,[double(mn) double(mx)], [double(mn_lg) double(mx_lg)]);
dinfo.minmax = info.minmax;  % current plot range = default range.
nvertsubfig = round(sqrt(totalnumsubfigs-.1));
nhorisubfig = ceil(totalnumsubfigs/nvertsubfig);
m=0;
rngAx = rngIMG + [-.5;.5]*(diff(rngIMG)./(sIMG-1));
if isempty(info.labels)
    info.labels = cell(1, numel(sIMG) );
end;
for k = 1 : numel(sIMG)
    if isempty(info.labels{k}) && ~ischar(info.labels{k})
        % if label not specified provide default. Allow specification of
        % empty label by ''.
        info.labels{k} = ['Dimension ' num2str(k)];
    end;
end;
for k=1:ndims-1
    for l=k+1:ndims
        hpos = mod(m,nhorisubfig);
        vpos = floor(m/nhorisubfig);
        m=m+1;
        info.subfigs(m) = axes('OuterPosition',[hpos/nhorisubfig 1-(vpos+1)/nvertsubfig 1/nhorisubfig 1/nvertsubfig],'parent',info.panel);
        set(info.subfigs(m),'xlim',rngAx(:,tileDims(l))','ylim',rngAx(:,tileDims(k)'));
        zoom reset;
        info.imgs(m) = image(rngIMG(:,tileDims(l)),rngIMG(:,tileDims(k)), zeros(sIMG(tileDims(k)),sIMG(tileDims(l)),'uint8'));
        colormap(info.colormap); % some version of MATLAB uses different colormap per axes.
        set(info.imgs(m),'ButtonDownFcn',{@imageclick,{gcf, m}});
        xlabel(info.subfigs(m), info.labels{ tileDims(l) });
        ylabel(info.subfigs(m), info.labels{ tileDims(k) });
        info.subfigdims(m,:) = [tileDims(k) tileDims(l) colordim];
    end;
end;

% Create notes box. Content filled in updatesliders:
info.notes_box = [];
if doaddnotes_box
    hpos = mod(m,nhorisubfig);
    vpos = floor(m/nhorisubfig);
    borderspace_box = .05;
    info.notes_box = uicontrol('Style','edit','Min',1,'Max',10,...
                               'String','',...
                               'BackgroundColor',[1 1 1],...
                               'HorizontalAlignment','left',...
                               'Units', 'normalized',...
                                'Position',[(hpos+borderspace_box)/nhorisubfig 1-(vpos+1-borderspace_box)/nvertsubfig 1-(hpos+2*borderspace_box)/nhorisubfig (1-2*borderspace_box)/nvertsubfig],...
                                'parent',info.panel);
end;
dinfo.activeslider = info.sliders(1);

if ~isempty( info.ROIs ) 
    info = prepareROIs( info );
end;
% all objects created. 
% Create callbacks on the figure and store the userdata : 

set(info.fighnld,'ResizeFcn',{@resizefig,info.fighnld});
set(info.fighnld,'userdata',info);
set(info.panel,'userdata',dinfo);
% Final updates:
update_colorbar( info.fighnld, 0);
resizefig(info.fighnld,[],info.fighnld);
sliderupdate(info.fighnld,[],[info.fighnld 0]);
set(info.fighnld,'Toolbar','figure');
set(info.fighnld,'WindowButtonMotionFcn',{@sliderupdate,info.fighnld})

    % unfortunately, the zoom, pan, ... modes do prevent setting WindowScrollWheelFcn
    % and return an warning that is not interesting for someone calling imagebrowse. 
    % Therefore, turn that warning off. When applicable, update state that gets restored when the zoom, pan ... mode exits.
    % preserve old warning state.
    [oldwarn, oldmsgid] = lastwarn('','');
    ws = warning('off','MATLAB:modes:mode:InvalidPropertySet');
    set(info.fighnld,'WindowScrollWheelFcn',{@scrollweelupdate,info.fighnld});
    warning(ws);
    [tst,tstid] = lastwarn(oldwarn, oldmsgid);
    if isequal(tstid, 'MATLAB:modes:mode:InvalidPropertySet')
        % apparently, some ui mode was selected.
        % set WindowScrollWheelFcn in state that gets restored when ui-mode ends:
        hManager = uigetmodemanager(info.fighnld);
        currMode = get(hManager,'CurrentMode');
        currMode.FigureState.WindowScrollWheelFcn = {@scrollweelupdate,info.fighnld};
%         warning('IMAGEBROWSE:cantSetWindowScrollWheelFcn', 'Please deselect all ui modes (zoom, pan, ...) when starting imagebrowse, since now scrolling through the data is disabled.');
    end;
    
if info.linkzoom
    zoomhndl = zoom;
    set(zoomhndl, 'ActionPostCallback', @afterzoomupdate)
end;


% Figure now created. Return to caller, and wait for some callback.....

function set_minmax( info, mn, mx)
set(info.colorbarimage,'YData',[mn mx]);
set(info.colorbar,'Ylim',[mn mx]);
dinfo = get(info.panel,'userdata');
imgid = dinfo.activeColormapImage;
if imgid==0
    midx = 1;
    switch dinfo.drawmethod % reduce complex data:
        case 'P'
            midx = 2;
        case 'R'
            midx = 3;
        case 'I'
            midx = 4;
    end;
    dinfo.minmax( midx, 1:2, dinfo.dolog + 1 ) = [mn mx];
else
    dinfo.overlays(imgid).min = mn;
    dinfo.overlays(imgid).max = mx;
end;
set(info.panel,'userdata', dinfo);
       
function sliderupdate(src,event,infosrc)%#ok: event
% Called when any (can be multiple) of the sliders is updated. 
% Adjust the texts above the sliders, possibly the notes_box, and redraw
% the subfigures (that may have been changed). 
fighndl = infosrc(1);
info = get(fighndl,'userdata');
if ~ishandle(info.panel)
    stopimagebrowse(fighndl); % happens when clf is used.
    return;
end;
dinfo = get(info.panel,'userdata');
sliderpos = zeros(size(info.sliders));
for k=1:numel(sliderpos);
    if ~ishandle(info.sliders(k))
        stopimagebrowse(fighndl); % happens when clf is used.
        return;
    end;
    sliderpos(k) = get(info.sliders(k),'value');
    
end;
if isequal(get(src,'Type'),'uicontrol') && isequal(get(src,'Style'),'slider')
    if dinfo.activeslider~=src
        dinfo.activeslider=src;
        set(info.panel,'userdata',dinfo);
    end;
end;
sliderpos(1:numel(info.tileDims)) = round(sliderpos(1:numel(info.tileDims)));
oldsliderpos = dinfo.sliderpos;
if isequal(oldsliderpos,sliderpos)
    return;
end;
dinfo.sliderpos = sliderpos;
set(info.panel,'userdata',dinfo);
sliderchanged = false(1,numel(sliderpos));
if numel(sliderpos)==numel(oldsliderpos)
    sliderchanged(sliderpos~=oldsliderpos) = true;
else
    sliderchanged(1:end+1) = true;
end;
% baseindx = repmat({1},numel(info.sIMG),1);
for k=1:numel(sliderpos);
%     if sliderchanged(k)
        if k==numel(info.tileDims)+2
            [mn, mx] = getminmax(info);
            adj = (mx-mn)*2^(-1-sliderpos(end));
            newmn = sliderpos(end-1)-adj;
            newmx = sliderpos(end-1)+adj;
            set_minmax( info, newmn, newmx);
            sliderval = newmx-newmn;
%             [sliderval, mn, mx] = getcontrastval(info);
%             adj = (mx-mn)*2^(-1-sliderpos(end));
%             [mn ,mx] = deal( sliderpos(end-1)-adj,sliderpos(end-1)+adj);
%             colrval = [1 0;.75 .25;.5 .5; .25 .75;0 1]*[mn;mx];
%             set(info.colorbar,'YtickLabel',num2str(colrval))
        else
            sliderval = sliderpos(k);
        end;
        if k<=numel(info.tileDims)
            rng = info.rngIMG(:,info.tileDims(k));
            sliderval = (sliderval-1)/(info.sIMG(info.tileDims(k))-1) * (rng(2)-rng(1))+rng(1);
        end;
        set(info.texts(k), 'String' ,  sprintf(get(info.texts(k),'UserData'),sliderval)); 
%     end;
%     baseindx(info.tileDims(k)) = {sliderpos(k)};
end;
if ~isempty(info.notes_box) && any(sliderchanged) 
    % update note box:
    baseidx = sliderpos(1:numel(info.tileDims));
    strs = cell( 2, numel(baseidx));
    lastcolidx = 0;
    if ~isempty(info.notes)
        for k=1:numel(baseidx);
            if ~isempty(info.notes{k})
                lastcolidx = lastcolidx + 1;
                strs{1,lastcolidx} = info.labels{ k };
                strs{2,lastcolidx} = info.notes{k}{ baseidx(k) };
            end;
        end;
    end;
    if isempty( info.description )
        note_str = sprintf(repmat('%s:\n  %s\n',1,lastcolidx),strs{:,1:lastcolidx});
    else
        note_str = sprintf(['%s\n' repmat('\n%s:\n  %s',1,lastcolidx)],info.description,strs{:,1:lastcolidx});
    end;
    set(info.notes_box,'String',note_str);
end;
% info.baseindx = baseindx;
% for efficiency do not redraw those that do not change.
sliderchangedd = false(1,numel(info.sIMG));
sliderchangedd(info.tileDims) = sliderchangedd(1:numel(info.tileDims));
updatelist = find(sum(sliderchanged) ~= sum(sliderchangedd(info.subfigdims),2)');
redrawsubplots(fighndl, updatelist)

function scrollweelupdate(src,event,infosrc)
% Called when using the scroll wheel (of the mouse)
% Update the 'active' slider.
if event.VerticalScrollCount ~= 0 
fighndl = infosrc(1);
info = get(fighndl,'userdata');
dinfo = get(info.panel,'userdata');
mn = get(dinfo.activeslider,'Min');
mx = get(dinfo.activeslider,'Max');
stp = get(dinfo.activeslider,'SliderStep'); % step relative to max-min (ARG!!)
newpos = get(dinfo.activeslider,'value') + event.VerticalScrollCount * stp(1) * (mx-mn);
set(dinfo.activeslider,'value' , max( mn ,min( mx, newpos )));
clear info;
sliderupdate(src,event,infosrc);
end;

function [contrastval, mn, mx] = getcontrastval(info)
[mn, mx] = getminmax(info);
contrastval = (mx-mn)*2^(-get(info.sliders(numel(info.tileDims)+2),'value'));

function [mn, mx] = getminmax( info , imgid )
% get minimum and maximum 
if isfield(info,'panel')
    dinfo = get(info.panel,'userdata');
else
    dinfo = info;
end;
if nargin<2
    imgid = dinfo.activeColormapImage;
end;
if imgid==0
    dolog = dinfo.dolog;
    mn = info.minmax(:,1,dolog+1);
    mx = info.minmax(:,2,dolog+1);
    if numel(mn)>1
        midx = 1;
        switch dinfo.drawmethod % reduce complex data:
            case 'P'
                midx = 2;
            case 'R'
                midx = 3;
            case 'I'
                midx = 4;
        end;
        mn = mn(midx);
        mx = mx(midx);
    end;
else
    mn = info.overlays( imgid ).min;
    mx = info.overlays( imgid ).max;
end;

function updateminmax( info )
if info.doscale==2
    return;
end;
[mn, mx] = getminmax( info );
centerslider = numel(info.tileDims)+1;
centermin = get(info.sliders(centerslider),'Min');
centermax = get(info.sliders(centerslider),'Max');
if mn == centermin && mx == centermax
    return;
end
curcenter = get(info.sliders(centerslider),'value');
newcenter = ((curcenter-centermin)/(centermax-centermin))*(mx-mn)+mn;
set(info.sliders(centerslider),'Min',mn,'Max', mx,'Value',newcenter);
dinfo = get(info.panel,'userdata');
sliderpos = dinfo.sliderpos;
sliderpos(centerslider) = nan; % undefined value to force update. %newcenter;
dinfo.sliderpos = sliderpos;
set(info.panel,'userdata',dinfo);
% % update text only:
% k = numel(info.tileDims)+1; 
% set(info.texts(k), 'String' ,  sprintf(get(info.texts(k),'UserData'),newcenter)); 
% k = numel(info.tileDims)+2; 
% set(info.texts(k), 'String' ,  sprintf(get(info.texts(k),'UserData'),getcontrastval(info)));



function redrawsubplots(fighndl, updatelist)
% The drawing function. Updates all or a selection of the subplots based on
% the current focus location and contrast settings. 
info = get( fighndl, 'userdata');
dinfo = get(info.panel,'userdata');
if nargin<2
    updatelist = 1:numel(info.subfigs);
end;
doscale = info.doscale;
dolog = dinfo.dolog;
[mn, mx] = getminmax( dinfo , 0 );
sliderpos = dinfo.sliderpos;
baseindx = repmat({1},numel(info.sIMG),1);
baseindx(info.tileDims) = mat2cell(sliderpos(1:numel(info.tileDims)),ones(1,numel(info.tileDims)),1);
% if doscale~=2
%     adj = (mx-mn)*2^(-1-sliderpos(end));
%     [mn ,mx] = deal( sliderpos(end-1)-adj,sliderpos(end-1)+adj);
% end;
if ~isempty(info.overlays)
    overlays_checked = find( dinfo.activeOverlays );
end;
for m = updatelist
    if ~ishandle(info.imgs(m))
        stopimagebrowse(fighndl);
        return;
    end;
    
    curindx = baseindx;%info.baseindx;
    curdims = info.subfigdims(m,:);
    curindx(curdims) = {':'};
    curslice = info.IMG(curindx{:}); curselsize = size(curslice);
    curslice = squeeze(curslice);
    needpermute = any(diff(curdims)<0);
    if needpermute
        [dummy, sdimindx] = sort(curdims);
        curslice = ipermute(curslice,sdimindx);
    end;
    if dolog
        curslice = log(double(curslice));
    end;
    if doscale==2
        [mn,mx] = advancedminmax(curslice(:));
    end;
    
    switch dinfo.drawmethod % reduce complex data:
        case 'M'
            curslice = abs(curslice);
        case 'P'
            curslice = angle(curslice);
        case 'R'
            curslice = real(curslice);
        case 'I'
            curslice = imag(curslice);
    end;
    if ~isreal(curslice) 
        curslice = real(curslice); % avoid problems with log of negative real values
    end;

    switch doscale
    case {1,2}
        if size(curslice,3)==1
            % use colormap:
            sc = (size(info.colormap,1)-1e-11)./(mx-mn);
        else
            % true color:
            sc = (256-1e-11)./(mx-mn);
        end;
        curslice = uint8(floor((double(curslice)-mn).*sc));
    end;
    if ~isempty(info.overlays)
        if size(curslice,3)==1 
            % not true color, so apply colormap:
            curslice = applycolormap( curslice, info.colormap );
        else
            curslice = double(curslice)/256;
        end;
        for overlaynr = overlays_checked
            o = info.overlays( overlaynr );
            sz_o = size(o.image);sz_o(end+1:numel(info.sIMG))=1;
            scalardims = sz_o==1;
            curindx_o = curindx;
            if any(scalardims)
                curindx_o(scalardims) ={ 1 };
%                 repdims = curdims( scalardims( curdims(1:2) ) );
%                 for k = 1 :numel(repdims)
%                     curindx_o{ repdims(k) } = ones(info.sIMG(repdims(k)),1);
%                 end;
            end;
            perm = 1:numel(info.sIMG); perm(curdims)=[];
            overlayimg = permute( o.image( curindx_o{:} ) , [curdims perm]);
            overlaymsk = permute( o.mask( curindx_o{:} ), [curdims perm]);
            [mno, mxo] = getminmax( dinfo , overlaynr );
            sc = 1./(mxo-mno);
            overlayimg = (overlayimg-mno).*sc;
            if ~isempty( o.colormap )
                % if colormap, then overlayimg should not have color, else
                % it should
                overlayimg =  applycolormap( overlayimg, o.colormap );
            end;
            curslice = curslice + bsxfun(@times,bsxfun(@minus, overlayimg,curslice),overlaymsk);
        end;
        curslice = uint8(curslice*256);
    end;
    set(info.imgs(m),'CData',curslice);
    if ~isempty( info.ROIs ) && ~isempty( info.ROIs_subplot{m} )
        chldrn = get(info.subfigs(m), 'Children');
        delete(chldrn(1:end-1)); % last child is the image itself.
        curtiledims = find( ~ismember( info.tileDims, curdims)) ;
        selsubplots = all( bsxfun(@eq, info.ROIs_subplotidx{m}, sliderpos(curtiledims)' ) | isnan(info.ROIs_subplotidx{m}), 2) ;
        plotROIidx = info.ROIs_subplot{m}( selsubplots ) ;
        if ~isempty(plotROIidx)
            curdim1 = curdims(1);
            curdim2 = curdims(2);
            for ROIidx = plotROIidx
                hline = line( info.ROIs{ ROIidx }(:, curdim2 ), info.ROIs{ ROIidx }(:, curdim1 ) , 'parent', info.subfigs(m));
                if ~isempty(info.ROIcolors)
                    set(hline, 'color', info.ROIcolors{ ROIidx } );
                end;
            end;
        end;
    end;
end;

function info = prepareROIs( info )
ROIs_subplot = cell(1, size(info.subfigdims,1) );
ROIs_subplotidx = cell(1, size(info.subfigdims,1) );
k=1;
info.ROIs = info.ROIs(:);
while k <= numel(info.ROIs)
    if isa(info.ROIs{k},'polygon')
        pol = info.ROIs{k};
        ndimpol = size(pol.points,1);
        spatsz = info.sIMG(1:ndimpol);
        newroi = cell(10,1);
        nn = 0 ; 
        for sf = 1 : size(info.subfigdims,1)
            if all(info.subfigdims(sf,:)<=ndimpol)
                sl = pol.slice( info.subfigdims(sf,:), spatsz );
                for sli = 1: numel(sl)
                    curves = sl{sli}.get_2Dcontours;
                    ne = nn + numel(curves);
                    if ne>numel(newroi)
                        newroi(ne*2)={[]};
                    end;
                    newroi( nn+1 : ne ) =curves;
                    nn= ne;
                end;
            end;
        end;
        newroi = newroi(1:nn);
        for ne = 1:nn
            newroi{ne} = newroi{ne}([1:end 1],:); % make contour closed.
            % this is implicit in contours extracted from polygon objects,
            % but since imagebrowse allows to draw non closed contours, we
            % have to close them manually.
        end;
        if ndimpol<max(info.tileDims)
            md = max(info.tileDims);
            for ne = 1:nn
                newroi{ne}(:,end+1:md) = nan;
            end;
        end;
        info.ROIs = [info.ROIs(1:k-1);newroi;info.ROIs(k+1:end)];
        if ~isempty(info.ROIcolors)
            info.ROIcolors = info.ROIcolors([1:k-1 k*ones(1,nn) k+1:end]);
        end;
    else
        % matrix with ROI:
        extend = max(info.ROIs{k},[],1)- min(info.ROIs{k},[],1);
        extend(isnan(extend))=0; % NAN means same in all planes in that dimension.
        [dummy, srtidx] = sort( extend( info.tileDims ) ); % do not evaluate extend in 'colordim'
        subfigdims = sort( info.tileDims( srtidx(end-1:end) ) );
        subfigidx = find( all(bsxfun(@eq, info.subfigdims(:,1:2), subfigdims),2) );
        % numel( subfigidx ) ==1 , as each pair of tiledims is present
        % exactly once.
        not_inplane_dims = sort( info.tileDims( srtidx(1:end-2) ) );
        not_inplane_pos = info.ROIs{k}(: , not_inplane_dims ) ;
        rng = info.rngIMG(:, not_inplane_dims );
        not_inplane_idx = bsxfun(@times, bsxfun(@minus, not_inplane_pos, rng(1,:)) , (info.sIMG(not_inplane_dims)-1)./(rng(2,:)- rng(1,:)) )+1; 
        not_inplane_idxm = round(not_inplane_idx) ; not_inplane_idxm(isnan(not_inplane_idxm)) =0.5;
        plotidx = unique( not_inplane_idxm, 'rows');plotidx(plotidx==0.5)=nan;
        ROIs_subplotidx{subfigidx}( end+(1:size(plotidx,1)),:) = plotidx;
        ROIs_subplot{subfigidx}(end+(1:size(plotidx,1))) = k;
        k=k+1;
    end;
end;
info.ROIs_subplotidx = ROIs_subplotidx;
info.ROIs_subplot    = ROIs_subplot;

function resizefig(src,event,infosrc)%#ok : = event
% Called when figure is resized. Updates size of the uipanel with all
% orthogonal views.
info = get(infosrc,'userdata');
figsz = get(info.fighnld,'position');
if ~ishandle(info.panel)
    stopimagebrowse(infosrc);
    return;
end;
set(info.panel,'position',[110/figsz(3) 0 (figsz(3)-110)/figsz(3) 1]);

function imageclick(src,event,infosrc) %#ok: = event
% Called when clicking anywhere in the subpanels. 
% Changes the focuspoint in the dimension of the subfigure to the closest
% clicked point. First the slider positions are updated, followed by the
% slider update callback (which is not fired by the set function)
info = get(infosrc{1},'userdata');
if ~isequal(get(info.fighnld,'SelectionType'),'normal')
    return; % only process left clicks.
end;
cursorpos =  get(info.subfigs(infosrc{2}),'CurrentPoint');
cursorpos = cursorpos(1,[2 1]);
curdims = info.subfigdims(infosrc{2},1:2);
cursorpos = round((cursorpos-info.rngIMG(1,curdims))./ ([-1 1]*info.rngIMG(:,curdims)).* (info.sIMG(curdims)-1)+1);
if all(cursorpos>0 & cursorpos<=info.sIMG(curdims))
    set(info.sliders(info.tileDims==curdims(1)),'Value',cursorpos(1,1));
    set(info.sliders(info.tileDims==curdims(2)),'Value',cursorpos(1,2));
end;
sliderupdate(info.fighnld,[],info.fighnld);

function stopimagebrowse(infosrc)
% Try to 'gracefully' shutdown the imagebrowse figure so that we wont throw any
% errors after this call.
%
% Do not clear the figure entirely, since typically this is called after
% (accidentely) part of the figure is painted over (by a call to plot) and
% we don't want to destroy that content. 
disp('Imagebrowse figure corrupt; stopping imagebrowse. This is usually due to ''external'' plotting or clf.');
info = get(infosrc,'userdata');
set(info.fighnld,'ResizeFcn',[]);
set(info.fighnld,'userdata',[]);
set(info.fighnld,'WindowButtonMotionFcn',[]);
for k=1:numel(info.sliders)
    if ishandle(info.sliders(k))
        set(info.sliders(k),'callback',[]);
    end;
end;
for m=1:numel(info.imgs)
    if ishandle(info.imgs(m))
        set(info.imgs(m),'ButtonDownFcn',[]);
    end;
end;

function selcbk(source, eventdata, infosrc)%#ok
% Called when changing the selection Magnitude/real/imag/phase:
if ~isequal(eventdata.OldValue,eventdata.NewValue)
    info = get(infosrc,'userdata');
    dinfo = get(info.panel,'userdata');
    dinfo.drawmethod = get(eventdata.NewValue,'Tag');
    set(info.panel,'userdata',dinfo);
    updateminmax(info);
    redrawsubplots(info.fighnld);
end;

function togglelog(source,eventdata, infosrc)%#ok
% Called when the 'log' toggle button changed status (checked<->unchecked)
% if ~isequal(eventdata.OldValue,eventdata.NewValue)
    info = get(infosrc,'userdata');
    dinfo = get(info.panel,'userdata');
    dinfo.dolog = (get(source,'Value') == get(source,'Max'));
    set(info.panel,'userdata',dinfo);
    updateminmax(info);
    redrawsubplots(info.fighnld);
% end;

function setzoom(source,eventdata, infosrc)%#ok
info = get(infosrc,'userdata');
for m=1:numel(info.subfigs)
    if ishandle(info.subfigs(m))
        ax = info.subfigs(m);
        set(ax,'units','pixels');
        loc = get(ax,'Position');
        set(ax,'units','normalized');
        xl = get(ax,'xlim');
        yl = get(ax,'ylim');
        set(ax,'xlim',mean(xl)+[1 loc(3)]-(1+loc(3))/2,'ylim',mean(yl)+[1 loc(4)]-(1+loc(4))/2);
    end;
end;    

function [bool] = isvalidImageBrowsefigure( curfig )
% Try to determine if curfig is created/filled by imagebrowse and it is
% still valid.
bool = ~isempty(curfig);
if ~bool, return;end;
%  bool = all(isfield(get(curfig,'UserData'), {'IMG','sIMG','rngIMG'}));  % isfield does not support cell list of strings on old MATLAB, so use replacement code below:
reqfields = {'IMG','sIMG','rngIMG'}; 
ud = get(curfig,'UserData');
for k=1:numel(reqfields)
    bool = bool && isfield( ud, reqfields{k} ); % assume true imagebrowse figure when these fields are present
end;
if ~bool, return;end;
ud = get(curfig,'UserData');
bool  = ishandle(ud.fighnld) && isequal(ud.fighnld,curfig) && ishandle(ud.panel) && all(ishandle([ud.subfigs;ud.imgs]));
%if ~bool, return;end;

function afterzoomupdate( obj, event_obj) 
% After a zoom action, update the range in all subfigures that contain one
% of the dimensions of the zoomed subfigure. 
curfig = obj;
ud = get(curfig,'UserData');
updatedaxisindx = find( ud.subfigs == event_obj.Axes );
updateddims = ud.subfigdims( updatedaxisindx ,1:2);
for hr = 1 :2
    if hr==1
        newlim = get( event_obj.Axes, 'Ylim');
    else
        newlim = get( event_obj.Axes, 'Xlim');
    end;
    for updax = find(ud.subfigdims(:,1)==updateddims(hr))'
        if updax~=updatedaxisindx
            set( ud.subfigs( updax) ,'Ylim',newlim);
        end;
    end;
    for updax = find(ud.subfigdims(:,2)==updateddims(hr))'
        if updax~=updatedaxisindx
            set( ud.subfigs( updax) ,'Xlim',newlim);
        end;
    end;
    
end;

function [img ] = applycolormap( img, cmap )
% [img_out] = applycolormap( img_in, cmap);
%
% img_in = n x m, double between 0 and 1 or uint8 (0..255)
% cmap = n_colors x 3
% img_out = n x m x 3

if isa( img,'uint8')
    indx = min(double(img)+1, size(cmap,1)) ;
else
    indx = max(1,min( size(cmap,1), round(img * size(cmap,1))+1));
end;
img = reshape( cmap( indx , :), [size(img) size(cmap,2)]);

function update_colorbar(figno, activeboxnr)
info = get(figno,'userdata');
dinfo = get(info.panel,'userdata');
activeOverlays = dinfo.activeOverlays;
colorbarorder = [1:numel(info.overlays) 0];
if nargin>=2 
    if activeboxnr==0
        colorbarorder = colorbarorder([end 1:end-1]);
    elseif activeOverlays( activeboxnr )
        colorbarorder(1:activeboxnr) = colorbarorder([activeboxnr 1:activeboxnr-1]);
    end;
end;
for k = colorbarorder
    if k==0 || activeOverlays(k)
        if k==0
%             cmap = info.colormap;
            cname = '';
            set(info.colorbarimage,'CData',info.colorbarimagemats{1});
        else
            cmap = info.overlays(k).colormap;
            cname = info.overlays(k).name;
            set(info.colorbarimage,'CData',applycolormap( info.colorbarimagemats{2}, cmap ) );
        end;
%         colormap( cmap );
        dinfo.activeColormapImage = k;
        set( info.panel,'userdata', dinfo );
        if isempty(cname)
            set( info.colorbarname, 'String', cname,'visible','off');    
        else
            set( info.colorbarname, 'String', cname,'visible','on');
        end;
        break;
    end;
end;

% 

function checkboxchange( obj, dummy, checkboxnumber, checkboxgroup)
figno = get(checkboxgroup,'parent');
info = get(figno,'userdata');
dinfo = get(info.panel,'userdata');
dinfo.activeOverlays(checkboxnumber) = get(obj,'Value');
set(info.panel,'userdata',dinfo);
update_colorbar( figno, checkboxnumber );
updateminmax( info );
sliderupdate(info.fighnld,[],[info.fighnld 0]);
% redrawsubplots( figno );