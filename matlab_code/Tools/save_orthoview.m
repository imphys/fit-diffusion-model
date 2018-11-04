function save_orthoview( img , focuspoint , upsampfact , filenamebase )
% save_orthoview( img , focuspoint , upsampfact , filenamebase );
% Saves 3 orthogonal view through focuspoint of the 3D image img.
% as png's, 'neirest neighbor' upsampled in each dimension by 
% an integer factor upsampfact to the files [filenamebase num2str(i) '.png'], 
% for the views i=1:3
%
% If img is integer, it is passed as is, otherwise it should be scaled
% between 0 (Black) and 1 (white).
% A 4D image with the fourth dimension of size 3 is output as color image.
%
% If image a scalar, and ishandle(img), it is assumed to be an
% imagebrowse figure handle. When focuspoint is non empty it is set and the
% orthogonal views are stored, including the colorbar of the imagebrowse figure.
%
% Created by Dirk Poot, ErasmusMC, 28-2-2012

if numel(img)==1 && ishandle(img)
    h = img;
    % imagebrowse figure handle
    info = get( h, 'UserData');
    if ~isempty(focuspoint)
        for k=1:numel(info.tileDims)
            set( info.sliders( k ) , 'Value', focuspoint( info.tileDims(k) ) );
        end;
        if numel(info.tileDims)>1
            cb = get(info.sliders(1),'Callback');
            cb{1}([], [], cb{2});
        end;
    end;
    drawnow;pause(.1);
    if numel(upsampfact)==1
        upsampfact = upsampfact(ones(1,max(info.subfigdims(:))));
    end
    sel = cell(1,size(info.subfigdims,2));sel{end}=':';
    for k=1:size(info.subfigdims,1)
        tmp = get( info.imgs(k),'CData' );
        for d = 1:2
            td = info.subfigdims( k, d );
            sel{d} = floor( (0:upsampfact(td)*size(tmp,d)-1)/upsampfact(td))+1;
        end;
        fn = [filenamebase sprintf('%d',info.subfigdims(k,:)) '.png'];
        if size(tmp,3)==1
            imwrite(tmp(sel{:}), get(h,'colormap') , fn, 'png');
        else
            imwrite(tmp(sel{:}), fn, 'png');
        end;
    end;
    lf = figure(); % create new figure;
    oldpos = get(info.colorbar,'position');
    oldunits = get(info.colorbar,'units');
    bottom = 0.02;
    top = .98;
    set(info.colorbar,'parent',lf,'Units','normalized','Position',[.01 bottom .3 top-bottom]);
    set(lf,'colormap',get(h,'colormap'));
    
    ytick = get(info.colorbar,'Ytick');ytick = ytick(:);
    yticklbl = str2num(get(info.colorbar,'YTickLabel'));
    nonzero = isfinite(ytick)& isfinite(yticklbl) & yticklbl~=0 & ytick~=0;
    tickscale = log( median( ytick(nonzero)./yticklbl(nonzero)));
    ylim_val = get(info.colorbar,'Ylim');
    if abs(tickscale)>1
        top = .89; % x 10^(b)  included in axis.
    else
        % reduce hight when largest tick is close to top
        top  = min(1, .97- (max(ytick)-max(ylim_val))/diff(ylim_val)); 
    end;
    % increase bottom spacing when smallest tick is close to bottom:
    bottom = max(.02, .04-(min(ytick)-min(ylim_val))/diff(ylim_val));
    set(info.colorbar,'Position',[.01 bottom .3 top-bottom]);

    fn = [filenamebase 'legend'];
    drawnow;pause(.1); % make sure all drawing is complete.
    PrintToPdf( fn, [1.5 5]);
    drawnow;pause(.1); % make sure all drawing is complete.
    set(info.colorbar,'parent',h,'Units',oldunits,'Position',oldpos);
    close(lf);
    return;
end;

sel = repmat({':'},1,4);
if numel(upsampfact)==1
    upsampfact = [1 1 1 1]*upsampfact;
end;
for k=1:3
   sel{k} = floor( (0:upsampfact(k)*size(img,k)-1)/upsampfact(k))+1;
end;
for k=1:3
    fn = [filenamebase num2str(k) '.png'];
    selk = sel;
    selk{k} = focuspoint(k);
    selimg = permute( img(selk{:}) ,[1:k-1 k+1:4 k]);
    imwrite(selimg, fn, 'png');
end;


