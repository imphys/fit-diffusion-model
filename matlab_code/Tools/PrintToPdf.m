function PrintToPdf(IMGName,size,fmts)
% PrintToPdf(IMGName [, size [,fmts]])
% Saves the current figure to a pdf named IMGName. Resolution of the saved
% figure is fixed to 800 dpi. Size [width height] is the size of the figure
% in cm (even when that is not the default). When not specified the size
% defaults to 10*8 cm, and the papersize is set accordingly.
% fmts is an optional string or cell array with the formats in which it
% needs to be saved, both print 
 
% Copyright 2/2/2007: Dirk Poot, University of Antwerp

if nargin<2 || isempty(size)
    size = [10 8];
end;
if nargin<3
    fmts = {'pdf'};
end;
if ischar(fmts)
    fmts = {fmts};
end
translationtbl = {'pdf' , '-dpdf' 
                  'emf' , '-dmeta'
                  'eps' , '-depsc2'};% left column entries should be unique.
if any(IMGName=='.')
    IMGName = strrep(IMGName,'.','_');
    disp(['image name contained ''.''s , replaced to get the new name:' IMGName]);
end;
curfig = gcf;
set(curfig,'PaperUnits','centimeters');
set(curfig,'papersize',size);
set(curfig,'PaperPosition',[0 0 size]);
drawnow;pause(0.1); % allow updates and possibly interrupt.
modifications = 0;
for k=1:numel(fmts)
    try 
        fmt = fmts{k};
        tr = find(strcmpi({translationtbl{:,1}},fmt));
        if ~isempty(tr)
            fmt = translationtbl{tr,2};
        end;
        if strcmpi(fmt,'-dmeta') && modifications==0
            % meta file (.emf) behaves strangely when no resize function is present (ARG!!) bugreport: 231161.
            rf = get(curfig, 'resizeFcn');
            resizeWarn = '';
            % put in dummy resize function and disable CustomResizeFcnInPrint warning
            if isempty(rf) 
              set(curfig, 'resizefcn', 'disp('''');');
              resizeWarn = warning('off', 'MATLAB:Print:CustomResizeFcnInPrint');
            end
        end;
        print('-f',fmt,'-r800', IMGName);
%         print('-f','-dpdf','-r800', IMGName);
%     print('-f','-dpsc2','-r800', IMGName);
%     print('-f','-dmeta','-r800', IMGName);
    catch
        a= lasterror;
        if strcmp(a.identifier,'MATLAB:Print:CannotCreateOutputFile')
            answ = input('File cannot be created, probably locked. Press <Enter> to retry (''N'' to stop). ','s');
            if ~(isequal(answ,'N') | isequal(answ,'n'))
                if ishandle(curfig)
                    print('-f',fmt,'-r800', IMGName);
                else
                    warning('PrintToPDF:figgone','figure destroyed, so cannot print it anymore');
                end;
            end;
        else
            if modifications
                % restore previous state   
                if isempty(rf)
                  set(curfig, 'resizeFcn', '');
                  warning(resizeWarn);
                end
            end;            
            error(a);
        end;
    end;
end;
if modifications
    % restore previous state   
    if isempty(rf)
      set(curfig, 'resizeFcn', '');
      warning(resizeWarn);
    end
end;