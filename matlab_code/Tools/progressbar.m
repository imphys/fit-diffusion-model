function [pos] = progressbar(varargin)
% Displays a multi leveled progressbar. This makes it easy to nest
% computations and still have a easy interpretable progress indication.
%
% The normal usage is the subsequent calling of these 3 functions:
%
% progressbar('start', [firstVal lastVal], description , ... ); 
%        alternative: progressbar('descend', ... );  
%     Creates a new progressbar (when none exists), or adds a new level to
%     the current progressbar. 
%   [lastVal] or [firstVal lastVal]:
%     This is the range that the progressbar should span, 
%     firstVal (default = 0) implies 0% progress, 
%     lastVal (default = 1) implies 100% progress. 
%   description:
%     Set description to display the current action in the progressbar title.
%   remaining inputs are treated as option value pairs, see further down.
%
% progressbar(newPos, caption)
%   newPos: a scalar numerical argument.
%     Sets the current progressbar level to the desired position, the
%     percentage is computed with respect to the firstVal-lastVal value pair
%     given in the start of this level.
%   caption: optional string argument
%     Set caption to the axes title, display specific progess information.
%     (eg: How many loops have been done)
%
% progressbar('ready' , ... );
%        alternative: progressbar( 'ascend' , ... ); 
%   Completes current level of computation, when this was the last level
%   the progressbar is closed. Use this function at the end of a loop to
%   finish the progressbar.
%
%
% Other/ more advanced options, described below:
% 'parfor_start', 'position', 'close', 'EstTimeLeft', 'minTimeInterval',
% 'startStatistics', 'getStatistics', 'clearTimeEstimation', 'DebugButton',
% 'hasGUI'
% These options can be combined; also with 'start' or 'ready' arguments
%
% progressbar('parfor_set_progressdir', progressfilepath )
%     Set path before starting parfor progress indication ('parfor_start')
%   progressfilepath: Path that is accessibly by all workers and the
%     client. Used to store the progress of each of the workers, using at
%     least 2 second intervals. The file system is used since the parallel
%     computing toolbox allows no custom comunication between workers and
%     client.
% progressbar('parfor_start', numCalls [, description [,..]])
%     Works (almost) as 'start', but now for parfor loops. 
%     You should set the progressbar file path ('parfor_set_progressdir') first.
%   numCalls: The integer number of times progressbar is called within
%     the loop (i.e. by the workers), usually the upper limit of the parfor
%     loop variable.  
%   Just use 'ready' (as usual) when the parfor loop is finished.
%   Note that parfor has quite some initialisation and finalisation
%   overhead (copy of data), therefore estimating the remaining computation time
%   might not be accurate in all cases. Also set the 'minTimeInterval' to
%   (approximately) at least the time required for 1 worker to compute 1
%   iteration of the parfor loop. 
%   When you are not in paralel mode, a normal progressbar is started,
%   since the parfor loop is then evaluated as normal for loop.
%
% [pos] = progressbar('position');
%   Last position the current progressbar was set too, as fraction completed 
%   NaN is returned when no progressbar is available.
%
% progressbar('close');
%   Resets and closes the progressbar, use only to reset after a manual
%   break of the computations.
%
% progressbar('EstTimeLeft', 'on'/'off');
%   Display an estimate of the time remaining. This is only a good estimate
%   when (on the average) the progress is linear in time. An uncertainty
%   estimate is included as well.
%
% progressbar('minTimeInterval', interval);
%   Sets the minimum amount of clock time (in seconds) that has to have
%   passed since the last update of the progressbar, sets it for the
%   current level this value is the default for all levels below. 
%
%   Default value = 0 (all updates are displayed);
%
%   This option is provided since updating the progressbar may take a
%   significant amount of time when it's done very often (in a fast loop)
%   and more than (say) 10 updates per second are (usually) not useful for
%   the user.  
%
% progressbar('UseCallCount', tf);
%   Boolean option which, if true, specifies that I should use the number of times
%   progressbar is called with a value (progressbar(newPos, caption)), instead of
%   using newPos to update the position. (So, newPos is completely ignored after it
%   is verified that it is a number)
%
% progressbar('startStatistics', count);
%  Attempts to get count samples with statistics data from the current
%  progressbar. Call this at start time to get the statistics sampled over
%  the entire progress. Clears previous statistic gathering of the current
%  progressbar.
%
% progressbar('getStatistics');
%  Reads the gathered statistics of the current progressbar (so call it
%  before calling progressbar('close');). The results is a structure array
%  with interesting statistic fields about time and memory usage.
%
% progressbar('clearTimeEstimation');
%   Clears the time estimation history. Use this when the previous points
%   do not give an good indication of the computation time of the remaining
%   part. (for example when you restart the computation at a more accurate
%   level)
%
% progressbar('DebugButton','show'); / progressbar('DebugButton','hide');
%   Shows/hides a button with a caption 'Debug'. When this button is
%   pressed execution is stopped (dbstop). This allows you to dynamically
%   stop to debug your code, which is not possible otherwise. Note that
%   click events are only processed during pause or drawnow commands.
%   (progressbar uses drawnow to update the progressbar, but you might want
%   to add additional drawnow's in your code to more frequently allow a
%   break, as well as allowing screen updates)
%
% progressbar('debugstop');
%   Issues a breakpoint after processing all progressbar commands.
%
% progressbar('hasGUI','no');
%   Disables all gui work, effectively disables all progress indication.
%   re-enable with: progressbar('hasgui','yes');
%
% Copyright Dirk Poot (Dirk.Poot@ua.ac.be)
%           University of Antwerp
% 

% Revision history, original version 20-6-2006
% 22- 6-2006 : lots of small improvements/debugging.
%  4- 8-2006 : added 'getStatistics' 
% 11- 9-2006 : added 'DebugButton', improved timeleft estimation
% 20- 9-2006 : enabled using in MATLAB 6.5 (might be broken again?)
% 30-10-2006 : improved handling of deleted object
% 13- 8-2007 : don't divide by zero on interval [0 0].
%  6- 9-2007 : improved help text.
% 13- 9-2007 : removed MLint warnings, also use short circuit operators (so progressbar doesn't work on MATLAB < 6.5 anymore)
%  3- 1-2008 : set default 'debugbutton' state to 'show'
% 17- 7-2008 : Added parfor (distributed computing) support. 
% 22- 7-2008 : added 'hasgui' option to disable graphical output.
% 25-11-2008 : cleaned up code and added support for multiple segments in
%              the progressbar
% 23- 9-2008 : added that workers can start sub progressbar's.

persistent ProgBData
persistent GeneralData
if nargin<1
    error('progressbar needs arguments');
end;

if ischar(varargin{1})
  remInp = varargin;
  if isempty(GeneralData)
      if exist('timer','file')==2
          try
              GeneralData.timer = timer('TimerFcn','progressbar(''Redraw'');', 'StartDelay', .5 , 'TasksToExecute',1);
          catch %#ok : catch any error during timer creation
              % Most likely cause: no Java (activated)
              % When error during timer creation: Don't use a timer.
              GeneralData.timer = [];
          end;
      else
          GeneralData.timer = [];
      end;
      GeneralData.CurrIndex = 0;
      GeneralData.NeedDraw = 0;
      GeneralData.ProgB = [];
      GeneralData.doDebugStop = 0;
      GeneralData.DebugButton = [];
      GeneralData.ParallelState = 0;
      GeneralData.ParallelBarIndx = nan;
      GeneralData.parfor_progressdir = '';
%       GeneralData.CallCountAtLevel = 0;
      GeneralData.hasGUI = usejava('jvm'); % detect if we have a gui.
      GeneralData.NextUpdate = now;
  end;
  while length(remInp)>=1
    index = GeneralData.CurrIndex;
    usedInputs = 2;
    switch lower(remInp{1})
        case {'descend', 'start'}
            % Note that the default figure settings of the progressbar are set
            % when the figure is created (which is done when needed in the
            % 'redraw' section.
            if isempty(ProgBData)
%                 ProgBData=struct('StartStop',[],'minTimeInterval',[],'AxHndl',[],'PatchHndl',[],'title',[]);
%                 ProgBData=struct('StartStop',{},'minTimeInterval',{},'AxHndl',{},'PatchHndl',{},'title',{},'EstTimeHndl',{},'Times',{},'LastUpd',{},'UpdHist',{},'nUpdHist',{},'Statistics',{},'StatisticsCount',{});
                ProgBData=struct('StartStop',{},'minTimeInterval',{},'AxHndl',{},'PatchHndl',{},'title',{},'EstTimeHndl',{},'UpdHist',{},'nUpdHist',{},'Statistics',{},'StatisticsCount',{},'UseCallCnt',{},'CallCount',{});
                % StartStop : [start stop] indices, start = 0%; stop = 100%
                % minTimeInterval: minimum time interval, in units of 'now'.
            end;
            index = index+1;
            GeneralData.CurrIndex = index;
            if index <= 1
                ProgBData(1).minTimeInterval = 0;
                ProgBData(1).title = 'Progressbar';
            else
                ProgBData(index).minTimeInterval = ProgBData(index-1).minTimeInterval;
                % length ProgBData is increased by line above!
                ProgBData(index).title = ProgBData(index-1).title; 
            end;
            if length(remInp)>=3
                if ~isempty(remInp{3})
                    ProgBData(index).title = remInp{3}; 
                end;
                usedInputs = 1;

%                 if rem(nargin-3,2)
%                     error('option value pairs do not match');
%                 end;
                remInp{3} = remInp{2};
            else
                usedInputs = 0;
                if length(remInp)==1
                    remInp{2} = 1;
                end;
            end;
            remInp{usedInputs+1} = 'adjustlimits';
            ProgBData(index).StartStop = [0 1];
%             ProgBData(index).LastUpd = [value(1) now];
            ProgBData(index).UpdHist = [0; now];
            ProgBData(index).nUpdHist = 1;
            ProgBData(index).UseCallCnt = false;
            ProgBData(index).Statistics = [];
            ProgBData(index).CallCount  = 0;
            
            GeneralData.NextUpdate  = ProgBData(index).UpdHist(2,1) + ProgBData(index).minTimeInterval;
            
            if GeneralData.hasGUI %GeneralData.ParallelState(1)~=2 % when not worker in parallel mode:
                GeneralData.NeedDraw = GeneralData.NeedDraw + 1;
    %             ProgBData(index).EstTimeHndl = [];
    %             progressbar('reposition');
    %             remInp{end+1} = 'reposition';
                if ishandle(ProgBData(index).EstTimeHndl)
                    set(ProgBData(index).EstTimeHndl,'String','');
                end;
                if ~isempty(GeneralData.timer)  && ishandle(GeneralData.timer(1))
                    if ~strcmp(get(GeneralData.timer(1),'Running'),'on')
                        start(GeneralData.timer(1));
                    end;
                else
                    remInp{end+1} = 'Redraw'; %#ok: expand inside loop: we dont know how large it will get.
                end;
            end;
            pos = index; % return current index;
        case 'adjustlimits'
            %progressbar('AdjustLimits',progbbase + 2*nnz(needtodo));
            if length(remInp)<2 || length(remInp{2})<1
                value = 1;
                if length(remInp)<2
                    usedInputs = 1;
                end;
            else
                value = remInp{2};
            end;
            if size(value,1)==2 && size(value,2)==1
                warning('PROGRESSBAR:Ambiguous_input','Ambiguous start stop. Old version would interpret this differently');
            end;
            if size(value,2)<2
                value = [zeros(size(value,1),1) value]; %#ok : expand when input provided has length 1
            end;
            adj = abs(value(:,1) - value(:,2)) < eps* max(abs(value),[],2);
            if any(adj)
                value(adj,2) = value(adj,1) + 200*(realmin+eps * max(abs(value(adj,:)),[],2)); %#ok : lint remark not true;
            end;
            ProgBData(index).StartStop = value;
%             if isnan(ProgBData(index).UpdHist(1))
%                 ProgBData(index).UpdHist(1) = value(1);
%             end;
        case 'parfor_set_progressdir'
            storeName = remInp{2};
            if ~isempty(storeName)
                if storeName(end)~=filesep
                    storeName = [storeName filesep];%#ok : just add filesep, not grown as implied by mlint.
                end;
                if exist(storeName,'dir')~=7
                    error('PROGRESSBAR:dirnotfound','Client cannot find directory specified for parallel progressbar');
                end;
                storeName = [storeName 'progressbar'];%#ok : just add 'progressbar', not grown as implied by mlint.
            end;
            GeneralData.parfor_progressdir = storeName;
        case {'parfor_start'}
            if ~inpmode
                warning('PROGRESSBAR:NotParallel','Not inside parallel computing environment, starting normal progressbar');
                if numel(remInp)<3
                    remInp(end+1:3) = {[]};
                end;
                remInp(end+1:end+2) = {'UseCallCount', true};
            else
                GeneralData.ParallelState(1) = 1; % indicate client
                for f_idx=2:numel(GeneralData.ParallelState)
                    if GeneralData.ParallelState(f_idx)>0
                        try 
                            fclose(GeneralData.ParallelState(f_idx));
                        catch %#ok not interested in a close error.
                        end;
                    end;
                end;
                GeneralData.ParallelState = GeneralData.ParallelState(1);
                storeName = GeneralData.parfor_progressdir;
                if isempty(storeName)
                    warning('PROGRESSBAR:NotParallelDir','Directory in which the workers store their progress not initialized. \nThe progressbar is started (to avoid errors), but it will not update. \nCall progressbar(''parfor_set_progressdir'', dir) to set the directory.');
                end;
                pctRunOnAll(['progressbar(''parfor_start_worker'',''' storeName ''');']);
                if prod(size(GeneralData.timer))==1  %#ok<PSIZE>
                    GeneralData.timer(2) = timer('TimerFcn','progressbar(''parfor_updClient'');','StartDelay', 2 , 'TasksToExecute',inf,'BusyMode','drop','ExecutionMode','fixedSpacing','Period',2);
                end;
                GeneralData.progrStoreName = storeName;
                if prod(size(GeneralData.timer))<2 %#ok: bug in numel of timer objects.
                    warning('Progressbar:NoTimer','Timer not found, client displays no updates');
                else
                    if ~strcmp(get(GeneralData.timer(2),'Running'),'on')
                        start(GeneralData.timer(2));
                    end;
                end;
                GeneralData.ParallelBarIndx = index+1;
            end;
            remInp{1} = 'start';
            usedInputs = 0;
        case 'parfor_start_worker'
            if ~exist('numlabs','builtin') || numlabs==1 
                if GeneralData.ParallelState(1)~=1
                    warning('PROGRESSBAR:NotParallel','Cannot start worker progressbar when not in parallel mode');
                end;
            elseif GeneralData.ParallelState(1)~=1
                for f_idx=2:numel(GeneralData.ParallelState)
                    if GeneralData.ParallelState(f_idx)>0
                        try 
                            fclose(GeneralData.ParallelState(f_idx));
                        catch %#ok not interested in a close error.
                        end;
                    end;
                end;
                GeneralData.ParallelState = [2 0]; % indicate: worker, no open file, 
                GeneralData.hasGUI = 0; % worker has no gui.
                if GeneralData.CurrIndex~=0
%                     GeneralData % debug
%                     ProgBData   % debug
                    warning('PROGRESSBAR:WrongLevel','Worker progressbar started not at first level; probably a previous bar was not finished properly');
                end;
                GeneralData.ParallelBarIndx = index+1;
                if isempty(remInp{2})
                    storeName = mfilename('fullpath') ;
                else
                    storeName = remInp{2};
                end;
                storeName = [storeName '_worker_' num2str(labindex) '.prb']; %#ok: is not growing
                GeneralData.ParallelState(2) = fopen( storeName,'w');
                if GeneralData.ParallelState(2)<0
                    GeneralData.ParallelState(:) = 0;
                    error('Cannot open worker progress output file');
                end;
                fwrite(GeneralData.ParallelState(2),0,'uint32');                
                if index>0
                    % close all progressbar indicators.
%                     remInp(end+1:end+index) = repmat({'ready'},index,1);
                    GeneralData.CurrIndex = 0;
                end;
                remInp(end+1:end+7) = {'start',[],[],'mintimeinterval',2,'UseCallCount', true};
            end;
        case {'parfor_updclient'}
%             disp('update client'); % debug
            if length(GeneralData.ParallelState)==1
                f_idx = 1; dolp = true;
                storeName = GeneralData.progrStoreName;
                while f_idx<4 || dolp 
                    WrkrstoreName = [storeName '_worker_' num2str(f_idx) '.prb'];
                    if exist(WrkrstoreName,'file')
                        GeneralData.ParallelState(f_idx+1) = fopen(WrkrstoreName,'r');
                    else
                        dolp =false;
                    end;
                    f_idx = f_idx +1 ;
                end;
            end;
            curpos = zeros(numel(GeneralData.ParallelState)-1,1);
            for f_idx=2:numel(GeneralData.ParallelState)
                if GeneralData.ParallelState(f_idx)>0
                    fseek(GeneralData.ParallelState(f_idx),0,'bof');
                    try
                        tmp = fread(GeneralData.ParallelState(f_idx), inf, 'single');
                        tmpel = find(isnan(tmp),1,'first')-1;
                        curpos(f_idx-1,1:tmpel)  = tmp(1:tmpel);
                    catch %#ok<CTCH>
                        % dont want to warn user repeatively that a file could not be
                        % read.
                    end;
                end;
            end;
            barcntdiff = size(curpos,2) - (index-GeneralData.ParallelBarIndx+1);
            if barcntdiff>0
                args = repmat({'start',ones(size(curpos,1),1),[]}, 1, barcntdiff);
                progressbar(args{:});
            elseif barcntdiff<0
                GeneralData.CurrIndex = GeneralData.ParallelBarIndx+size(curpos,2)-1;
                remInp{end+1} = 'redraw'; %#ok<AGROW>
            end;
            for barindx=1:size(curpos,2)
                lpidx = GeneralData.ParallelBarIndx+barindx-1;
                GeneralData.NextUpdate  = ProgBData(lpidx).UpdHist(2,1) + ProgBData(lpidx).minTimeInterval;
                progressbar(curpos(:,barindx),[num2str(sum(GeneralData.ParallelState>0)-1) ' workers.'], lpidx );
            end;
            usedInputs = 1;
        case {'ascend','ready'}
            usedInputs = 1;
            if GeneralData.ParallelState(1)==1 || (GeneralData.ParallelState(1)==2 && GeneralData.CurrIndex==GeneralData.ParallelBarIndx)
                isclient = GeneralData.ParallelState(1)==1;
%                 if isclient % debug
%                     disp('progressbar ready on client')
%                 else
%                     disp('progressbar ready on worker')
%                 end;
                if isclient && prod(size(GeneralData.timer))>=2 %#ok:bug in numel of timer objects
                    stop(GeneralData.timer(2));
                end;
                for f_idx=2:numel(GeneralData.ParallelState)
                    if GeneralData.ParallelState(f_idx)>0
                        fclose(GeneralData.ParallelState(f_idx));
                    end;
                end;
                nfiles = numel(GeneralData.ParallelState)-1;
                GeneralData.ParallelState = 0; % Wanted: calls 'ready' in client with not parallel set!
                if isclient
                    % Close bars created by workers: (redraw by subsequent pctRunOnAll('ready'))
                    GeneralData.CurrIndex = GeneralData.ParallelBarIndx;
                    % Close client created bar on workers (and client as well)
                    if inpmode
                        % warn the workers to quit, but only when matlabpool still
                        % running.
                        pctRunOnAll('progressbar(''ready'');');
                    else
                        usedInputs = 0; % still need to do the actual 'ready' 
                        % note closing files will probably not work since I can't
                        % notify the workers to close it.
                    end;
                    % Cleanup files:
                    storeName = GeneralData.progrStoreName;
                    for f_idx=1:nfiles
                        WrkrstoreName = [storeName '_worker_' num2str(f_idx) '.prb'];
                        if exist(WrkrstoreName,'file')
                            delete(WrkrstoreName);
                        end;
                    end;
                else
                    GeneralData.CurrIndex = index-1; % clear worker progressbar.
                end;
            else
%                 disp('normal progressbar ready'); %debug
                if index<=0
                    warning('PROGRESSBAR:cannotascend','Cannot ascend to higher level when no higher level is present.');
                    return;
                end;
                if ishandle(ProgBData(index).EstTimeHndl)
                    delete(ProgBData(index).EstTimeHndl);
                end;
                index = index-1;
                GeneralData.CurrIndex = index;
                if index>0
                    GeneralData.NextUpdate  = ProgBData(index).UpdHist(2,1) + ProgBData(index).minTimeInterval;
                end;
                if ~isempty(GeneralData.timer) 
                    if ~strcmp(get(GeneralData.timer(1),'Running'),'on')
                        start(GeneralData.timer(1));
                    end;
                else
                    remInp{end+1} = 'Redraw';  %#ok: expand inside loop: initially we dont know how large it will get.
                end;
                GeneralData.NeedDraw = max(0,GeneralData.NeedDraw -1);
            end;
        case 'redraw'
            usedInputs = 1;
          if GeneralData.hasGUI
            if isempty(GeneralData.ProgB) || ~ishandle(GeneralData.ProgB)
                % create figure when there is none.
%                 oldFig = [];
%                 if length(findobj('Type','figure'))>0
%                     oldFig = gcf;
%                 end;
                GeneralData.ProgB = 335965207; % random but constant figure number, high number to minimizes figure number clashes;
                while ishandle(GeneralData.ProgB)
                    GeneralData.ProgB = floor( rand * 2030832096 +1);
                end;
                GeneralData.ProgB = figure(GeneralData.ProgB);
                set(GeneralData.ProgB,'MenuBar', 'none', ...
                                'NumberTitle','off' ...
                                ,'HandleVisibility','off'...
                                ,'Visible','off'...
                  );
                % set default values, allow to be overridden by inputs that
                % follow.
                if ~isdeployed
                    remInp(2:end+2) = [{'DebugButton','show'}, remInp(2:end)];
                end;
%                 % set focus back.
%                 if ~isempty(oldFig) 
%                     figure(oldFig);
%                 end;
            end;
            if GeneralData.NeedDraw>0
            
                for k= GeneralData.NeedDraw-1:-1:0
                    isnotvalid = isempty( ProgBData(index-k).AxHndl );
                    if ~isnotvalid 
                        isnotvalid = ~ishandle(ProgBData(index-k).AxHndl);
                    end;
                    if isnotvalid
%                       if ~activatedFig  
%                           set(GeneralData.ProgB,'HandleVisibility','on');
%                           figure(GeneralData.ProgB);
%                           activatedFig = 1;
%                       end;
                      ProgBData(index-k).AxHndl = axes('XLim',[0 1],...
                        'Box','on', ...
                        'YLim',[0 1],...
                        'XTickMode','manual',...
                        'YTickMode','manual',...
                        'XTick',[],...
                        'YTick',[],...
                        'XTickLabelMode','manual',...
                        'XTickLabel',[],...
                        'YTickLabelMode','manual',...
                        'YTickLabel',[]...
                        ,'parent',GeneralData.ProgB ...
                        );
                    end;
                    xpatch = [0 0 0 0]';
                    ypatch = [0 0 1 1]';
                    if ishandle(ProgBData(index-k).PatchHndl) % does handle empty correctely.
                        set( ProgBData(index-k).PatchHndl,'XData',xpatch);
                    else
%                         set(GeneralData.ProgB,'HandleVisibility','on');
%                         axes(ProgBData(index-k).AxHndl);
%                         activatedFig = 1;
                        ProgBData(index-k).PatchHndl(2) = patch(xpatch,ypatch,[1 .5 0],'parent',ProgBData(index-k).AxHndl);
                        ProgBData(index-k).PatchHndl(1) = patch(xpatch,ypatch,'y','parent',ProgBData(index-k).AxHndl);
                    end;
                    if ProgBData(index-k).EstTimeHndl == -1
                        remInp(end+(1:2)) = {'MakeTimeLeftText', index-k}; 
                    end;
                end;
                GeneralData.NeedDraw = 0;
            end;
            if length(ProgBData)>index
                for k=length(ProgBData):-1:index+1
                    if ishandle(ProgBData(k).AxHndl)
                        delete(ProgBData(k).AxHndl);
                        ProgBData(k).AxHndl = [];
                        ProgBData(k).PatchHndl = [];
                        ProgBData(k).EstTimeHndl =[];
                    end;
%                     if ishandle(ProgBData(k).EstTimeHndl)
%                         delete(ProgBData(k).EstTimeHndl)
%                     end;
                    ProgBData(k)=[];
                end;
            end;
            if index<=0 %| isempty(ProgBData(index).AxHndl)
                if ishandle(GeneralData.ProgB)
                    set(GeneralData.ProgB,'visible','off');
                else 
                    remInp{1}= 'close';
                    usedInputs = 0;
                end;
            else
                remInp{end+1} = 'reposition'; %#ok: expand inside loop: initially we dont know how large it will get.
            end;
          end; % end has gui
        case 'reposition'
          if GeneralData.hasGUI  
            if isempty(GeneralData.ProgB) || ~ishandle(GeneralData.ProgB)
                warning('PROGRESSBAR:notfound','No progressbar, so cannot modify it.');
                return;
            end;
            pointsPerPixel = 72/get(0,'ScreenPixelsPerInch');
            pos = get(GeneralData.ProgB,'Position');   
            
            width = 400 * pointsPerPixel; 
            height = (45+60*index) * pointsPerPixel;
            if 0
                pos(1:2) = pos(1:2) - ([width height]-pos(3:4))/2;
                pos(3:4) = [width height];   
            else
                pos(2) = pos(2) - (height-pos(4))/2;
                pos(4) = height;
            end;
            set(GeneralData.ProgB,'Position',pos);        

            axPosStp = [0 -60*pointsPerPixel/height 0 0];
            axPos = [.1 1-5*pointsPerPixel/height .8 20* pointsPerPixel/height];
            if ~strcmp(get(GeneralData.ProgB,'visible'),'on')
                set(GeneralData.ProgB,'visible','on');
            end;
            for k=1:index
                if ishandle(ProgBData(k).AxHndl)
                    set(ProgBData(k).AxHndl,'Position',axPos+k*axPosStp);
                end;
            end;
            set(GeneralData.ProgB,'name',ProgBData(index).title);
            drawnow;
          end; %end has gui
          usedInputs = 1;
        case 'position'
            if length(ProgBData)>=1 && ishandle(GeneralData.ProgB)
                pos = ProgBData(index).UpdHist(1,1);    
            else
                pos = nan;
            end;
            usedInputs = 1;
        case 'close'
            if ishandle(GeneralData.ProgB)
%                 close(GeneralData.ProgB);
                delete(GeneralData.ProgB);
            end;
            GeneralData = [];
            ProgBData =[];
%             usedInputs = 1;
            return;
        case 'mintimeinterval'
            ProgBData(index).minTimeInterval =  remInp{2} / 86400;
            GeneralData.NextUpdate  = ProgBData(index).UpdHist(2,1) + ProgBData(index).minTimeInterval;
            if prod(size(GeneralData.timer))>1 %#ok: bug in numel of timer objects
                didruntimer = strcmp(get(GeneralData.timer(2),'Running'),'on');
                if didruntimer
                    stop(GeneralData.timer(2));
                end;
                set(GeneralData.timer(2),'Period',max(2,remInp{2}));
                if didruntimer
                    start(GeneralData.timer(2));
                end;
            end;
        case 'cleartimeestimation'
            ProgBData(index).UpdHist = [ProgBData(index).UpdHist(1,1);now];
            usedInputs =1;
        case 'esttimeleft'
            if strcmpi(remInp{2},'on')
                if GeneralData.NeedDraw>0
                    if ishandle(ProgBData(index).EstTimeHndl) % treat empty correctely.
                    else
                        ProgBData(index).EstTimeHndl = -1;
                    end;
                else
                    remInp{1} = 'MakeTimeLeftText';
                    remInp{2} = index;
                    usedInputs = 0;
                end;
                ProgBData(index).nUpdHist = 40;
            elseif strcmpi(remInp{2},'off')
                if ishandle(ProgBData(index).EstTimeHndl)
                    delete(ProgBData(index).EstTimeHndl);
                end;
                ProgBData(index).EstTimeHndl =[];
                ProgBData(index).nUpdHist = 0;
            else
                error('set EstTimeLeft to [on/off]');
            end;
        case 'maketimelefttext'
            index = remInp{2};
            if ishandle(ProgBData(index).EstTimeHndl) % takes care of empty handles.
                set(ProgBData(index).EstTimeHndl,'String',''); % clear time left.
            elseif GeneralData.hasGUI
                ProgBData(index).EstTimeHndl = text(0.01,.5,'','parent',ProgBData(index).AxHndl);
            end;
%             set(ProgBData(index).EstTimeHndl,'String','test')
%         case 'title'
%             if ishandle(ProgBData(index).AxHndl)
%                 set(get(ProgBData(index).AxHndl,'title'),'string',remInp{2});
%             end;
        case 'startstatistics'
            count = remInp{2};
            if ~numel(count)==1 || round(count)~=count
                error('number of statistic points should be integer scalar');
            end
            if exist('imaqmem') %#ok<EXIST>
                val = imaqmem;
            elseif exist('memory') %#ok<EXIST>
                val = memory;
            end;
            val.position = ProgBData(index).UpdHist(1,1);
            val.time = ProgBData(index).UpdHist(2,1);
            names = fieldnames(val);
            % allocate result matrices:
            for k=1:length(names)
                eval(['val.' names{k} '(1,count)=0;']); % required for versions < R13
%                 val.(names{k})(count)=0;
            end;
            ProgBData(index).Statistics = val;
            ProgBData(index).StatisticsCount = 1;
        case 'getstatistics'
            pos = ProgBData(index).Statistics;
            if ~isempty(pos)
                names = fieldnames(pos);
                for k=1:length(names)
                    eval(['pos.' names{k} '(ProgBData(index).StatisticsCount+1:end)=[];']); % required for versions < R13
%                     pos.(names{k})(ProgBData(index).StatisticsCount+1:end)=[];
                end;
            end;
            usedInputs = 1;
        case 'debugbutton'
          if GeneralData.hasGUI && ~isdeployed
            if isempty(GeneralData.DebugButton) || ~any(ishandle(GeneralData.DebugButton)) %use any to evaluate empty to 0.
                if ishandle(GeneralData.ProgB) % deal with empty handle
                else
                    progressbar('redraw');
                end;
                GeneralData.DebugButton = uicontrol(GeneralData.ProgB,'Style', 'pushbutton', 'String', 'Debug',...
                    'Position', [30 3 60 20],'Callback', 'progressbar(''DebugStop'');');
            end;
            if strcmpi(remInp{2},'show')% || strcmpi(remInp{2},'immediatestop')
                set(GeneralData.DebugButton,'visible','on');
%                 if strcmpi(remInp{2},'immediatestop')
%                     GeneralData.doDebugStop=-1;
%                 end;
            else
                set(GeneralData.DebugButton,'visible','off');
            end;
          end;
        case 'usecallcount'
            if ~islogical(remInp{2}) || numel(remInp{2})~=1
                warning('PROGRESSBAR:invalidvalue','invalid value for option ''UseCallCount''');
            else
                ProgBData(index).UseCallCnt = remInp{2};
            end;
        case 'debugstop'
            GeneralData.doDebugStop = 1;
            usedInputs = 1;
        case 'hasgui'
            switch lower(remInp{2})
                case 'yes'
                    GeneralData.hasGUI = true;
                    remInp{end+1} = 'redraw'; %#ok: I want to expand.
                case 'no'
                    GeneralData.hasGUI = false;
                otherwise
                    warning('PROGRESSBAR:invalidvalue','Unknown value for ''hasgui''');
            end;
        otherwise
            chr = remInp{1};
            if ~ischar(chr)
                chr = [num2str(size(chr)) ' ' class(chr)];
            end;
            warning('PROGRESSBAR:invalidargument',['invallid argument in progressbar (' chr ')']);
            
%             error('invallid argument in progressbar');
    end;
    remInp = remInp(1+usedInputs:end);
  end; % end while inputs;
else
    if isempty(ProgBData)
        error('initialize progressbar first');
    end;
    if nargin>=3 
        if numel(varargin{3})~=1 || ~isnumeric(varargin{3}) || varargin{3}~=round(varargin{3})
            error('Third input argument should be integer scalar');
        end;
        index = varargin{3};
    else
        index = GeneralData.CurrIndex;
    end;
    ProgBData(index).CallCount  = ProgBData(index).CallCount +1;
    if GeneralData.NextUpdate > now
%         ProgBData(index).LastUpd(2) + ProgBData(index).minTimeInterval 
        % no update needed yet 
        return;
    end;
    
    if ProgBData(index).UseCallCnt
        newPos = ProgBData(index).CallCount;
    else
        newPos = varargin{1}(:);
    end
    cur = ProgBData(index);
    ststp = cur.StartStop;
    pos = (newPos-ststp(:,1))./(ststp(:,2)-ststp(:,1));
    
    if size(ststp,1)>1
        fpos = max(pos);
    else
        fpos = sum(pos);
    end;
    % update history; DO
    % - Prepend current. 
    % - When new position is equal to the previous position, remove the last entered entry. 
    % - When old history is longer than nUpdHist, truncate history.
    cur.UpdHist = [[fpos;now] cur.UpdHist(:,(1+(fpos==cur.UpdHist(1,1))):min(size(cur.UpdHist,2),ProgBData(index).nUpdHist))];
	GeneralData.NextUpdate  = cur.UpdHist(2,1) + cur.minTimeInterval;
    
    if GeneralData.ParallelState(1)==2
        ProgBData(index) = cur;   
        % Worker progress store
        fseek(GeneralData.ParallelState(2),0,'bof');
        tmp = nan(index+1,1,'single');
        for k=1:index;
            tmp(k) = ProgBData(k).UpdHist(1,1);
        end;
%         [tmp'  fpos] % debug
        fwrite(GeneralData.ParallelState(2),tmp,'single');
        return;
    end;
    
    if GeneralData.NeedDraw>0
        progressbar('Redraw');
    end;
    if isempty(GeneralData.ProgB)
        % figure closed 
        return;
    end;
    if ~ishandle(GeneralData.ProgB)
        warning('PROGRESSBAR:closed','Progressbar closed, no further progress displayed until a new level is created.');
        GeneralData.ProgB = [];
        return;
    end;

   
    if GeneralData.hasGUI
        if ~ishandle(cur.PatchHndl)
            warning('PROGRESSBAR:lostbar','Do not draw to the progressbar');
            GeneralData.NeedDraw = max(GeneralData.NeedDraw,1);
            progressbar('redraw');
            cur = ProgBData(index);
        end;
        pos = pos';
        if size(ststp,1)>1
            yh = (0:size(pos,2))/size(pos,2);
            s1 = ceil(size(pos,2)/2);
            s2 = ceil((size(pos,2)-1)/2);
            set( cur.PatchHndl(1),'XData',[zeros(1,s1); pos(1:2:end); pos(1:2:end); zeros(1,s1)],'YData',[yh(1:2:end-1);yh(1:2:end-1);yh(2:2:end);yh(2:2:end)]);
            set( cur.PatchHndl(2),'XData',[zeros(1,s2); pos(2:2:end); pos(2:2:end); zeros(1,s2)],'YData',[yh(2:2:end-1);yh(2:2:end-1);yh(3:2:end);yh(3:2:end)]);
        else
            pos = cumsum([0 pos]);
            set( cur.PatchHndl(1),'XData',[pos(1:2:end-1); pos(2:2:end); pos(2:2:end); pos(1:2:end-1)],'YData',[0;0;1;1]*ones(1,ceil((size(pos,2)-1)/2)));
            set( cur.PatchHndl(2),'XData',[pos(2:2:end-1); pos(3:2:end); pos(3:2:end); pos(2:2:end-1)],'YData',[0;0;1;1]*ones(1,ceil((size(pos,2))/2-1)));
            pos = pos(2:end);
        end;

        if nargin>1 && ~isempty(varargin{2})
            value = varargin{2};
            set(get(cur.AxHndl,'title'),'string',value);
        end;
        
        if ishandle(cur.EstTimeHndl)
            if size(ststp,1)>1
                pos = max(pos);
            else
                pos = pos(end);
            end;
            progress = cur.UpdHist(1,[1 end])*[1;-1];

            if progress <= eps 
                timeLeft = nan;
            else
                DC = diff(cur.UpdHist,1,2)';
                prgr = DC(:,1)\DC(:,2);
                resid = DC(:,2)-DC(:,1)*prgr;
                resVar = (resid'*resid)./max(eps,size(DC,1)-1);
                if size(DC,1)>3
                    % do weighting of exceptionally large intervals.
                    % (to make it much more robust to large time steps, as you
                    %  get by hibernating)
    %                 oldResVar = resVar*6;
    %                 while oldResVar>resVar*5
                        W = min(1 , (9.*resVar)^1.5.*abs(resid).^(-3));
                        scaler = inv(DC(:,1)'*(W.*DC(:,1)));
                        prgr = scaler*(DC(:,1)'*(W.*DC(:,2)));
                        resid = DC(:,2)-DC(:,1)*prgr;
    %                     oldResVar = resVar;
                        resVar = (resid'*(W.*resid))./(size(DC,1)-1);
    %                 end;
                else
                    scaler = 1./(DC(:,1)'*DC(:,1));
                    if size(DC,1)==1
                        resVar = inf;
                    end;
                end;
                timeLeftDelt =  sqrt(resVar * scaler) * (1-fpos) *86400 * 2.5; 
                timeLeft = (1-fpos) .* prgr .* 86400;

    %             X = [ones(size(cur.Times,1),1) cur.Times(:,1)];
    %             prgr = X\ cur.Times(:,2);
    %             timeLeft = (ststp(2)-newPos) .* prgr(2) .* 86400;
    %             if size(X,1)>=3
    %                 resNrm = sum((cur.Times(:,2) - X*prgr).^2) ./ (size(X,1)-2);
    %                 % cov(prgr) = inv(X'*X)
    %                 % cov of a prediction at x=X0: X0*cov(prgr)*X0'
    %                 % don't use inv(X'*X) to avoid problems with bad scaling. (inv(X'*X) * a = (X\(X'\a)) )
    %                 % scale by 2 to get ~96% confidence interval; 
    %                 % scale by 86400 to get from days to seconds.
    %                 timeLeftDelt = 2*sqrt( resNrm * ([1 ststp(2)] * (X \(X'\[1 ststp(2)]'))))*86400;
    %             else
    %                 timeLeftDelt = nan;
    %             end;
    %             timeLeft =  (ststp(2)-newPos)* (cur.Times(end,2)-cur.Times(1,2)) / progress *86400;
            end;
            if timeLeft>=0 && timeLeft<1e9
                % make timeleft part:
                if timeLeft<60
                    timeLeftStr = sprintf('%4.1fs',timeLeft) ;
                elseif timeLeft<3600
                    timeLeftStr = sprintf('%2.0fm:%4.1fs',floor(timeLeft/60),rem(timeLeft,60));
                else
                    if timeLeftDelt>=100
                        timeLeftStr = sprintf('%dh:%2.0fm',floor(timeLeft/3600),rem(floor(timeLeft/60),60));
                    else
                        timeLeftStr = sprintf('%dh:%2.0fm:%2.0fs',floor(timeLeft/3600),rem(floor(timeLeft/60),60),rem(timeLeft,60));
                    end;
                end;
                % make error part:
                if isfinite(timeLeftDelt) 
                    if timeLeftDelt>=100
                        timeLeftStr = [timeLeftStr sprintf(' +/- %2.0fm',round(timeLeftDelt/60))];
                    elseif timeLeftDelt>=2
                        timeLeftStr = [timeLeftStr ' +/- ' num2str(round(timeLeftDelt)) 's'];
                    else
                        timeLeftStr = [timeLeftStr sprintf(' +/- %2.1fs',timeLeftDelt)];
                    end;
                else
                    timeLeftStr = [timeLeftStr ' +/- ?'];
                end;
            else
                timeLeftStr = 'unknown';
            end;
            set(cur.EstTimeHndl,'String',['Time left: ' timeLeftStr]);
        end; % end estTimeLeft
        drawnow;
    end; % end hasGUI

    if ~isempty(cur.Statistics)
        count = cur.StatisticsCount;
        if count/length(cur.Statistics.time) <= fpos
            count = min(count+1,length(cur.Statistics.time));
            cur.StatisticsCount = count;
            if exist('imaqmem') %#ok<EXIST>
                val = imaqmem;
            elseif exist('memory') %#ok<EXIST>
                val = memory;
            end;
            val.position = varargin{1};
            val.time = cur.UpdHist(2,1);
            names = fieldnames(val);
            for k=1:length(names)
%                 eval(['cur.Statistics.' names{k} '(1:count,:)=val.' names{k} ';']);% required for versions < R13
                cur.Statistics.(names{k})(1:numel(val.(names{k})),count) =val.(names{k})(:);
            end;
        end;
    end; 
    ProgBData(index) = cur;    
end;

if GeneralData.doDebugStop
    [st]=dbstack;
    if length(st)<=1 || ~strcmp(st(1).name,st(min(length(st),2)).name) % when inside another progressbar function call, step to the outer version.
        GeneralData.doDebugStop =0;
        str=['dbstop in ''' st(1).name ''' at ' num2str(st(1).line+6) ];% last line of progressbar
        eval(str);
% PRESS F10 (dbstep) twice to step out to the calling function.     Happy debugging!
        eval(['dbclear in ''' st(1).name ''' at ' num2str(st(1).line+6)]); return;
    end;
end;
