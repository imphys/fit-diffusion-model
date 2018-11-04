function [Px, gradifo] = project1( x, projectParameters, project, option)
% call: [Px, gradifo] = project1( x, projectParameters, project, 1)
% or  : [Px, gradifo] = project1( x, gradifo          , project, option); % with option = 2, 3, 4, 5
% Option specifies which 'column' of project is used.
%
% Function intended to 'project' one image to the data domain or the adjoint transform.
% Correctly handles cell array project's
% 
% INPUTS:
%  x : input to be projected
%  projectParameters : the projectParameters of the current row of project
%  gradifo : output of previous project{1} of this row.
%  project : one row of the entire project array; the row corresponding to x.
% 
% Created by Dirk Poot, Erasmus MC, 2-2-2012
projgrad = cell(1,min(1,nargout-1));
if iscell(project{1})
    if ismember(option,[1 2 4])
        Px = cell(size(project{1}));
    else
        Px = 0;
    end;
    if nargout>=2
        gradifo = cell(size(project{1}));
    end;
    for k=1 : numel(project{1})
        if iscell( project{option} )
            projfun = project{option}{k};
        else
            projfun = project{option};
        end;
        if ismember(option,[1 2 4])
            [Px{k}, projgrad{:}] = projfun( x   , projectParameters{k} );
        else
            [Pxk, projgrad{:}] = projfun( x{k}, projectParameters{k} );
            Px = Px + Pxk;
        end;
        if nargout>=2
            gradifo{k} = projgrad{1};
        end;
    end;
else
    [Px, projgrad{:}] = project{option}(   x , projectParameters );
    if nargout>=2
        gradifo = projgrad{1};
    end;
end;