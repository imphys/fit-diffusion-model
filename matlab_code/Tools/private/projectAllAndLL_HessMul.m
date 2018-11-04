function [Hx, projHx] = projectAllAndLL_HessMul(hessinfo, x_theta, opts, x_projectParameters)
% [Hx] = projectAllAndLL_HessMul(hessinfo, x_theta, opts, x_projectParameters)
% Multiply x_theta, which is size(data) with hessian of likelihood function, ignoring hessian of project.
% d2 LL/d img  d img * x_theta= d2 LL/ d A d A  * dA/ d img * dA/ d img * x_theta
%
% Created by Dirk Poot, Erasmus MC, 2-2-2012

hasProjPar = nargin>=4 && ~isempty( x_projectParameters );
if hasProjPar 
    projHx = cell(size(x_projectParameters));
end;
szData =  opts.spatialSize;
sel = repmat({':'},1,numel(szData)+1);
xp = permute( reshape(x_theta,[size(x_theta,1) szData]) , [2:numel(szData)+1 1]);
Hx = cell(size(opts.project,1),1);
for k=1:size(opts.project,1)
    sel{end} = opts.projectSelectSource{k};
    [Px , hessinfo.projectGradinfo{k}]= project1( xp(sel{:}) , hessinfo.projectGradinfo{k} , opts.project(k,:), 2 );

    if hasProjPar
        [PxProj, hessinfo.projectGradinfo{k} ] = project1( x_projectParameters{k} , hessinfo.projectGradinfo{k} , opts.project(k,:), 4 );
        Px = Px+PxProj;
    end;    
    
    if iscell(Px)
        for k2 = 1:numel(Px)
            Px{k2} = -Px{k2} .* hessinfo.lgpdfhess{k}{k2};
        end;
    else
        Px = -Px .* hessinfo.lgpdfhess{k};
    end;
    Hx{k} = project1(Px  , hessinfo.projectGradinfo{k} , opts.project(k,:), 3 );
    if numel(sel{end})==1
        Hx{k} = reshape(Hx{k},1,[]);
    else
        Hx{k} = reshape(Hx{k},[], numel(sel{end})).';
    end
    if hasProjPar  
        [projHx{k}, hessinfo.projectGradinfo{k} ] = project1( Px , hessinfo.projectGradinfo{k} , opts.project(k,:), 5 );
    end;
end;

if isequal([opts.projectSelectSource{:}],1:size(x_theta,1))
    Hx = vertcat(Hx{:});
else
    tmp = zeros(size(x_theta));
    for k=1:numel(Hx)
        tmp(opts.projectSelectSource{k},:) = tmp(opts.projectSelectSource{k},:) + Hx{k};
    end;
    Hx = tmp;
end;