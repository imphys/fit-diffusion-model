function [Hx ] = fullThetaAndProjectParCostfun_HessMul(hessinfo, x, opts)
% Multiply with hessian of full cost function
% d2 LL/d theta d theta * x= d2 LL/ d img d img  * dimg/ d theta * dimg/ d theta * x
%                            + d LL/d img * d2 img / d theta d theta * x
%                            + d2 regul/ d theta d theta *x
%                            + d2 prior/d theta d theta *x
if size(x,2)>1
    Hx = cell(1,size(x,2));
    for k=1:size(x,2)
        Hx{k} = fullThetaAndProjectParCostfun_HessMul(hessinfo, x(:,k), opts);
    end;
    Hx = [Hx{:}];
    return;
end;

sztheta = [opts.numParam opts.spatialSize];
curst = prod(sztheta);
x_theta = reshape(x(1:curst),sztheta);
x_projectParameters = cell( 1 , numel(opts.data_in) );
for k=1:numel(x_projectParameters)
    szpp = size(opts.projectParameters{k});
    cured = curst + prod(szpp);
    x_projectParameters{k} = reshape( x(curst+1:cured),szpp);
    curst = cured;
end;

gfx = opts.function_jacMul( hessinfo.gfun  , x_theta );
[gfHx , gProjHx] = projectAllAndLL_HessMul(hessinfo.projGradHess, gfx, opts, x_projectParameters);
if ~isempty(hessinfo.gfun)
    gfHx = reshape(gfHx, size(hessinfo.gfun,1),[]);
end;
Hx = opts.function_jacMulAdj( hessinfo.gfun , gfHx ) ;
if ~isempty(hessinfo.HPfun)
    Hx = Hx + permute( sum( bsxfun(@times, hessinfo.HPfun , x_theta ), 1) ,  [3 2 1]) ;
end;
Hx = Hx(:);
if ~isempty(hessinfo.Hregul)
    Hxr = opts.spatialRegularizer{2}(hessinfo.Hregul, x_theta );
    Hx = Hx+Hxr(:);
end;
for k=1:numel(x_projectParameters)
    if size(opts.projectParameterPrior,1)==1
        projectParameterPrior_k = opts.projectParameterPrior;
    else
        projectParameterPrior_k = opts.projectParameterPrior(k,:);
    end;
    if ~isempty(projectParameterPrior_k{1})
        gProjHx{k} = gProjHx{k} - projectParameterPrior_k{2}( hessinfo.projectprior{k} , x_projectParameters{k});
    end;
end;

Hx = vertcat(Hx,gProjHx{:});
