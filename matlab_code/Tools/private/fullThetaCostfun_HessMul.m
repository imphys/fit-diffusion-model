function [Hx ] = fullThetaCostfun_HessMul(hessinfo, x, opts)
% Multiply with hessian of full cost function
% d2 LL/d theta d theta * x= d2 LL/ d img d img  * dimg/ d theta * dimg/ d theta * x
%                            + d LL/d img * d2 img / d theta d theta * x
%                            + d2 regul/ d theta d theta *x
%                            + d2 prior/d theta d theta *x
if size(x,2)>1
    Hx = cell(1,size(x,2));
    for k=1:size(x,2)
        Hx{k} = fullThetaCostfun_HessMul(hessinfo, x(:,k), opts);
    end;
    Hx = [Hx{:}];
    return;
end;
x = reshape(x, opts.numParam ,[]);
gfx = opts.function_jacMul( hessinfo.gfun  , x );
gfHx = projectAllAndLL_HessMul(hessinfo.projGradHess, gfx, opts);
if ~isempty(hessinfo.gfun)
    gfHx = reshape(gfHx, size(hessinfo.gfun,1),[]);
end;
Hx = opts.function_jacMulAdj( hessinfo.gfun , gfHx ) ;
if ~isempty(hessinfo.HPfun)
    Hx = Hx + permute( sum( bsxfun(@times, hessinfo.HPfun , x ), 1) ,  [3 2 1]) ;
end;
Hx = Hx(:);
if ~isempty(hessinfo.Hregul)
    Hxr = opts.spatialRegularizer{2}(hessinfo.Hregul, x );
    Hx = Hx+Hxr(:);
end;