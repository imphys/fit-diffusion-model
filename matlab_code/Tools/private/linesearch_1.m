function alpha = linesearch_1(funfcn, f0, g0, f1, g1, maxpolyorder, bracket)
% alpha = linesearch_1(funfcn, f0, g0, f1, g1, maxpolyorder [, bracket])
% Custom line search algorithm. 
% uses, [f0, g0] = funfcn( 0 ) 
%       [f1, g1] = funfcn( 1 ) 
%
% interpolates points within bracket with polynome of at most order 'maxpolyorder'
% and locates extrema points within bracket. Default initial bracket = [0 3]
%
% Created by Dirk Poot, Erasmus MC. 

sigma = .1; % acceptable gradient fraction (this is the maximum allowed gradient for an acceptable point, as fraction of g0)
rho = .1; % acceptable function value improvement (function value should be < f0 + rho * g0 * alpha)

% Find a bracket of acceptable points
x = [0 f0 g0;
     1 f1 g1];
if nargin<7
    bracket = [0 3] ; 
end;
satisfied = false;
iteridx = 1;
while ~satisfied
    selx = (x(:,1)>=bracket(1)) & (x(:,1)<=bracket(2));
    if numel(selx)*2-1>maxpolyorder
        % reduce bracket.
        [dum, b] = min(x(:,2)); % locate current global minimum.
        [dum,ind]  = sort(x(:,1));
        fi = find(ind==b); 
        if ((x(b,3)>0) && (fi>1)) || (fi==size(x,1))
            selx = ind([-1 0]+fi);
        else
            selx = ind([0 1]+fi);
        end;
        if x(selx(1),3)>0 || x(selx(end),3)<0
            warning('One or both edge points of interval have unexpected gradient. Minimum may be outside bracket.');
        end;
        bracket = x(selx,1)';
    end;
    sx = x(selx,:);
    alpha_c = mean(sx(:,1)); % center value around which polynomial is expanded.
    rx = sx(:,1)-alpha_c;
    polord = min(maxpolyorder, size(sx,1)*2-1);
    X = zeros(size(sx,1)*2,polord+1);
    for k=0:polord
        X(:,k+1) = [rx.^(polord-k);(polord-k).*rx.^max(0, polord-k-1 )];
    end;
    poly = X\[sx(:,2);sx(:,3)]; % matrix division (or least squares fit)
    dpoly = poly(1:end-1)'.*(polord:-1:1);
    r = roots(dpoly);
    r(abs(imag(r)) > eps) = []; % remove complex roots.
    r( ((r+alpha_c) < bracket(1)) | ((r+alpha_c) > bracket(2)) ) = [];
    if (iteridx==1)
        r(end+1) = bracket(2)- alpha_c; % also consider edge point.
    end;
    if isempty(r)
        % function has no extrema (within bracket) ????
        % return alpha with lowest value:
        break;
    end;
    fextrema_pred = poly(1);
    for k=2:numel(poly)
        fextrema_pred = fextrema_pred .* r + poly(k);
    end;
    [ropt , ind] = min(fextrema_pred); % locate lowest extemum.
    newtestalpha = r( ind ) + alpha_c;
    if ropt > min(x(:,2))
        % Cant find better candidate within interval.
        break;
    end;

    [ft gt] = funfcn( newtestalpha );
    x(end+1,:) = [newtestalpha ft gt];
    
    satisfied = (iteridx>3) || ((abs(gt) <= -sigma*g0) && (ft <= f0+ newtestalpha * g0 * rho )); % test #iterations and Wolfe conditions.
    
    iteridx = iteridx +1;
    
end;
% return lowest value found:
[dum,b]=min(x(:,2));
alpha = x(b,1);