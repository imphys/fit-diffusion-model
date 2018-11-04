function [f, g, H] = totalVariationVecRegularizer(x , spacing , weights, offset , mask, outputhessmul, derivdir, cyclic)
% [f, g, H] = totalVariationVecRegularizer(x , spacing , weights [, offset, [mask, [outputhessmul , [derivdir, cyclic]]]])
% Computes the total variation cost, gradient and hessian for an image x 
% in which each voxel is vector valued.
% 
% INPUTS:
% x       : N-dimensional image, the first dimension is the 'vector' dimension.
%           So the vector of values of the first voxel is x(:,1,1,..).
%           n = size(x,1) : the number of elements in each parameter vector.
% spacing : N element vector that specifies the spacing between the voxels 
%           in the corresponding dimensions. The first element is not used as
%           it has no corresponding spatial dimension. The remaining elements 
%           are used to correct the gradient computations for voxel spacing.
% weights : n element vector or [n x n] positive definite symmetric matrix 
%           with which the gradient norm is computed. A vector weights 
%           is slightly faster, but equivalent to using diag(weights).
%           The 2-norm of a gradient vector ga = (x(:,a)-x(:,b)) is computed by
%           Fab = ga'*weights*ga.
% offset  : Default = .1. A scalar offset >0 that smooths the total variation 
%           function so optimization with standard routines is easier.
% mask    : optional, logical matrix that selects the voxels of which the gradient 
%           (and hessian) should be computed. Default = [] => all voxels.
%           size(mask)==[1 size(x,[2..N])]
% outputhessmul : integer that if 0 (=false and default) indicates that the third output 
%           should be a sparse matrix with the explicit hessian.
%           if 1 (=true) a hessian multiply structure, so hessmulfun(H, x) multiplies the hessian with x
%           This form uses less memory and typically is faster (when optimizing). 
%           Don't use this form to explicitly get the full hessian (since that's slow and consumes a lot of memory).
%           if 2 : return hessmulfun (single output)
% derivdir : 1 (default) : forward and backward derivatives
%            0           : only forward derivatives.
% cyclic : default false; if true makes x cyclic
%
% OUTPUTS:
% f   : Total variation penalty value: sum_i sqrt( offset^2 + sum_n |x(:,i)-x(:,n)|^2_weights ) -offset
%       i : all voxels
%       n : direct neighbors of i (or forward neighbors if derivdir==0)
%       Note that both forward and backward derivatives are used, this in contrast
%       to most standard total variation penalties. While this introduces a 
%       very slight smoothing effect, it makes the total variation symmetric.
%       So the total variation of x and x flipped over any of the spatial dimensions is equal.
% g   : gradient of f with respect to x. (size(g)==size(x))
% H   : Hessian of f with respect to x. (size(H) == [numel(x) numel(x)], sparse matrix)
%       or, if outputhessmul==true: structure so that hessmulfun( H, x) == Hessian * x;
%       with hessmulfun = totalVariationVecRegularizer([] , [], [], [], [], 2);
%
% Created by Dirk Poot, Erasmus MC, 23-3-2011

if nargin<4 || isempty(offset)
    offset = .1; % makes function smoother => easier to optimize.
end;
if nargin < 7 || isempty(derivdir)
    derivdir = 1;
end;
if nargin<8 
    cyclic = false;
else
    if cyclic && derivdir~=1
        error('can currently do cyclic when derivdir==1');
    end;
end;
useMexVersion = true; % set to false if you don't want to use the mex file.

if nargin<6 || isempty(outputhessmul)
    outputhessmul = false;
elseif outputhessmul==2
    if useMexVersion
        f = @TVhessmul_mex;
    else
        f = @TVhessmul;
    end;
    return;
end;
if nargout>=3 && ~outputhessmul
    useMexVersion = false;
    warning('TotalVariationVecRegularizer:Slow','Using slower MATLAB code instead of mex routine, since explicit hessian is requested');
end;
if useMexVersion
    argo = cell(1,nargout);
    rspacing2 = 1./spacing.^2;
    if derivdir==1 % both forward and backward derivatives under sqrt.
        if cyclic
            sel = cell(1,ndims(x));
            sel{1} = ':';
            for dim = 2:numel(sel);
                sel{dim} = [size(x,dim) 1:size(x,dim) 1];
            end;
            [argo{:}] = TotalVariation_Vector_cSkipB(x(sel{:}),rspacing2,weights, offset);
            if nargout>=2
                sela = repmat( {':'},1,ndims(x));
                for dim = 2:numel(sel);
                    selr = sela; selw= sela;
                    selr{dim} = [1 size(x,dim)+2];
                    selw{dim} = [size(x,dim)+1 2];
                    argo{2}(selw{:}) = argo{2}(selw{:}) + argo{2}(selr{:});
                    selc = sela;
                    selc{dim} = [2:size(x,dim)+1];
                    argo{2} = argo{2}(selc{:});
                end;
            end;
        else
            [argo{:}] = TotalVariation_Vector_c(x,rspacing2,weights, offset);
        end;
    elseif derivdir==2 % central differences.
        gradvec = [-.5 0 .5];gradvecoffset = -1;
        [argo{:}] = TotalVariation_Vector_gradfilt_fc(x,rspacing2,weights, offset, gradvec, gradvecoffset);
    elseif derivdir==3 % central differences2.
        gradvec = [-.125 .75 0 -.75 .125];gradvecoffset = -2;
        [argo{:}] = TotalVariation_Vector_gradfilt_fc(x,rspacing2,weights, offset, gradvec, gradvecoffset);
    else% only forward derivatives under sqrt.
        [argo{:}] = TotalVariation_Vector_fc(x,rspacing2,weights, offset);
    end;
    f = argo{1};
    if nargout>=2
        g = argo{2};
        if nargin>=5 && ~isempty(mask)
            g = g(mask);
        else
            mask = [];
        end;
        if nargout>=3
            H.img = x;
            H.rspacing2 = rspacing2;
            H.weights = weights;
            H.offset = offset;
            H.H = argo{3};
            H.mask = mask;
            H.derivdir = derivdir;
            H.cyclic = cyclic;
            if ~outputhessmul
                error('mex version cannot output explicit hessian efficiently. set input argument ''outputhessmul'' to true.');
                % actually, TotalVariation_Vector_c can be rebuild to output a sparse hessian matrix, but that's relatively slow. 
            end;
        end;
    end;
    return;
end;

% non mex version:
ndim = ndims(x);
sx = size(x);
sumsq = zeros([1 sx(2:end)]);
indx = repmat({':'},1,ndim);
ismatweight = all(size(weights)==sx([1 1]));
if ~ismatweight
    weights = weights(:)'; % enforce row column.
end;
for dim = 2:ndim
    indxu = indx;
    indxu{dim} = 1:size(x,dim)-1;
    dx = diff(x,1,dim)/spacing(dim);
    dxsz = size(dx);
    if ismatweight
        dx_sq = sum(dx(:,:) .* (weights * dx(:,:)),1);
    else
        dx_sq = weights * dx(:,:).^2;
    end;
    dx_sq = reshape(dx_sq, [1 dxsz(2:end)]);
    sumsq(indxu{:}) = sumsq(indxu{:}) + dx_sq;
    indxu{dim} = 2:size(x,dim);
    sumsq(indxu{:}) = sumsq(indxu{:}) + dx_sq;
end;
% f_beforesum = sqrt( sumsq + offset^2*ones([1 sx(2:end)]));
f_beforesum = sqrt( sumsq + offset^2 );
f = sum( f_beforesum(:) - offset );

if nargout>=2
    hasmask = nargin>=5 && ~isempty(mask);
    % derivative of f w.r.t. x :
    % f = sqrt( sum ( w * dx_spatial ^2  + offset^2)  ) - offset
    % dfdxi = .5 / sqrt( sum ( w * dx_spatial ^2  + offset^2)  ) * w * 2 * dx_spatial * d(dx_spatial )/d(xi)
    %       = ( sum ( w * dx_spatial ^2  + offset^2)  ).^(-.5) * w  * dx_spatial * sgn
    % d(dx_spatial )/d(xi) = 1 for first term and -1 for second term.
    % Using :
    % d( sqrt( f(x) ) )/d(x) = d( (f(x))^(0.5) )/d(x) = 0.5 * f(x)^(-.5) * {d(f)/d(x)}|x

    % x = [a b c]
    %                     f_1        f_2          f_3
    % f = sum( sqrt( w*[(b-a)^2 (c-b)^2+(b-a)^2 (c-b)^2] + offset^2) - offset )
    % dfda =   1/sqrt( w*[(b-a)^2        ] + offset^2 ) * w * (b-a) * (-1)    
    %        + 1/sqrt( w*[(c-b)^2+(b-a)^2] + offset^2 ) * w * (b-a) * (-1)
    % dfdb =   1/sqrt( w*[(b-a)^2        ] + offset^2 ) * w * (b-a) * (1) 
    %          1/sqrt( w*[(c-b)^2+(b-a)^2] + offset^2 ) * w * ((c-b) * (-1) + (b-a)*(1)) 
    %        + 1/sqrt( w*[(c-b)          ] + offset^2 ) * w * (c-b) * -1
   
    rf = 1./ f_beforesum;
    g = zeros(size(x));
    if nargout>=3
        sumwdx = zeros(size(x));
    end;
    if nargout>=3
        linindx = reshape(1:numel(x),size(x));
        wdxc = cell(1,ndim);
    end;
    for dim = 2:ndim
        indxu = indx;
        indxu{dim} = 1:size(x,dim)-1;
        indxv = indx;
        indxv{dim} = 2:size(x,dim);
        if ismatweight
            dx  = diff(x,1,dim);
            szdx= size(dx);
            wdx =               (weights/spacing(dim)^2) * dx(:,:);
            wdx = reshape(wdx, szdx);
        else
            wdx = bsxfun(@times, weights'/spacing(dim)^2  , diff(x,1,dim));
        end;
        tmpdx = bsxfun(@times, wdx, rf(indxu{:}) + rf(indxv{:}));
        g(indxu{:}) = g(indxu{:}) - tmpdx;
        g(indxv{:}) = g(indxv{:}) + tmpdx;

        if nargout>=3
            wdxc{dim} = wdx;
            sumwdx(indxu{:}) = sumwdx(indxu{:}) - wdx;
            sumwdx(indxv{:}) = sumwdx(indxv{:}) + wdx;
        end;
            
    end;
    if hasmask
        g = g(:,mask);
    end;
    if nargout>=3
        if nargin<6 || isempty(outputhessmul)
            outputhessmul = false;
        end;
        rf3 = rf.^3;
        
        if outputhessmul 
            rfiprfn =0 ;
            for dim = 2:ndim
%                 filt = reshape([1 2 1]/spacing(dim)^2 , [ones(1,dim-1) 3 1]);
%                 rfiprfn = rfiprfn + convn(rf, filt, 'same');
                filt = reshape([1 1]/spacing(dim) , [ones(1,dim-1) 2 1]);
%                 filt2 = reshape([1 0 1]/spacing(dim)^2 , [ones(1,dim-1) 3 1]);
                rfiprfn = rfiprfn + convn(convn(rf, filt, 'valid'),filt,'full');%-convn(rf,filt2,'same'); % handle edge cases properly.
            end;
            hessinf.sizex   = size(x);
            hessinf.weights = weights;
            hessinf.rfiprfn = rfiprfn;
        	hessinf.rf      = rf;
        	hessinf.rf3     = rf3;
            hessinf.wdxc    = wdxc;
            hessinf.spacing = spacing;
            hessinf.sumwdx  = sumwdx;
            hessinf.derivdir = derivdir;
            H = hessinf;
%            H = @(y) TVhessmul(hessinf, y);
            return;
        end;
        % Vector weight:
        % f_i = sqrt( offset^2 + sum_j wj * sum_n ((xj_i-xj_n)/sp_n)^2 )       % ignored '- offset' since that is constant
        % i : voxel location
        % j : different element of vector
        % n : neigbors of i;   n~=i
        % d f_i /dxj_i =  1/f_i * (wj * sum_n ((xj_i - xj_n)/sp_n^2) )
        % d f_i /dxj_n = -1/f_i * (wj *       ((xj_i - xj_n)/sp_n^2) )
        % d2 f_i /dxj_i dxk_i =  1/f_i   * (wj * sum_n (  delta_kj   /sp_n^2) )
        %                      - 1/f_i^3 * (wj * sum_n ((xj_i - xj_n)/sp_n^2) ) * (wk * sum_n ((xk_i - xk_n)/sp_n^2) ) 
        % d2 f_i /dxj_i dxk_n =  1/f_i   * (wj *       ( -delta_kj   /sp_n^2) )
        %                      + 1/f_i^3 * (wj * sum_n ((xj_i - xj_n)/sp_n^2) ) * (wk *       ((xk_i - xk_n)/sp_n^2) )
        % d2 f_i /dxj_n dxk_m = +1/f_i   * (wj * delta_kj * delta_nm /sp_n^2  )
        %                      + 1/f_i^3 * (wj *       ((xj_i - xj_n)/sp_n^2) ) * (wk *       ((xk_i - xk_m)/sp_n^2) )

        % Matrix weight:
        % f_i = sqrt( offset^2 + sum_n sum_jk wjk * ((xj_i-xj_n)/sp_n)*((xk_i-xk_n)/sp_n) )       % ignored '- offset' since that is constant
        % i   : voxel location
        % j,k : different element of vector
        % n   : neigbors of i;   n~=i
        %
        % d f_i /dxl_i =  1/2 * 1/f_i * ( sum_n sum_jk wjk * ( ((xj_i - xj_n)/sp_n) * delta_kl /sp_n + delta_jl/sp_n * ((xk_i-xk_n)/sp_n)) )
        %              =  1/f_i * sum_n sum_j  wjl *  (xj_i - xj_n)/sp_n^2            % if wjl==wlj
        % d f_i /dxl_n = -1/(2*f_i) * ( sum_jk wjk *(((xj_i - xj_n)/sp_n  )*delta_kl/sp_n) + delta_jl/sp_n*((xk_i-xk_n)/sp_n) ))
        %              = -1/(2*f_i) * ( sum_j  wjl *  (xj_i - xj_n)/sp_n^2 + sum_j wlj * (xj_i-xj_n)/sp_n^2 )
        %              = -1/f_i *       sum_j  wjl *  (xj_i - xj_n)/sp_n^2            % if wjl==wlj
        %
        % d2 f_i /dxl_i dxk_i =  1/f_i   * ( sum_n sum_j wjl * (   delta_jk  /sp_n^2) )
        %                      - 1/f_i^3 * ( sum_n sum_j wjl * ((xj_i - xj_n)/sp_n^2) ) * ( sum_n sum_j wjk * ((xj_i - xj_n)/sp_n^2) )
        %                     =  1/f_i   * ( sum_n       wkl                 /sp_n^2  ) 
        %                      - 1/f_i^3 * ( sum_n sum_j wjl * ((xj_i - xj_n)/sp_n^2) ) * ( sum_n sum_j wjk * ((xj_i - xj_n)/sp_n^2) )
        % d2 f_i /dxl_i dxk_n =  1/f_i   * (       sum_j wjl * (  -delta_jk  /sp_n^2) )
        %                      + 1/f_i^3 * ( sum_m sum_j wjl * ((xj_i - xj_m)/sp_m^2) ) * (       sum_j wjk * ((xj_i - xj_n)/sp_n^2) )
        %                     =- 1/f_i   * (             wkl                 /sp_n^2  )
        %                      + 1/f_i^3 * ( sum_m sum_j wjl * ((xj_i - xj_m)/sp_m^2) ) * (       sum_j wjk * ((xj_i - xj_n)/sp_n^2) )
        % d2 f_i /dxl_n dxk_m =- 1/f_i   * (       sum_j wjl * (  -delta_jk  /sp_n^2) ) *  delta_nm
        %                      + 1/f_i^3 * (       sum_j wjl *  (xj_i - xj_n)/sp_n^2  ) * (       sum_j wjk * (xj_i - xj_m)/sp_m^2 )
        %                     =  1/f_i   * (             wkl                 /sp_n^2) ) *  delta_nm
        %                      + 1/f_i^3 * (       sum_j wjl *  (xj_i - xj_n)/sp_n^2  ) * (       sum_j wjk * (xj_i - xj_m)/sp_m^2 )

        neighborcnt = sum(2./spacing(2:end).^2)*ones(size(rf));
        rf3sumw = bsxfun(@times, rf3, sumwdx);
        Hi = cell(13,ndim);
        Hj = cell(size(Hi));
        Hv = cell(size(Hi));
        [J, K] = find(ones(size(x,1)));
%         repI = ones(numel(I,1));
        for dim = 2:ndim
            indxu = indx;
            indxu{dim} = [1 size(x,dim)];
            neighborcnt(indxu{:}) = neighborcnt(indxu{:}) - 1/spacing(dim)^2;
            indxu{dim} = 1:size(x,dim)-1;
            indxv = indx;
            indxv{dim} = 2:size(x,dim);        
            wdx = wdxc{dim};
            if hasmask
                % apply mask in computation. 
                indxtmp = indxu;
                indxtmp{dim} = [indxv{dim} 1];
                mu = mask;
                mv = mask(indxtmp{2:end});
                indxtmp{dim} = size(x,dim);
                mu(indxtmp{2:end}) = false;
                mv(indxtmp{2:end}) = false;
                lumask = mu & mv; 
                indxtmp{dim} = [size(x,dim) indxu{dim}];
                lvmask = lumask(indxtmp{2:end}); 
                
                lu = linindx(:,lumask);
                lv = linindx(:,lvmask);
                indxu{2} = lumask;indxu(3:end)=[];
                indxtmp{dim} = [size(x,dim) 1:size(x,dim)-1];
                indxv{2} = lvmask;indxv(3:end)=[];
                indxtmp{dim} = 1:size(x,dim)-1;
                wdx = wdx(:,lumask(indxtmp{2:end}));

                % lud & lvd : masked versions of lu and lv, containing the linear index of u and v, which are separated by 1 step in dimension dim.
                indxtmp{dim} = [size(x,dim) 1:size(x,dim)-1];
                
                indxvd = {':',mu(indxtmp{2:end})}; % specifies i==v, i might be out of mask, but m==n===u should be inside mask
                lud = linindx(:,mu); % specifies m==n==u when i == v

                indxud = {':',mv}; % specifies i==u, i might be out of mask, but m==n===v should be inside mask
                lvd = linindx(:,mv(indxtmp{2:end})); % specifies m==n==v when i == u
            else
                lu = linindx(indxu{:});
                lv = linindx(indxv{:});
                lud = lu;
                lvd = lv;
                indxud = indxu;
                indxvd = indxv;
            end;
            
            % d2 f_i /dxj_i dxk_n
            % n = next in dim :
            if ismatweight
                Hi{1,dim} = lu(J,:);
                Hj{1,dim} = lv(K,:);
            else
                Hi{1,dim} = lu;
                Hj{1,dim} = lv;
            end;
            Hv{1,dim} = bsxfun(@times, rf(indxu{:}) ,- weights(:)/spacing(dim)^2) ; %1/f_i   * (wj *       ( -delta_kj   /sp_n^2) )
            Hi{2,dim} = Hj{1,dim};
            Hj{2,dim} = Hi{1,dim};
            Hv{2,dim} = Hv{1,dim};
            % n = previous in dim :
            Hi{3,dim} = Hi{1,dim};
            Hj{3,dim} = Hj{1,dim};
            Hv{3,dim} = bsxfun(@times, rf(indxv{:}) ,- weights(:)/spacing(dim)^2) ; %1/f_i   * (wj *       ( -delta_kj   /sp_n^2) )
            Hi{4,dim} = Hj{3,dim};
            Hj{4,dim} = Hi{3,dim};
            Hv{4,dim} = Hv{3,dim};
            
            % n = next in dim :
            Hi{5,dim} = lu(J,:);
            Hj{5,dim} = lv(K,:);
            Hv{5,dim} = -reshape(rf3sumw(J,indxu{2:end}),numel(J),[]).*wdx(K,:) ; % + 1/f_i^3 * (wj * sum_n (xj_i - xj_n) ) * (wk *       (xk_i - xk_n) )
            Hi{6,dim} = Hj{5,dim};
            Hj{6,dim} = Hi{5,dim};
            Hv{6,dim} = Hv{5,dim};
            % n = previous in dim :
            Hi{7,dim} = lv(J,:);
            Hj{7,dim} = lu(K,:);
            Hv{7,dim} = reshape(rf3sumw(J,indxv{2:end}),numel(J),[]).*wdx(K,:) ; % + 1/f_i^3 * (wj * sum_n (xj_i - xj_n) ) * (wk *       (xk_i - xk_n) )
            Hi{8,dim} = Hj{7,dim};
            Hj{8,dim} = Hi{7,dim};
            Hv{8,dim} = Hv{7,dim};
            
            % d2 f_i /dxj_n dxk_m 
            if ismatweight
                Hi{ 9,dim} = lvd(J,:);
                Hj{ 9,dim} = lvd(K,:);
            else
                Hi{ 9,dim} = lvd;
                Hj{ 9,dim} = lvd;
            end;
            Hv{ 9,dim} = bsxfun(@times, rf(indxud{:}) , weights(:)/spacing(dim)^2) ; %+1/f_i   * (wj * delta_kj * delta_nm /sp_n^2  )
            if ismatweight
                Hi{10,dim} = lud(J,:);
                Hj{10,dim} = lud(K,:);
            else
                Hi{10,dim} = lud;
                Hj{10,dim} = lud;
            end
            Hv{10,dim} = bsxfun(@times, rf(indxvd{:}) , weights(:)/spacing(dim)^2) ; %+1/f_i   * (wj * delta_kj * delta_nm /sp_n^2  )
            hcrow = 11;
            wdx = wdxc{dim};
            for dim2 = 2:dim 
                wdx_dim2 = wdxc{dim2};
                if dim==dim2
                    prevornextsel = [1 2 4]; % do not do both  n - i - m  and  m - i - n, if dim==dim2
                else
                    prevornextsel = 1:4;
                end;
                for prevornext=prevornextsel
                    % select if i = u or v:
                    ioffsetindim = mod(prevornext,2)~=1;
                    if ioffsetindim==0
                        % i is lower than n in dim
                        indxi = indxud;
                        indxn = indxvd;
                    else
                        % i is higher than n in dim
                        indxi = indxvd;
                        indxn = indxud;
                    end;
                    
                    ioffsetindim2 = prevornext>2;
                    % create selection of m and n
                    if hasmask
                        % both m and n should be within mask, i should be valid, but not nececarily within mask.
                        mmask = mask;
                        indxtmp = repmat({':'},1,ndims(x)-1);
                        if dim~=dim2
                            if ioffsetindim==0
                                % n is higher in dim
                                % so first element of n has no corresponding i, but in dim m corresponds to i, so last element of m has no corresponding n 
                                indxtmp{dim-1} = size(x,dim);
                            else
                                % n is lower in dim
                                % so last element has no corresponding i, but in dim m corresponds to i, so first element of m has no corresponding n 
                                indxtmp{dim-1} = 1;
                            end;
                            mmask(indxtmp{:}) = false;
                        end;
                        indxtmp{dim-1} = ':';
                        if ioffsetindim2==0
                            % m is higher in dim2
                            % so first element has no corresponding i
                            % add extra space when n is lower than i and m moves in same dimension
                            indxtmp{dim2-1} = 1 + (0:double(((dim==dim2)&&(ioffsetindim==1)))); 
                        else
                            % m is lower in dim2
                            % so last element has no corresponding i
                            % add extra space when n is higher than i and m moves in same dimension
                            indxtmp{dim2-1} = size(x,dim2)- (0:double(((dim==dim2)&&(ioffsetindim==0))));
                        end;
                        mmask(indxtmp{:}) = false;
                        mmask = find(mmask);
                        cumszmask = cumprod([1 size(mask)]);
                        % if ioffsetindim2==0 : m>i in dim2, so step backward in dim2 to move from m to i
                        % if ioffsetindim==0 : n>i in dim, so step forward in dim to move from i to n
                        stepmTon = (1-2*ioffsetindim) * cumszmask(dim-1) + (2*ioffsetindim2-1) * cumszmask(dim2-1);
                        mmask = mmask(mask(mmask + stepmTon)); % keep only those for which also n is inside mask.
                        nmask = mmask + stepmTon;
                        imask = mmask + (2*ioffsetindim2-1) * cumszmask(dim2-1);

                        indxm = indxi;
                        indxm{2} = mmask;
                        indxi{2} = imask;
                        indxn{2} = nmask;

                        wsel_n = cell(1,ndim-1);
                        [ wsel_n{:} ] = ind2sub(size(mask), nmask);
                        if ioffsetindim==0
                            wsel_n{dim-1} = wsel_n{dim-1}-1;
                        end;
                        szn = size(wdx);
                        nmask_wdx = sub2ind(szn(2:end), wsel_n{:});
                        wdx_n = wdx(J , nmask_wdx);
                        
                        wsel_m = cell(1,ndim-1);
                        [ wsel_m{:} ] = ind2sub(size(mask), mmask);
                        if ioffsetindim2==0
                            wsel_m{dim2-1} = wsel_m{dim2-1}-1;
                        end;
                        szm = size(wdx_dim2);
                        mmask_wdx = sub2ind(szm(2:end), wsel_m{:});
                        wdx_m = wdx_dim2(K,mmask_wdx);
                    else
                        % expand selection in dim2 if 'select all':
                        if dim2<=numel(indxi) && isequal(indxi{dim2},':')
                            indxi{dim2} = 1:size(x,dim2);
                            indxn{dim2} = 1:size(x,dim2);
                        end;
                        indxm = indxi;
                        wsel_n = repmat({':'},1,ndims(x));
                        wsel_m = repmat({':'},1,ndims(x));
                        % create selector for i, m and n:
                        if dim==dim2 && (prevornext==1 || prevornext==4)
                            indxm{dim2} = indxn{dim2};
                        else
                            indxm{dim2} = indxm{dim2}(2-ioffsetindim2:end  -ioffsetindim2);
                            indxi{dim2} = indxi{dim2}(1+ioffsetindim2:end-1+ioffsetindim2);
                            indxn{dim2} = indxn{dim2}(1+ioffsetindim2:end-1+ioffsetindim2);
                            wsel_n{dim2} = 1+ioffsetindim2:size(wdx,dim2)-1+ioffsetindim2;
                            wsel_m{dim} = 2-mod(prevornext,2):size(wdx_dim2,dim)-mod(prevornext,2);
                        end;
                        wsel_n{1} = J;
                        wsel_m{1} = K;
                        wdx_n = wdx(wsel_n{:});
                        wdx_m = wdx_dim2(wsel_m{:});
                    end;
                        
%                     dim,dim2
%                     indxi{dim}
%                     indxn{dim}
%                     indxm{dim}
                    if ioffsetindim==0
                        wdx_n = -wdx_n;
                    end;                    
                    if ioffsetindim2
                        wdx_m = -wdx_m;
                    end;
                    indxn{1} = J;
                    indxm{1} = K;
                    
                    Hi{hcrow,dim} = linindx(indxn{:});
                    Hj{hcrow,dim} = linindx(indxm{:});
                    Hv{hcrow,dim} = bsxfun(@times, rf3(indxi{:}),  wdx_n .* wdx_m) ; %+ 1/f_i^3 * (wj * (xj_i - xj_n) ) * (wk * (xk_i - xk_m) )
                    % DEBUG:
%                     col = Hi{hcrow,dim}(1,:)==Hj{hcrow,dim}(1,:) & Hi{hcrow,dim}(1,:)==23;
%                     if any(col)
%                         [dim dim2 prevornext]
%                         [Hi{hcrow,dim}(:,col) Hj{hcrow,dim}(:,col) Hv{hcrow,dim}(:,col)]
%                     end;
                    hcrow = hcrow +1;
                    if dim~=dim2 || ioffsetindim~=ioffsetindim2 %~isequal(indxn{dim},indxm{dim})
                        Hi{hcrow,dim} = Hj{hcrow-1,dim};
                        Hj{hcrow,dim} = Hi{hcrow-1,dim};
                        Hv{hcrow,dim} = Hv{hcrow-1,dim};
                        hcrow = hcrow+1;
                    end;
                end;
            end;
        end;
        
        if hasmask
            linindxm = linindx(:,mask);
            rfm = rf(:,mask);
            rf3m = rf3(:,mask);
            neighborcntm = neighborcnt(:,mask);
            sumwdxm = sumwdx(:,mask);
        else
            linindxm = linindx;
            rfm = rf;
            rf3m = rf3;
            neighborcntm = neighborcnt;
            sumwdxm = sumwdx;
        end;
        % d2 f_i /dxj_i dxk_i =  
        if ismatweight
            Hi{1,1} = linindxm(J,:);
            Hj{1,1} = linindxm(K,:);
        else
            Hi{1,1} = linindxm;
            Hj{1,1} = linindxm;
        end;
        Hv{1,1} = bsxfun(@times, rfm .* neighborcntm , weights(:)); % 1/f_i   * (wj * sum_n   delta_kj    )
        Hv{1,1} = Hv{1,1}(:,:);
        Hi{2,1} = linindxm(J,:);
        Hj{2,1} = linindxm(K,:);
        Hv{2,1} = bsxfun(@times, -rf3m(:,:) , sumwdxm(J,:).*sumwdxm(K,:)); % - 1/f_i^3 * (wj * sum_n (xj_i - xj_n) ) * (wk * sum_n (xk_i - xk_n) ) 
        
        
        for k=1:numel(Hi)
            % DEBUG test 
%             if ~isequal(size(Hi{k}),size(Hj{k})) ||~isequal(size(Hi{k}),size(Hv{k}))
            if ~isequal(numel(Hi{k}),numel(Hj{k})) ||~isequal(numel(Hi{k}),numel(Hv{k}))
                error('incorrect size of Hi, Hj and/or Hv');
            end;
            Hi{k} = reshape(Hi{k},[],1);
            Hj{k} = reshape(Hj{k},[],1);
            Hv{k} = reshape(Hv{k},[],1);
        end;
        H = sparse(vertcat(Hi{:}),vertcat(Hj{:}),vertcat(Hv{:}),numel(x),numel(x));
        if 0
            q= [vertcat(Hi{:}),vertcat(Hj{:}),vertcat(Hv{:})];
            qind = find(q(:,1)<=3 & q(:,2)==1);
            q(qind,:)
        end;
        if hasmask
            % integrated mask in hessian computation, so we avoid computation of the full hessian.
            % But used original indices, so select the part within the mask here:
            rmask = repmat(mask(:)',[size(x,1), 1]);
            %DEBUG: check if any part is unnececarily filled:
%             if any(any(H(:,~rmask)))
%                 error('For efficiency should not fill in part that is removed.');
%             end;
            H = H(rmask,rmask);
        end;
    end;
end;    

%%
function [Hy] = TVhessmul_mex(hessinf, y)

if size(y,1)==numel(hessinf.img) && size(y,2)>1
    % multiple y
    Hy = cell(1,size(y,2));
    if hessinf.derivdir==1
        if hessinf.cyclic
            for k=1:numel(Hy);
                ndim = ndims(hessinf.img);
                szimg = size(hessinf.img);
                yk = reshape( y(:,k), szimg );
                sel = cell(1,ndim);
                sel{1} = ':';
                for dim = 2:ndim;
                    sel{dim} = [szimg(dim) 1:szimg(dim) 1];
                end;
                Hyk = TotalVariation_Vector_c(hessinf.img(sel{:}) ,hessinf.rspacing2,hessinf.weights,hessinf.offset, hessinf.H, yk(sel{:}));
                sela = repmat( {':'},1,ndim);
                for dim = 2:numel(sel);
                    selr = sela; selw= sela;
                    selr{dim} = [1 szimg(dim)+2];
                    selw{dim} = [szimg(dim)+1 2];
                    Hyk(selw{:}) = Hyk(selw{:}) + Hyk(selr{:});
                    selc = sela;
                    selc{dim} = [2:szimg(dim)+1];
                    Hyk = Hyk(selc{:});
                end;
                Hy{k} = Hyk(:);
            end;
        else
            for k=1:numel(Hy);
                Hy{k} = TotalVariation_Vector_c(hessinf.img,hessinf.rspacing2,hessinf.weights,hessinf.offset, hessinf.H, y(:,k));
            end;
        end;
    elseif hessinf.derivdir==2 || hessinf.derivdir==3 % central differences.
        if hessinf.derivdir==2
            gradvec = [-.5 0 .5];gradvecoffset = -1;
        else
            gradvec = [-.125 .75 0 -.75 .125];gradvecoffset = -2;
        end;
        % analytic hessian not implmented yet, computing numerical hessian
        % from gradient:
        for k=1:numel(Hy);
            yt = y(:,k);
            nrmsc = .0001/sqrt(sum(yt.^2));
            img_f = hessinf.img + nrmsc * reshape(yt,size(hessinf.img));
            [ff, gf] = TotalVariation_Vector_gradfilt_fc(img_f,hessinf.rspacing2,hessinf.weights,hessinf.offset, gradvec, gradvecoffset);
            img_b = hessinf.img - nrmsc * reshape(yt,size(hessinf.img));
            [fb, gb] = TotalVariation_Vector_gradfilt_fc(img_b,hessinf.rspacing2,hessinf.weights,hessinf.offset, gradvec, gradvecoffset);
            Hy{k} = (gf(:)-gb(:))/(2*nrmsc);
        end;
    else
        for k=1:numel(Hy);
            Hy{k} = TotalVariation_Vector_fc(hessinf.img,hessinf.rspacing2,hessinf.weights,hessinf.offset, hessinf.H, y(:,k));
        end;
    end;
    Hy = [Hy{:}]; % horizontal concatenation.
else
    % single y:
    if hessinf.derivdir==1 % forward and backward gradients under sqrt
         if hessinf.cyclic
            ndim = ndims(hessinf.img);
            szimg = size(hessinf.img);
            yk = reshape( y, szimg );
            sel = cell(1,ndim);
            sel{1} = ':';
            for dim = 2:ndim;
                sel{dim} = [szimg(dim) 1:szimg(dim) 1];
            end;
            Hyk = TotalVariation_Vector_c(hessinf.img(sel{:}) ,hessinf.rspacing2,hessinf.weights,hessinf.offset, hessinf.H, yk(sel{:}));
            sela = repmat( {':'},1,ndim);
            for dim = 2:numel(sel);
                selr = sela; selw= sela;
                selr{dim} = [1 szimg(dim)+2];
                selw{dim} = [szimg(dim)+1 2];
                Hyk(selw{:}) = Hyk(selw{:}) + Hyk(selr{:});
                selc = sela;
                selc{dim} = [2:szimg(dim)+1];
                Hyk = Hyk(selc{:});
            end;
            Hy = Hyk(:);
         else
            [Hy] = TotalVariation_Vector_c(          hessinf.img,hessinf.rspacing2,hessinf.weights,hessinf.offset,                         hessinf.H, y);
         end;
    elseif hessinf.derivdir==2 || hessinf.derivdir==3 % central differences.
        if hessinf.derivdir==2
            gradvec = [-.5 0 .5];gradvecoffset = -1;
        else
            gradvec = [-.125 .75 0 -.75 .125];gradvecoffset = -2;
        end;
        % analytic hessian not implmented yet, computing numerical hessian
        % from gradient:
        
        yt = y;
        nrmsc = .0001/sum(yt.^2);
        img_f = hessinf.img + nrmsc * reshape(yt,size(hessinf.img));
        [ff, gf] = TotalVariation_Vector_gradfilt_fc(img_f,hessinf.rspacing2,hessinf.weights,hessinf.offset, gradvec, gradvecoffset);
        img_b = hessinf.img - nrmsc * reshape(yt,size(hessinf.img));
        [fb, gb] = TotalVariation_Vector_gradfilt_fc(img_b,hessinf.rspacing2,hessinf.weights,hessinf.offset, gradvec, gradvecoffset);
        Hy = reshape( (gf-gb)/(2*nrmsc), size(yt));

    else % forward dderivatives only.
        [Hy] = TotalVariation_Vector_fc(         hessinf.img,hessinf.rspacing2,hessinf.weights,hessinf.offset,                         hessinf.H, y);
    end;
    Hy = Hy(:); % reshape into vector.
end;



function [Hy] = TVhessmul(hessinf, y)
% Multiplies vector image y with hessian of total variation
% This code is basically a copy of the explicit hessian create code
% with modification to not explicitly create hessian, but multiply immediately with x 

% sum_k  d2 f_i /dxl_i dxk_i yk_i  =  1/f_i   *                                                ( sum_n       (sum_k wkl * yk_i)                /sp_n^2  ) 
%                                   - 1/f_i^3 * ( sum_n sum_j wjl * ((xj_i - xj_n)/sp_n^2) ) * ( sum_n sum_j (sum_k wjk * yk_i) *((xj_i - xj_n)/sp_n^2) )
% sum_kn d2 f_i /dxl_i dxk_n yk_n  =- 1/f_i   *                                                ( sum_n       (sum_k wkl * yk_n)                /sp_n^2  )
%                                   + 1/f_i^3 * ( sum_m sum_j wjl * ((xj_i - xj_m)/sp_m^2) ) * ( sum_n sum_j (sum_k wjk * yk_n) *((xj_i - xj_n)/sp_n^2) )
% sum_km d2 f_i /dxl_n dxk_m yk_m  =  1/f_i   *                                                (             (sum_k wkl * yk_n)                /sp_n^2) ) 
%                                   + 1/f_i^3 * (       sum_j wjl *  (xj_i - xj_n)/sp_n^2  ) * ( sum_m sum_j (sum_k wjk * yk_m) * (xj_i - xj_m)/sp_m^2  )

% Hyl_i =   sum_k  d2 f_i /dxl_i dxk_i yk_i        + sum_kn d2 f_i /dxl_i dxk_n yk_n       + sum_kn d2 f_n /dxl_i dxk_n yk_n 
%         + sum_n sum_km d2 f_n /dxl_i dxk_m yk_m  
%       =       1/f_i   * (  ( sum_n       (sum_k wkl * yk_i)                /sp_n^2  ) - ( sum_n       (sum_k wkl * yk_n)                /sp_n^2  ) )
%         sum_n 1/f_n   *    (             (sum_k wkl * yk_i)                /sp_n^2) ) 
%         sum_n-1/f_n   *    (             (sum_k wkl * yk_n)                /sp_n^2  )
%             - 1/f_i^3 * ( sum_n sum_j wjl * ((xj_i - xj_n)/sp_n^2) ) * ( sum_n sum_j (sum_k wjk * yk_i) *((xj_i - xj_n)/sp_n^2) )
%             + 1/f_i^3 * ( sum_m sum_j wjl * ((xj_i - xj_m)/sp_m^2) ) * ( sum_n sum_j (sum_k wjk * yk_n) *((xj_i - xj_n)/sp_n^2) )
%       + sum_n 1/f_n^3 * ( sum_m sum_j wjl * ((xj_i - xj_m)/sp_m^2) ) * (       sum_j (sum_k wjk * yk_n) *((xj_i - xj_n)/sp_n^2) )
%       + sum_n 1/f_n^3 * (       sum_j wjl *  (xj_n - xj_i)/sp_n^2  ) * ( sum_m sum_j (sum_k wjk * yk_m) * (xj_n - xj_m)/sp_m^2  )  
%
%       = (sum_n ( 1/f_i + 1/f_n ) /sp_n^2 ) * ((sum_k wkl * yk_i)             ) 
%         -sum_n ( 1/f_i + 1/f_n ) /sp_n^2   * ((sum_k wkl * yk_n)             ) 
%       + sum_n 1/f_i^3 * ( sum_m sum_j wjl * ((xj_i - xj_m)/sp_m^2) ) * (      sum_j ((sum_k wjk * yk_n) - (sum_k wjk * yk_i)) *((xj_i - xj_n)/sp_n^2) )
%       + sum_n 1/f_n^3 * ( sum_m sum_j wjl * ((xj_i - xj_m)/sp_m^2) ) * (      sum_j  (sum_k wjk * yk_n)                       *((xj_i - xj_n)/sp_n^2) )
%       - sum_n 1/f_n^3 * (       sum_j wjl *  (xj_i - xj_n)/sp_n^2  ) * (sum_m sum_j  (sum_k wjk * yk_m)                       * (xj_n - xj_m)/sp_m^2  )  

if prod(hessinf.sizex)<numel(y)
    % handle multi vector case:
    Hy = cell(1,size(y,2));
    for k=1:size(y,2)
        Hy{k} = TVhessmul(y(:,k), hessinf);
    end;
    Hy = [Hy{:}];
    return;
end;

szx     = hessinf.sizex;
weights = hessinf.weights;
rfiprfn = hessinf.rfiprfn;
rf      = hessinf.rf;
rf3     = hessinf.rf3;
wdxc    = hessinf.wdxc;
spacing = hessinf.spacing;
sumwdx  = hessinf.sumwdx;
ndim = numel(szx);

rf3sumw = bsxfun(@times, rf3, sumwdx);
origszy = size(y);
y = reshape(y,szx); % also checks size
ismatweight = size(weights,1)==size(weights,2)&& size(weights,1)>1;
if ismatweight
    wy = reshape(weights * y(:,:), szx);
else
    wy = bsxfun(@times, weights(:), y);
end;
% Compute
% sum_n ( 1/f_i + 1/f_n ) * (      (sum_k wkl * yk_i)                /sp_n^2  ) 
% with rfiprfn = sum_n ( 1/f_i + 1/f_n ) /sp_n^2
Hy = bsxfun(@times, rfiprfn , wy );
T1 = zeros([1 szx(2:end)]);
T3 = zeros(szx);
for dim = 2:ndim
    selhi = repmat({':'},1, ndim);
    sello = selhi;
    selhi{dim} = 2:szx(dim);
    sello{dim} = 1:szx(dim)-1;
    zero = zeros([szx(1:dim-1) 1 szx(dim+1:end)]);
    % Compute
    % - sum_n  1/f_i * (     (sum_k wkl * yk_n)                /sp_n^2  ) 
    wykn = (cat(dim, wy(selhi{:}), zero) + cat(dim, zero, wy(sello{:})));
    rfkn = (cat(dim,  bsxfun(@times, rf(selhi{:})/spacing(dim)^2 , wy(selhi{:})), zero) + cat(dim, zero, bsxfun(@times, rf(sello{:})/spacing(dim)^2, wy(sello{:}) )));
    Hy = Hy - bsxfun(@times, rf/spacing(dim)^2, wykn)-rfkn;
    dwy = diff(wy,1,dim);
    % Compute 
    % (sum_n sum_j ((sum_k wjk * yk_n) - (sum_k wjk * yk_i)) *((xj_i - xj_n)/sp_n^2) )
    st = sum( dwy .* wdxc{dim},1);
%     T1(sello{:}) = T1(sello{:}) - st;
%     T1(selhi{:}) = T1(selhi{:}) - st;
    T3(sello{:}) = T3(sello{:}) + bsxfun(@times, rf3sumw(sello{:}),  -sum( y(selhi{:}) .* wdxc{dim} ,1));
    T3(selhi{:}) = T3(selhi{:}) + bsxfun(@times, wdxc{dim},  -sum( y(sello{:}) .* rf3sumw(sello{:}) ,1));

    T3(selhi{:}) = T3(selhi{:}) + bsxfun(@times, rf3sumw(selhi{:}),   sum( y(sello{:}) .* wdxc{dim} ,1));
    T3(sello{:}) = T3(sello{:}) + bsxfun(@times, wdxc{dim},   sum( y(selhi{:}) .* rf3sumw(selhi{:}) ,1));
    
    % with wdxc{dim} = ((-xj_i + xj_n)/sp_n^2        with n =i+1 in dim 
    
    for dim2 = 2:ndim
        for prevornext=1:4
            n_up = mod(prevornext,2)==1;
            m_up = prevornext<=2;
            if n_up
                seln = selhi;
                seli = sello;
                wdxseln = -wdxc{dim};
            else
                seln = sello;
                seli = selhi;
                wdxseln = wdxc{dim};
            end;
            if dim2~=dim
                seli{dim2} = 1:szx(dim2);
                seln{dim2} = 1:szx(dim2);
            end;
            if dim2==dim && n_up~=m_up
                selm = seli;
                if m_up
                    wdxselm = -wdxc{dim2};
                else
                    wdxselm =  wdxc{dim2};
                end;
            else
                selm = seln;
                if m_up
                    selm{dim2} = selm{dim2}(2:end);
                    seln{dim2} = seln{dim2}(1:end-1);
                    seli{dim2} = seli{dim2}(1:end-1);
                    wdxselm = -wdxc{dim2}(seln{:});
                else
                    selm{dim2} = selm{dim2}(1:end-1);
                    seln{dim2} = seln{dim2}(2:end);
                    seli{dim2} = seli{dim2}(2:end);
                    wdxselm =  wdxc{dim2}(selm{:});
                end;
                if n_up
                    wdxseln = wdxseln(seli{:});
                else
                    wdxseln = wdxseln(seln{:});
                end;                    
            end;
            Hy(seli{:}) = Hy(seli{:}) + bsxfun( @times, wdxseln, rf3(seln{:}) .* sum( y(selm{:}) .* wdxselm, 1 ) );
        end;
    end;
end;
% Hy = Hy + bsxfun(@times, sumwdx, T1.*rf3) + T3;
Hy = Hy - bsxfun(@times, sumwdx, rf3.*sum(sumwdx.* y,1)) + T3;
Hy = reshape(Hy,origszy);
%%
function valididatecode
%%
if 1 
    sz = [3 4 5];
else
    sz = [3 7];
end;
tstv = rand(prod(sz),1);sp = rand(1,numel(sz));
if 0 
    w = rand(sz(1),1);
else
    w = rand(sz(1),sz(1)+2); w = w * w';
end;
if 1 
    tstmask = [];
else
    % test with mask does not yet work:
    tstmask = true(sz);
    tstmask(:,:,1)=false;
end;
    
% % define functions:
fun     = @(q) totalVariationVecRegularizer(reshape(q,sz) , sp, w, [], tstmask);
funHmul = @(q) totalVariationVecRegularizer(reshape(q,sz) , sp, w, [], tstmask,1);
Hmulfun = totalVariationVecRegularizer( [] , sp, w, [], tstmask,2);

% % compute function value, gradient and hessian + compare with numerical result:
Hest = jacobianest( @(q) shuffleoutputs(fun, 2 , 2,{q}), tstv);
% %
figure(1)
[f, g, H] = fun(tstv);H = full(H);
% %
[ff, gg, Hifo] = funHmul(tstv);HH = Hmulfun( Hifo, eye(numel(gg)));
% %
imagebrowse(cat(3,full(H)-Hest, full(H)-HH),[-1 1])
set(gcf,'name','Difference between numerical and analytic hessian.')
disp([' max hessian difference: '  num2str(max(abs(HH(:)-H(:)))) ', maximum error : ' num2str(max(abs(Hest(:)-H(:))))])
figure(2);
gest = gradest(fun, tstv);
plot([g(:)-gest(:)])
title('Difference between numerical and analytic gradient.')
%%
sz = [1 3 4];
tstv = randn(sz);
sp  = [1 1 1];
w = 1;
tstmask = [];
fun     = @(q) totalVariationVecRegularizer(reshape(q,sz) , sp, w, [], tstmask);

%%
[f,g,h] = fun( tstv );
xt = linspace(-3,3,100)';
adjidx = ([1 2 2]-1)*[1;sz(1);sz(1)*sz(2)]+1;
ftrace = zeros(numel(xt),1);
for k=1:numel(xt)
    tstv2 = tstv;
    tstv2(adjidx) = xt(k);
    ftrace(k) = fun( tstv2 );
end;
%%
plot(xt, [ftrace f+g(adjidx)*(xt-tstv(adjidx))+h(adjidx,adjidx)*(xt-tstv(adjidx)).^2/2]);
function test_cyclic
r=randn(3,7,8);sp = 1+rand(1,ndims(r));weights = rand(size(r,1),1);offset =.01; mask =[]; outputhessmul = 1; derivdir = 1; cyclic = true
TVhessmulfun = totalVariationVecRegularizer(r , sp , weights, offset , mask, 2, derivdir, cyclic);
TVfunc = @(x) totalVariationVecRegularizer(x , sp , weights, offset , mask, outputhessmul, derivdir, cyclic);

[f, g, H] = TVfunc(r) ;

startidx = round(rand(1,ndims(r)).*size(r));
rot = cell(1,ndims(r));
for dim=1:ndims(r)
    rot{dim} = mod(startidx(dim)+(1:size(r,dim)),size(r,dim))+1;
end;
rot{1} = ':';
[fr, gr, Hr] = TVfunc( r(rot{:}) );
f-fr % should be zero as function should now be cyclic invariant.
%% test gradient and hessian:
dummy = validateDerivativeAndHessian( TVfunc, r, TVhessmulfun);
%%
imagebrowse(reshape( cat(3,dummy.H,dummy.Hest), [size(r) prod(size(r)) 2] ),[],'colordim',[])