function [image, L] = laplaceMulND( image, lambda, voxelspacing, cyclic)
% [ f , L ] = laplaceMulND( image, lambda [, voxelspacing [,cyclic] ])
% Function to multiply with a laplacian.
% Computes f such that 
%   f(:)' * image(:) = image(:)'* L' * L * image(:) = lambda * sumsum( laplacian( image ).^2 ) 
%                 =approx=  sumsum( ( sum_dim diff(image, 2, dim) / voxelspacing(dim) ).^2 ) % but with more attention to edges.
% 
% in which laplacian( image ) is computed with discrete second derivative, 
% taking the voxel spacing into account (when provided).
% Explicitly, for 2D images:
%   f(:)' * image(:) = sum(sum( (diff(image(:,2:end-1) ,2,1)/voxelspacing(1) + diff(image(2:end-1,:) ,2,2)/voxelspacing(2)).^2 ))
%                      + sum(sum( (diff(image(:,[1 end]) ,2,1)/voxelspacing(1)).^2 ))
%                      + sum(sum( (diff(image([1 end],:) ,2,2)/voxelspacing(2)).^2 ))
%
%
% INPUTS:
%  image  : N-dimensional image
%  lambda : scalar with which the result is scaled.
%  voxelspacing :  default = ones(1,ndims(image))
%                  Specifies the voxel spacing with which the laplacian is computed.
%  cyclic : default = false;
%               true  : input is cyclic (=wraps around) in each dimension
%               false : input not cyclic.
%
% OUTPUTS:
%  f      : reshape(L'*L * image(:), size(image))
%           the image multiplied by the square laplacian.
%           Note L is not computed explicitly.
%  L      : only if 2 output arguments are requested L is computed explicitly as sparse matrix. 
%           Note that L does contain a substantial number of non zeros, thus 
%           for large images not enough memory might be available.
% 
% Created by Dirk Poot
% Erasmus MC 23-3-2011

if nargin>=3 && ~isempty( voxelspacing )
    rvoxelspacing = {1./voxelspacing};
else
    rvoxelspacing = {};
end;
if nargin<4 
%     cyclic = false;
    % don't need to update rvoxelspacing as default is mex file is false.
else
    if isempty(rvoxelspacing)
        rvoxelspacing = {ones(1,ndims(image)), cyclic};
    else
        rvoxelspacing{2} = cyclic;
    end;
end;
% if cyclic;
%     % currently no optimization for cyclic version. L computed explicitly.
%     szim = size(image);
%     ndim = numel(szim);
%     if isempty(rvoxelspacing)
%         rvoxelspacing = ones(1,ndim);
%     else
%         rvoxelspacing = rvoxelspacing{1};
%     end;
%     stepdim = [0 kron(1:ndim,[-1 1])];
%     linindx = reshape(1:prod(szim),szim);
%     I = linindx(:)*ones(1,2*ndim+1);
%     J = I;
%     V = I;
%     mulf = [0 kron( reshape( rvoxelspacing(1:ndim) , 1, []),[ 1 1])];
%     mulf(1) = -sum(mulf);
%     sel = repmat({':'},1,ndim);
%     for dim = 1:numel(mulf);
%         seldim = sel;
%         workdim = abs(stepdim(dim));
%         stepdir = sign(stepdim(dim));
%         if workdim~=0
%             seldim{ workdim } = mod( (0:szim(workdim)-1)+stepdir , szim(workdim) ) +1;
%         end;
%         J(:,dim) = reshape( linindx( seldim{:}) , [] , 1);
%         V(:,dim) = mulf(dim);
% 
%     end;
%     L = sparse(I,J,V, prod(szim), prod(szim));
%     image = reshape( L'*(L*image(:)), szim);
%     return;
% else
    
newversion = true;
if newversion
    % new c++ routine.
    % applies only L, but this L is symmetrical:
    image = laplaceMulND_v2( image, lambda, rvoxelspacing{:});
    image = laplaceMulND_v2( image, 1     , rvoxelspacing{:});
    if nargout>=2 
        szim = size(image);
        ndim = numel(szim);
        if isempty(rvoxelspacing)
            rvoxelspacing = ones(1,ndim);
        else
            rvoxelspacing = rvoxelspacing{1};
        end;
        stepdim = [1 cumprod(szim(1:end-1))];
        linindx = reshape(1:prod(szim),szim);
        I = cell(ndim,1);
        J = I;
        V = I;
        sel = repmat({':'},1,ndim);
        for dim = 1:ndim;
            seldim = sel;
            seldim{dim} = 1:szim(dim)-1;
            basei = reshape( linindx(seldim{:}), 1 , []);
            I{dim} = bsxfun(@plus, [0;0;stepdim(dim);stepdim(dim)],basei);
            J{dim} = bsxfun(@plus, [stepdim(dim); 0; 0; stepdim(dim)], basei);
            V{dim} = ([1;-1;1;-1]*rvoxelspacing(dim)*sqrt(lambda)) * ones(size(basei));

            I{dim} = reshape(I{dim},[],1);
            J{dim} = reshape(J{dim},[],1);
            V{dim} = reshape(V{dim},[],1);
        end;
        I = vertcat(I{:});
        J = vertcat(J{:});
        V = vertcat(V{:});
        L = sparse(I,J,V, prod(szim), prod(szim));
    end;
    
    return;
end;
    
% Call c++ subroutine that does the hard work:
image = laplaceMulND_c(image, lambda, rvoxelspacing{:});

if nargout>=2 
    szim = size(image);
    ndim = numel(szim);
    if ndim>4
        warning('L is not computed exactly');
    end;
    if isempty(rvoxelspacing)
        rvoxelspacing = ones(1,ndim);
    else
        rvoxelspacing = rvoxelspacing{1};
    end;
    stepdim = [1 cumprod(szim(1:end-1))];
    linindx = reshape(1:prod(szim),szim);
    I = cell(ndim,1);
    J = I;
    V = I;
    sel = repmat({':'},1,ndim);
    for dim = 1:ndim;
        seldim = sel;
        seldim{dim} = 2:szim(dim)-1;
        I{dim} = linindx(seldim{:});
        I{dim} = reshape( I{dim} , 1, []);
        J{dim} = bsxfun(@plus, [-stepdim(dim); 0; stepdim(dim)], I{dim});
        I{dim} = [1;1;1]*I{dim};
        V{dim} = ([1;-2;1]*rvoxelspacing(dim)*sqrt(lambda)) * ones(1,size(I{dim},2));
        
        
        % fixup for low dimensional edges which 'laplaceMulND_c' handles kind of strangely (duplicates)
        % I keep it this way since adding extra weight to these edges probably is not a bad idea.
        remaindim = [1:dim-1 dim+1:ndim]; 
        for dim2 = remaindim(1:end-1)
            remaindim2 = remaindim(remaindim>dim2);
            for dim3= remaindim2
                seldim2 = sel;
                seldim2{dim2} = [1 szim(dim2)];
                seldim2{dim3} = [1 szim(dim3)];
                szimr = szim;szimr(dim)=max(0,szimr(dim)-2);
                V{dim} = reshape(V{dim},[3 szimr]);
                V{dim}(:,seldim2{:}) = sqrt(2)*V{dim}(:,seldim2{:});
                remaindim3 = remaindim2(remaindim2>dim3);
                for dim4 = remaindim3
                    seldim4 = seldim2;
                    seldim4{dim4} = [1 szim(dim4)];
                    V{dim}(:,seldim4{:}) = sqrt(3)/2*V{dim}(:,seldim4{:});
                end;
            end;
        end;
        I{dim} = reshape(I{dim},[],1);
        J{dim} = reshape(J{dim},[],1);
        V{dim} = reshape(V{dim},[],1);
    end;
    if ndim>=3
        for dim1=1:ndim
            for dim2 = dim1+1:ndim
                
            end;
        end;
    end;
    I = vertcat(I{:});
    J = vertcat(J{:});
    V = vertcat(V{:});
    L = sparse(I,J,V, prod(szim), prod(szim));
end;
% end;