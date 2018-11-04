function [curcoord, fval] = findlocalminima(mat, neigboursonly, cyclic, interpolate)
% [coords, val] = findlocalminima(mat , neigboursonly , cyclic, interpolate)
% computes the positions that are strict local minima of the n-dimensional
% matrix mat. 
% INPUTS:
% mat : N-dimensional matrix of which the minima are located
% neighboursonly : Default = 0; 
%       Specify 1 to only compare to the direct neigbours (4 in 2-D
%       matrices), by default all positions with an infinity norm<=1 from
%       each point are used (so all points with no more than 1 step in all
%       directions)  
% cyclic         : Default = 0;
%       Specify 1 or ones(size(size(mat))) to treat the matrix mat as
%       cyclic in all dimensions, or specify 1/0 for each dimension separately.
%       Specifying 1 (for the first dimension) means e.g. that mat(1,k) is a local
%       minimum only when mat(1,k)<mat(end,k) (in all dimensions)
% interpolate    : Default = 0;
%       Specify 1 to compute interpolated minima. A 2nd order polynomial (parabola)
%       is fitted to each local minimum (using the full 3^N neighbourhood), and the
%       minimum of that parabola is returned
%
% OUTPUTS:
% coords : matrix with in each row the location of a (interpolated) minimum of mat.
%          So mat(coords(k,1), ... coords(k,N)) is a local minimum. 
% val    : Column vector with the value of the matrix in the minima; might be
%          especially usefull when combined with the interpolate option. When the
%          output argument val is requested, the minima are sorted with respect to
%          val. When only 1 output argument is requested the results are unsorted!
%
% Created by Dirk Poot, University of Antwerp.

% original creation date unknown, added credits: 23-7-2008
% 24-7-2008: removed bug; used indloc instead of ind at last line.
% 7-8-2009 : added interpolate option.

if nargin<2 || isempty(neigboursonly)
    neigboursonly = false;
end;
if nargin<3 || isempty(cyclic)
    cyclic = 0;
end;
if nargin<4 || isempty(interpolate)
    interpolate = false;
end;
smat = size(mat);
if numel(smat)==2 && smat(2)==1
    smat = smat(1);
end;
if any(smat==1)
    error('All dimensions should have at least 2 elements (column vectors allowed)');
end;
if numel(cyclic)==1
    cyclic = cyclic*ones(1,numel(smat));
end;

ndims = numel(smat);
minim = true(size(mat));
indices = cell(ndims,1);
indices(:) = {':'};
for k=1:ndims
    der = diff(mat,1,k);
    curind = indices;
    curind(k) = {1:smat(k)-1};
    minim(curind{:}) = minim(curind{:}) & der>0;
    curind(k) = {2:smat(k)};
    minim(curind{:}) = minim(curind{:}) & der<=0;
    if cyclic(k)
        curinde = curind;
        curind(k) = {1};
        curinde(k) = {smat(k)};
        der = mat(curind{:})-mat(curinde{:});
        minim(curind{:}) = minim(curind{:}) & der <0;
        minim(curinde{:}) = minim(curinde{:}) & der >0;
    end;
    clear der
end;
curcoord = find(minim);
ind = cell(1,ndims);
[ind{:}] = ind2sub(smat,curcoord);
curcoord = horzcat(ind{:});
if ~neigboursonly && numel(smat)>1
    clear minim 
    % all direct neighbours are tested against, now do test all others;
    steps = zeros(1,0);
    for k=1:ndims
        steps = [kron([-1;0;1],ones(size(steps,1),1)) kron(ones(3,1),steps)];
    end;
    steps(sum(abs(steps),2)<=1,:)=[];
    cumsz = cumprod([1 smat(1:end-1)])';

    ind = curcoord-1;
    clear curcoord
    % process each step, remove indices that are not lowest.
    for m=1:size(steps,1);
        indloc = ind;
        indtst = ind + steps(m*ones(size(ind,1),1),:);
        for k_dim = find(cyclic)
            indtst(:,k_dim) = mod(indtst(:,k_dim),smat(k_dim));
        end;
        notdel = ~(any(indtst<0,2) | any(indtst>=smat(ones(size(ind,1),1),:),2));
        dosel = ~all(notdel);
        if dosel
            repr = find(notdel);
            indloc = indloc(repr,:);
            indtst = indtst(repr,:);
        end;
        notkeep = mat(indloc*cumsz+1)>=mat(indtst*cumsz+1);
        if any(notkeep)
            if dosel
                ind(repr(notkeep),:) =[];
            else
                ind(notkeep,:) =[];
            end;
        end;
    end;
    curcoord = ind+1;
end;
if interpolate
    deltax = cell(1,ndims);
    curcoord = curcoord';
    cumsz = cumprod([1 smat(1:end-1)]);
    ind = cumsz*curcoord - sum(cumsz(2:end));
    x = mat2cell(curcoord,ones(1,size(curcoord,1)),size(curcoord,2));
    selindxx =zeros(0,1);
    for k=1:ndims
        deltax{k} = [-1 0 1]'* ones(1,numel(x{k}));
        if any(x{k}==1)
            if cyclic(k)
                deltax{k}(1,x{k}==1) = smat(k)-1;
            else
                deltax{k}(:,x{k}==1) = deltax{k}(:,x{k}==1)+1;
                warning('findlocalminima:InterpolateOnBorder','speculative interpolation on border.');
%                 error('cannot yet interpolate on border');
            end;
        end;
        if any(x{k}(:)==smat(k))
            if cyclic(k)
                deltax{k}(3,x{k}==smat(k)) = -smat(k)+1;
            else
                deltax{k}(:,x{k}==smat(k)) = deltax{k}(:,x{k}==smat(k))-1;
                warning('findlocalminima:InterpolateOnBorder', 'speculative interpolation on border.');
%                 error('cannot yet interpolate on border');
            end;
        end;
        deltax{k} = deltax{k}*cumsz(k);
        selindxx = [repmat(selindxx,1,3);kron([1 2 3]+(k-1)*3,ones(1,size(selindxx,2)))];
    end;
    deltax = vertcat(deltax{:});
    deltax = reshape(sum(reshape(deltax(selindxx,:),[size(selindxx,1) size(selindxx,2) size(deltax,2)]),1),[size(selindxx,2) size(deltax,2)]);
    
    selmat = mat(ind(ones(1,size(deltax,1)),:)+deltax); % select neigborhood around all local minima.
    
    [indx,cnt] = symmetricTensorElements(2,ndims);
    % X * tht = value at all 3^ndims neigborhood points. We are looking for tht
    % 2D case: fun(x,y) = a + b x + c y + d x^2 + e x y + f y^2 with tht = [a b c d e f]
    % 
    stepx = selindxx'-3*ones(3^ndims,1)*(0:ndims-1)-2;
    X = [ones(3^ndims,1) stepx stepx(:,indx(:,1)).*stepx(:,indx(:,2))];
    tht = X\selmat;
    
    % now solve general 2nd order polynomial:
    grad = tht(2:ndims+1,:)/2;
    hesssel = diag(1./cnt')*tht(ndims+2:end,:);
    opt = zeros(ndims,size(tht,2));
    hess = zeros(ndims,ndims);
    for k=1:size(tht,2)
        hess(indx*[1;ndims]-ndims)=hesssel(:,k);
        hess(indx*[ndims;1]-ndims)=hesssel(:,k);
        opt(:,k) = hess \grad(:,k);
    end;
    if any(abs(opt(:)))>1
        % might be a warning instead?
        % We probably want to solve this by testing if this actually is another local
        % minimum we found. Not implemented yet.
        error('interpolated minimum outside of interpolation region.');
    end;
    curcoord = (curcoord-opt)';
    if nargout>1
        X = [ones(1,size(opt,2)); -opt; opt(indx(:,1),:).*opt(indx(:,2),:)];
        fval = sum(X.*tht,1)';
    end;
elseif nargout>1
    % not interpolate, fval requested
    cumsz = cumprod([1 smat(1:end-1)]);
    ind = curcoord*cumsz' - sum(cumsz(2:end));
    fval = mat(ind);
end;
if nargout>1
    % sort:
    [fval, ind] = sort(fval);
    curcoord = curcoord(ind,:);
end;