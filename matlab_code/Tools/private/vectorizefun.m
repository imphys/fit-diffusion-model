function [f, g, h] = vectorizefun( fun, x , fields )
% [f, g, h] = vectorizefun( fun, x [, fields] )
% vectorizes function 'fun' that only supports a single column vector parameter as input.
% Support function for fit_MRI
% 
%
% Created by Dirk Poot, Erasmus MC, 2-2-2012
args = cell(nargout, size(x,2));
if nargin==2
    % No fields case:
    for k=1:size(x,2);
        [args{:,k}] = fun( x(:,k) );
    end;
else%if nargin==3
    % fields also provided.
    for k=1:size(x,2);
        [args{:,k}] = fun( x(:,k), fields(:,k) );
    end;
end;
f = reshape(horzcat(args{1,:}), [],size(x,2));
if nargout>=2
    if ~isequal(size(squeeze(args{2,1})), [size(f,1) size(x,1)])
        error('incorrect gradient size, expected (nMRI x npar) or (nMRI x 1 x npar)');
    end;
    for k=1:size(x,2)
        args{2,k} = reshape(args{2,k},[size(f,1),1,size(x,1)]);
    end;
    g = cat(2,args{2,:});
    if nargout>=3
        nh = size(x,1)*(size(x,1)+1)/2; % number of hessian elements
        if ~isequal(size(squeeze(args{2,1})), [size(f,1) nh])
            error('incorrect hessian size, expected (nMRI x npar*(npar+1)/2) or (nMRI x 1 x npar*(npar+1)/2)');
        end;
        for k=1:size(x,2)
            args{3,k} = reshape(args{3,k},[size(f,1), 1, nh]);
        end;
        h = cat(2,args{3,:});
    end;
end;