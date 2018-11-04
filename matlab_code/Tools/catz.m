function [res]= catz(dim, varargin)
% function [res]= catz(dim, varargin);
% Concatenates arrays like CAT, but expands elements that are
% smaller with zero.
% 
% Created by Dirk Poot, University of Antwerp.
% 5-2-2007 : original version

numars = nargin-1;
classname = [];
allagree = true;
for k=1:numel(varargin)
    if isempty(classname)
        classname = class(varargin{k});
    else
        allagree = allagree && isequal(classname,class(varargin{k}));
    end;
end;
if ~allagree  
    classname = superiorfloat(varargin{:});
end;
sz = zeros(max(2,dim),numars);
for k=1:numars
    szk = size(varargin{k});
    sz(1:numel(szk),k) = szk(:)-1; % make sure all elements are 1 by default.
end;
sz = sz + 1; % make sure all elements are 1 by default.
maxsz = max(sz,[],2);
catsz = cumsum([1 sz(dim,:)]);
ressz = maxsz;
ressz(dim) = catsz(end)-1;
res = zeros(ressz',classname);
asgnc = cell(size(sz,1),1);
for k=1:numars
    for l=1:size(sz,1)
        asgnc{l} = 1:sz(l,k);
    end;
    asgnc{dim} = catsz(k):catsz(k+1)-1;
    res(asgnc{:}) = varargin{k};
end;