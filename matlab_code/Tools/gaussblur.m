function [out, filt] = gaussblur( img, sigmas, edges)
% [out] = gaussblur( img, sigmas [, edges])
% Performs a accurate and quick gaussian blurr of the N dimensional image img.
% 
% INPUTS
% img  : arbitrary ND image (single or double precision values)
% sigmas: scalar or vector that specifies the smoothing in each dimension. 
%         A sigma value of 0 indicates absolutely no smoothing in that dimension.
% edges : optional. Default = 'correct'
%         If 'correct' then the edge intensity is corrected for the reduced number of samples
%         that these have. 
%         Otherwise it is just a gaussian blur with (implicitly) extending the 
%         image with zeros in all dimensions.
%
% created by Dirk Poot (d.poot@erasmusmc.nl), Erasmus Medical center, Rotterdam
% 31-8-2011

if nargin<3
    edges = 'correct';
end;
docorrect = isequal(edges,'correct');
ndim = ndims(img);
if numel(sigmas)==1
    sigmas = sigmas(ones(1,ndim));
end;

if ~isreal(img)
    out = complex( gaussblur( real(img), sigmas, edges), gaussblur( imag(img), sigmas, edges));
    return;
end;

filt = cell(1,ndim);
for dim=1:ndim
    if sigmas(dim)==0
        filt{dim} = [];
    else
        tmp = exp(-(0:max(1,floor(3.5*sigmas(dim)))).^2/(2*sigmas(dim).^2));
        filt{dim} = [tmp(end:-1:2) tmp];
        filt{dim} = filt{dim} / sum(filt{dim});
    end;
end;

out = separableConvN_c( img , filt );    
if docorrect
    % compensate edges:
    for dim=1:ndim
        comp_dim = 1./separableConvN_c( ones(size(img,dim),1) , filt(dim) ); 
        comp_dim = reshape(comp_dim,[ones(1,dim-1) size(img,dim) 1]);
        out = bsxfun(@times, out, comp_dim);
    end;
end;
