function [data] = fftresample(data, newsz, dim, compensate)
% [dataout] = fftresample(data, newsize [, dim, [compensate]])
% [dataout] = fftresample(data, newsize [, dim, [shift]])
% Resamples data with the fourier transform. The data is treated as
% periodic signal and by using the fft an optimal sinc
% interpolation/resampling of this periodic data matrix is computed. 
% This function resamples, so it stretches/compresses the data in the
% dimensions where the size changes. 
% INPUTS:
% Data    : n dimensional image.
% newsize : n element vector specifiing the new size of the data matrix, or
%           scalar when combined with dim.
% dim     : the dimension in which the size of data is changed (when
%           newsize is a scalar)
% compensate : 'linear' : pre-compensate for subsequent linear
%                         interpolation. 
%                         -  Only usefull when inceasing the size.
%                         -  Apllied only on dimensions which change size.
%               a peacewise polynomial (pp) structure (from spline) 
%                       : compensation factor is ppval(pp,freq) where freq is the
%                         frequency in the resampled domain. 1 is the output sample
%                         frequency, so pp should be defined from 0 to 0.5 (= nyquist
%                         frequency).
%               a function handle to a function with 1 argument (the frequency), 
%               that evaluates the required compensation factor for each frequency.
% shift       : size(newsize) array that specifies shifts (implicitly) applied to the 
%               source image. The output is the same (except for roundoff errros) as 
%               when calling fftresample( fracshift(data, shift), newsz, dim)
%
% OUTPUTS:
% dataout : resampled version of data.
%
% Created by Dirk Poot, University of Antwerp.

% 25-2-2009 : original version
% 5-3-2009 : efficiently increase/decrease size in multiple dimensions.

if nargin<4
    compensate = [];
end;
if nargin>=3 && ~isempty(dim)
    ds = newsz;
    newsz = size(data);
    newsz(dim) = ds;
    if ~isempty(compensate) && ~ischar(compensate) && ~isstruct(compensate) && ~isa(compensate,'function_handle') && numel(compensate)==1
        tmp = compensate;
        compensate = zeros(size(newsz));
        compensate(dim) = tmp;
    end;
end;
datasz = size(data);
if numel(datasz)<numel(newsz)
    datasz(end+1:numel(newsz)) =1;
elseif numel(datasz)>numel(newsz)
    newsz(end+1:numel(datasz)) = datasz(numel(newsz)+1:end);
end;
isrealdata = isreal(data);
[ratio,order] = sort(newsz./datasz); % optimize order to maximize speed and minimize memory usage.
order(ratio==1)=[];

for dim = order%find( newsz ~= datasz )
    dataf = fft(data,[],dim);
    sel = cell(numel(datasz),1);
    sel(:)={':'};
    
    if ~isempty(compensate)
        f= min(0:datasz(dim)-1, datasz(dim):-1:1)/newsz(dim); % frequency between 0 and 1 (= sampling frequency; is first frequency aliased to 0)
        if ischar(compensate)
            switch lower(compensate)
                case 'linear'
                    adj = 1./(sinc(f).^2);
                otherwise
                    error('wrong compensation requested');
            end;
        elseif isstruct(compensate)
            % assume pp structure
            adj = ppval(compensate, f);
        elseif isa(compensate,'function_handle')
            adj = compensate( f );
        elseif isequal(size(compensate),size( newsz))
            sc = -[(0:floor(datasz(dim)/2)-1) (-ceil(datasz(dim)/2):-1)]'/datasz(dim)*2*pi*compensate(dim);
            adj = exp(1i * sc );
        else
            error('wrong compensation requested');
        end;
            
        adjsz = ones(size(datasz));
        adjsz(dim) = datasz(dim);
        dataf = bsxfun(@times, dataf , reshape(adj,adjsz));
    end;
    
    if newsz(dim)>datasz(dim)
        % increase size
        sel2 = sel;
        zerosz = datasz;
        halfsz = datasz(dim)/2 ;
        if halfsz ==ceil(halfsz )
            % even size; split center line
            sel{dim} = (halfsz+1);
            center = dataf(sel{:})/2; % select center line.

            sel{dim} = (1:halfsz);
            sel2{dim} = (halfsz+2:datasz(dim));
            zerosz(dim) = newsz(dim)-datasz(dim)-1;
        else
            % odd size; no center line.
            sel{dim} = (1:ceil(halfsz));
            sel2{dim} = (sel{dim}(end)+1:datasz(dim));
            zerosz(dim) = newsz(dim)-datasz(dim);
            center = [];
        end;
        dataf = cat(dim,dataf(sel{:}),center,zeros(zerosz),center,dataf(sel2{:}));
    else
        % reduce size:
        halfsz = newsz(dim)/2 ;
        if halfsz ==ceil(halfsz )
            % even new size, add 2 lines to create center line.
            sel2 = sel;
            sel{dim} = (halfsz+1);
            sel2{dim} = datasz(dim)-halfsz+1;
            center = dataf(sel{:})+dataf(sel2{:});
            
            sel{dim} = (1:halfsz);
            sel2{dim} = datasz(dim)+(-halfsz+2:0);
            
            dataf = cat(dim,dataf(sel{:}),center,dataf(sel2{:}));
        else
            % odd new size, just join.
            sel{dim} = [(1:ceil(halfsz)) datasz(dim)+(ceil(-halfsz)+1:0)];
            dataf = dataf(sel{:});
        end;
    end;
    data = ifft(dataf,[],dim)*(newsz(dim)/datasz(dim));
    datasz(dim) = newsz(dim);
end;
if isrealdata
    data = real(data);
end;

return;
%%
if 0
    % test cases:
    clf; hold all
    q=randn(6);q2=fftresample(q, [12 12]);q3=fftresample(q2, [6 6]);q-q3
    plot((0:size(q,1)-1)/size(q,1),q(:,1),(0:size(q2,1)-1)/size(q2,1),q2(:,1))
    q=2+randn(6);q2=fftresample(q, [11 11]);q3=fftresample(q2, [6 6]);q-q3
    plot((0:size(q,1)-1)/size(q,1),q(:,1),(0:size(q2,1)-1)/size(q2,1),q2(:,1))
    q=4+randn(5);q2=fftresample(q, [11 11]);q3=fftresample(q2, [5 5]);q-q3
    plot((0:size(q,1)-1)/size(q,1),q(:,1),(0:size(q2,1)-1)/size(q2,1),q2(:,1))
    q=6+randn(5);q2=fftresample(q, [10 10]);q3=fftresample(q2, [5 5]);q-q3
    plot((0:size(q,1)-1)/size(q,1),q(:,1),(0:size(q2,1)-1)/size(q2,1),q2(:,1))
    hold off
end;