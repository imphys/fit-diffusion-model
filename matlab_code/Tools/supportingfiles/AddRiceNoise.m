function [Data] = AddRiceNoise(Data, sigm)
% [Data] = AddRiceNoise(Data, sigm)
% Adds rice distributed noise to the (noise free) Data.
% This is done by adding complex noise (with standarddeviation sigm) to the
% data and then computing the magnitude of the resulting signal. 
% If Data is 0 the resulting distribution is a Rayleigh distribution.
%
% Created by Dirk Poot, University of Antwerp.

% last modified: 9-4-2008, renamed from AddRayleighNoise.

mx = numel(Data);
step = 1000;
szd = size(Data);
Data = reshape(Data, [],1);
for k=1:step:mx
    ed = min(mx,k+step-1);
    n = ed-k+1;
    Data(k:ed) = abs(sigm*(randn(n,1)+1i*randn(n,1)) + Data(k:ed));
end;
Data = reshape(Data, szd);