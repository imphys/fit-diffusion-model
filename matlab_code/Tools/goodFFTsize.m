function sz = goodFFTsize(n)
% sz = goodFFTsize( n )
% Finds a good size for an FFT with at least n samples.
% n might scalar or an array. 
%
% Created by Dirk Poot, University of Antwerp.
% Original version: ?
% 2-4-2009 : updated help

persistent tbl;
if isempty(tbl)
    fastprimes = [2 3 5];
    maxtbl = 10470;
    
    lgfp = log(fastprimes);
    lgmt = log(maxtbl);
    cnt = zeros(numel(fastprimes),1);
    cnt(1) =1;
    lp = true;
    while lp
        tbl(end+1) = prod(fastprimes.^(cnt'));
        k=1;
        while k<=numel(cnt)
            cnt(k) = cnt(k)+1;
            if lgfp*cnt>lgmt 
                cnt(k)=0;
                k=k+1;
            else
                break;
            end;
        end;
        lp = k<=numel(cnt);
    end;
    tbl = [-inf sort(tbl)];
end;

[cnt, ind] = histc(n,tbl);
sind = ind(ind~=0);

incr = n(ind~=0)>tbl(sind);
sind(incr) = sind(incr)+1;
sz = zeros(size(n));

sz(ind~=0) = tbl(sind);
for k=find(ind==0)
    lgn = log(n(k));
    n3 = 0:ceil(lgn/log(3)); % test all powers of 3 that potentially might give the best value.
    n2 = max(0,ceil( (lgn - n3*log(3)) / log(2))); % find lowest power of 2 that gives a value larger than n(k)
    z = 2.^n2.*3.^n3; % compute the test values (all are larger than n(k) )
    sz(k) = min(z); % find lowest test value and return that.
end;