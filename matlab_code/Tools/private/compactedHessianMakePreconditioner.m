function [Rs, mulfun] = compactedHessianMakePreconditioner(hessinfo, upperbandw, DM, DG )
% [R,PVEC] = HPRECON(H,UPPERBANDW,DM,DG) computes the 
% prepare R as preconditioner of 
%    M = DM*H*DM + DG
% so 
%    compactedHessianMulPreconditioner(x, R)
% approximates inv(M) * x
%
% Created by Dirk Poot (d.poot@erasmusmc.nl), Erasmus medical center, Rotterdam
if nargin<1
    Rs = @compactedHessianMulPreconditioner2;
    return;
end;
if isstruct(hessinfo) 
    if isfield(hessinfo,'Hfun')
        % for regularizedHessMulFun
        hessinfo = hessinfo.Hfun;
    end;
    if isfield(hessinfo,'psfdiag')
        [I,J] = find(triu(ones(size(hessinfo.dfdpar,2))));
        tmp = sum( sum(bsxfun(@times, hessinfo.psfdiag, hessinfo.dfdpar(:,I,:).*hessinfo.dfdpar(:,J,:)),1) , 3);
        hessinfo.hesscoeff = mean(hessinfo.hesscoeff(:,:),2) + tmp(:);
    end;
    hessinfo = hessinfo.hesscoeff;
end;
meanhess = mean(hessinfo(:,:),2);
npar = floor(sqrt(2*size(meanhess,1))); % size(hessinfo,1) = npar * (npar+1)/2
mat = triu(ones(npar));
mat(logical(mat)) = meanhess;
mat2 = mat'; 
if nargin<3
    mat2(1:npar+1:end)=0;
    rdDM = ones(npar, size(hessinfo(:,:),2));
else
    dDG = reshape(full(diag(DG)),npar,[]);
    rdDM = 1./reshape(full(diag(DM)),npar,[]);
    scdDG = dDG.*rdDM.^2;
    mat2(1:npar+1:end)=mean(scdDG,2);
end;
M = mat+mat2;

info = 1;

pvec = symamd(M);
ddiag = diag(M);
mind = min(ddiag);
lambda = 0;
if mind < 0, 
  lambda = -mind + .001; 
end
lpcnt=1;
while (info > 0) && (lpcnt<10)
  M = M + lambda*eye(npar);
  [R,info] = chol(M(pvec,pvec));
  lambda = lambda*2 + 10;
  lpcnt = lpcnt+1;
end 
if info>0
    % don't precondition when we cant find an acceptable preconditioner.
    [R,info] = chol( eye(numel(pvec)));
end;
Rs.R = R;
Rs.pvec = pvec;
Rs.rdDM = rdDM;
mulfun = @compactedHessianMulPreconditioner;

function [w] = compactedHessianMulPreconditioner(x, R)

x = reshape(x, size(R.R,1),[]);
w =  R.R'\(x(R.pvec,:).*R.rdDM);
w(R.pvec,:) = R.R \ w;
w = reshape(w.*R.rdDM ,[],1);

function [w] = compactedHessianMulPreconditioner2(R,x)

x = reshape(x, size(R.R,1),[]);
w =  R.R'\(x(R.pvec,:).*R.rdDM);
w(R.pvec,:) = R.R \ w;
w = reshape(w.*R.rdDM ,[],1);