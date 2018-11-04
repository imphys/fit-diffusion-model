function[RPCMTX,ppvec] = aprecon_full(A,upperbandw,DM,DG,varargin)
%APRECON Banded preconditioner function for least-squares problems.
%
%   [RPCMTX,PPVEC] = APRECON(A,UPPERBW,DM,DG) produces the sparse
%   nonsingular upper triangular matrix RPCMTX such that
%   RPCMTX'*RPCMTX is a preconditioner for the
%   matrix M = DM*(A'*A)*DM + DG, where DM is a positive
%   diagonal matrix, DG is a non-negative diagonal matrix,
%   and A is sparse rectangular matrix with more rows than columns.
%   PPVEC is the associated permutation (row) vector and UPPERBW
%   specifies the upperbandwidth of RPCMTX.

%   Default preconditioner for SLLSBOX and SNLS.

%   Copyright 1990-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2007/08/03 21:30:19 $

% Initialization

n = size(A,2);

% Form matrix M
TM = bsxfun(@times, A, full(diag(DM))');

% Determine factor of preconditioner.
   DDG = diag(sqrt(full(diag(DG)))); 
   TM = [TM;DDG];
   [dummy,RPCMTX] = qr(TM,0); 
%    RPCMTX = sparse(triu(RPCMTX)); % make sparse since then \ is faster (ARGH!!!)
   ppvec = 1:n;
   
   %    Modify for singularity?
   mdiag = min(abs(diag(RPCMTX)));
   lambda = 1;
   while mdiag < sqrt(eps);
      TM = [A*DM; DDG + lambda*eye(n)];
      lambda = 4*lambda;
      [dummy, RPCMTX] = qr(TM,0);
%       RPCMTX = RPCMTX(1:n,1:n);
%       ppvec = p;
      mdiag = min(abs(diag(RPCMTX)));
   end
   

