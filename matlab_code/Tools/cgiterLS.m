function [x, crit, y_fit, trace]=cgiterLS(mulA, mulAt, y, x_init, K , maxlps, saveiters, testA)
% [x, crit, y_fit, trace] = cgiterLS( A, At, y, x_init, K, maxlps, saveiters, testA)
% CG least squares iteration
% Solves y = A*x for x, in least squares sense:
%    x = arg min_x  (y-A*x)'*(y-A*x) =  arg min_x sum( abs(y-A*x).^2 )
% (computes x = inv(A'*A) * A' * y, without computing A'*A or the inverse)
% When a regularization matrix K is provided, the following regularized
% least squares problem is solved: 
%   x = arg min (A*x - y)'*(A*x - y) + x' * K * x
% (which has the general solution x = inv(A'*A + K) * A' * y)
%
% INPUTS:
% A     : sparse or full design matrix.
%         or function mulA to compute the matrix vector multiplication:
%             mulA(x) == A * x
% At    : should be empty when A is matrix. When A is a function, At should
%         be a function mulAt that computes the matrix vector multiplication:
%             mulAt(y) == (vec_y' * A)' == A'*vec_y
% y     : column vector with measurements
% x_init: optional intialisation of x, specify close to the optimal x for
%         faster convergence. As the conjugate gradient method is
%         guaranteed to converge to the true solution (when the condition
%         number of A'*A is low enough), no initialisation is required.
% K     : optional symmetric, positive (semi) definite regularisation matrix;
%         sparse, full, or a function that evaluates: 
%             mulK(x) = K * x
%            e.g.: mulK = @(vec) mulRegularisationnrm(vec, lambda, sz )
% maxlps: optional manual limit on the number of iterations.
% saveiters: default: false; optional save of the result in each iteration.
% testA : If positive integer: Test if At really is A' at testA
%         rows&columns, usefull for debugging when A and At are (different)
%         functions. (As bugs might intoduce differences). The testA rows
%         and columns are selected randomly, so the test might not be
%         appropriate for large sparse A (& At) unless testA is also large.
%
% NOTE: an old interface, currently supported for backward
%     compatibility only, allowed A to be a function in which the first 
%     argument selected multiplication with A or At:
%     mulfun(1, vec_x) = A * vec_x
%     mulfun(2, vec_y) = (vec_y' * A)' = A'*vec_y
%     please update call from:   cgiterLS( mulfun, ...)   % OLD INTERFACE
%                          to:   cgiterLS( @(x) mulfun(1, x), @(x) mulfun(2, x), ...)
%
% OUTPUTS:
% x     : approximation to the least squares solution of Y = A*x
% crit  : 4 element criterium vector, 
%          [ log([gApprx;x_updnrm/xnrm;rho]); -k]
% y_fit : Explained part of Y (so A*x)
% trace : matrix with 4 columns, stores crit at each iteration.
%
% Created by Dirk Poot, University of Antwerp.
% Modified 31-7-2009

% (A'*A)*x = (A'*y);    b = A'* y 
% solve for x
%
% This function computes: x = conjgrad(A'*A+K, A' * y, x_init) 
% with the general conjungated gradients function:
% function [x] = conjgrad(Q,b,x0)
%    r = b - Q*x0;
%    w = -r;
%    z = Q*w;
%    a = (r'*w)/(w'*z);
%    x = x0 + a*w;
%    B = 0;
%    for i = 1:size(Q)(1);
%       r = r - a*z;
%       if( norm(r) < 1e-10 )
%            break;
%       end if
%       B = (r'*z)/(w'*z);
%       w = -r + B*w;
%       z = Q*w;
%       a = (r'*w)/(w'*z);
%       x = x + a*w;
%    end
% end
%
% Which can be substituted to:
% (w -> -d; wrap iteration around; )
% function [x] = conjgradLS(A,y,x0)
%    r = A'*y - (A'*A+K)*x0;
%    d = r;
%    rho = r'*d;
%    x = x0;
%    B = 0;
%    for i = 1:size(A,2);
%       s = A * d;
%       ss = s'*s
%       a = rho/(ss +d'*K*d);
%       x = x - a*d;
%       z = -(A'*s + K*d);
%       r = r - a*z;
%       if( norm(r) < 1e-10 )
%            break;
%       end if
%       B = -(r'*z)/(d'*z);
%       d = r + B*d;
%       rho = r'*d;
%    end
% end
verbose = 1; % 0 : no progressbar
             % 1 : progressbar, no text
             
% Create functions to multiply with the specified matrix.
if ~(isa( mulAt,'function_handle') || (isempty(mulAt) && ~isa( mulA,'function_handle')))
    warning('CGiterLS:OldParam','Old order of input parameters detected, please update caller to new interface (see  "help cgiterLS")');
    % Shift all arguments to their correct name:
    if nargin>=7
        testA = saveiters; % 8
    end;
    if nargin>=6
        saveiters = maxlps; %7 
    end;
    if nargin>=5
        maxlps = K; %6
    end;
    if nargin>=4
        K = x_init; %5
    end;
    if nargin>=3
        x_init =y; %4
    end;
    if nargin>=2
        y = mulAt; %3
    end;
    if isa(mulA,'function_handle')
        A = mulA;
        mulA  = @(vect) A(1,vect);
        mulAt = @(vect) A(2,vect);
    else
        mulAt = [];
    end;
    nargshift = 1;
else
    nargshift = 0;
end;
if ~isa( mulA,'function_handle') && isempty(mulAt)
    A = mulA; mulA = [];
    hasAt = islogical(A) && issparse(A) && (nnz(A)<1e8); %false;
    mulA  = @(vect) A*vect;
    if hasAt
        Atransp = A';
        mulAt  = @(vect) Atransp*vect;
    else
        mulAt  = @(vect) (vect'* A)';
    end;
end;

% read inputs and set defaults:
if nargin<5-nargshift 
    K = [];
end;
if nargin<6-nargshift || isempty(maxlps)
    maxlps = 15;%min(3*size(r,1), ceil(100 * log(2+size(r,1))));
end;
if nargin<7-nargshift || isempty(saveiters)
    saveiters = false;
end;

% initialisation & constants:
k = 0;
k_lastbaseupd = 0;
toDbl = isa(y,'long');
hasK =  ~isempty(K);
save_y_fit = nargout>=3;
dotrace = nargout>=4;


if isa(K,'function_handle')
    mulK  = K;
else
    if isempty(K)
        mulK  = @(vect) 0;
    else
        mulK  = @(vect) K*vect;
    end;
end;
if nargin>=8-nargshift && testA>0
    szy = size(y);
    testAorig = testA;
    if numel(y)/2<testA
        ry = (1:min(testA,numel(y)))';
        testA = numel(ry);
    else
        ry = sort(ceil(rand(testA,1)*numel(y)));
    end;
    %ry = 100^3 + bsxfun(@plus,(20:25)',90*bsxfun(@plus,(20:25), 60*permute(20:25,[1 3 2])));numel(ry),ry= ry(:);
	storeAt = zeros(testA,testA);
    storeA = zeros(testA,testA);
    progressbar('start',numel(ry));
    for k=1:numel(ry)
        tst = zeros(size(y));
        tst(ry(k))=1;
        r = mulAt(tst);
        if k==1
            sx = size(r);
            if numel(r)/2<testAorig
                rx = min((1:testA)',numel(r));
            else
                rx = sort(ceil(rand(testA,1)*numel(r)));
            end;
            %rx = find(abs(r(:))>.53); % tuned to be testA
        elseif sx ~= size(r)
            error('size not constant');
        end;
        storeAt(:,k)= r(rx);
        tst = zeros(size(r));
        tst(rx(k))=1;
        r = mulA(tst);
        if size(r)~=szy
            error('sizes not consistant');
        end;
        storeA(:,k)= r(ry);
        progressbar(k);
    end;
    progressbar('ready');
    storeAt = storeAt';
    imagebrowse(abs(cat(3,storeA,storeAt,storeA-storeAt,storeA-storeAt)))
% imagebrowse(reshape(abs(cat(3,storeA,storeAt',storeA-storeAt',storeA-storeAt')),[18 12 12 testA 4]))    
    disp('inspect image, to see if A and At are equal enough.')
    disp(['Maximum absolute difference: ' num2str(max(max(abs(storeA-storeAt)))) ' maximum absolute values : ' num2str([max(max( storeA)) max(max(storeAt))])]);
    keyboard
%     storeA,storeAt
%     storeA-storeAt
end;
start = [ 0; 0; 0;   0.5   ];
limit = [-1;-1;-1; -maxlps ];
if verbose>0
    progressbar('start',[start limit],[],'mintimeinterval',1,'esttimeleft','on');
end;
% initialize residue (with x_init or zeros)
if nargin>=3 && ~isempty(x_init)
    start(4) = 1;
    if verbose>0
        progressbar('adjustlimits',[start limit]);
    end;
    s = mulA(x_init);
    if verbose>0
        progressbar([0;0;0;.5])
    end;
    if save_y_fit
        y_fit = s;
    end;
    r = mulAt( y - s );
    if hasK
        r = r - mulK(x_init);
    end;
    x = x_init;
    xnrm = sum(x_init.^2);
else
    r = mulAt(y);
    x = zeros(size(r,1),1);
    xnrm = 0;
    if save_y_fit
        y_fit = 0;
    end;
end;
if toDbl
	r = double(r);
end;
if verbose>0
    progressbar([0;0;0;0])
end;

d = r;
rho = r'*r;
if dotrace
    trace = zeros(min(ceil(maxlps/10),200),4);
    trace(1,:)=[nan nan log(r'*r) -k];
    traceInd = 2;
end;
x_base = zeros(size(r,1),1);
if isa(y,'long')
    x_base = long(x_base);
    gOrig = double(y);
    gOrig = gOrig'*gOrig;
    gEps = double(eps(y));
    gEps = gEps'*gEps;
    if nargin<5
        maxlps = maxlps * longprecision*(log(10)/log(2^60));
    end;
else
    gOrig = (y'*y);
    gEps = eps^2*gOrig*3*sqrt(size(y,1));
    rEps = eps^2*rho*3*sqrt(size(y,1));
end;
gApprx = gOrig;
start(1:3) = log([gOrig ;   1  ; rho ]);
limit(1:3) = log([gEps  ; eps^2; rEps]);
if verbose>0
    progressbar('adjustlimits',[start limit]);
end;
rhoThresh = rho/(256^6); % first time really large improvement needed.
kTresh = ceil(size(r,1)*.51);
crit = start;
while all(crit>limit)
    k = k+1;
    s = mulA( d );
    crit(4) = -k+.5;
    if verbose>0
        progressbar(crit);
    end;
    ss = s'*s;
    if hasK
        Kd = mulK(d);
        ss = ss + d'*Kd;
    end;
    if ss==0
        break;
    end;
    a = rho./ss;
    
    x = x + a*d;
    if save_y_fit
        y_fit = y_fit + a * s;
    end;
    if k==1
        xnrm = sum(x.^2);
    end;
    x_updnrm = a^2 .* sum(d.^2);
    if saveiters
        save(['cgiterLS_iter' num2str(k)],'x');
    end;
%     disp(['Loop ' num2str(k) ', rho = ' num2str(rho) ', gApprx = ' num2str(gApprx) ', |x_upd| = ' num2str(a*sqrt())]);
    if (k-k_lastbaseupd>kTresh) || (rho<rhoThresh*sqrt(k-k_lastbaseupd))
        k_lastbaseupd = k;
        disp('subtracting base for further precision improvement.');
        % recompute residual from original y (instead of A*y)
        if toDbl && islogical(A)
            sc = log2(max(abs(x))/bitmax*max(size(A)));
            x_ad = pow2(round(pow2(x,-sc)),sc);
            x_base = x_base + x_ad;
            x = x - x_ad;
            y = y - mulA( x_ad );
            ry = y;
        else
            % accurately update x_base (and set x to zero):
            x_b_old = x_base;
            x_base = x_base + x;
            x = x + (x_b_old - x_base) ;
            xnrm = sum(x_base.^2);
            if toDbl
                x = double(x);
                xnrm = double(xnrm);
            end;
            % compute unexplained part of y:
            y_fit = mulA( x_base );
            ry = (y - y_fit);
            if nnz(x)~=0
                s = mulA(x);
                if save_y_fit
                    y_fit = y_fit + s;
                end;
                ry = ry - s;
            end;
        end;
        if toDbl
            ry = double(ry);
        end;
        gApprx = ry'*ry;
        r = mulAt(ry);
        if hasK
            % r = A'*(y - A*x) - K*x;
            r = (r - mulK(x_base)) - mulK(x);
        end;
        if toDbl
            r = double(r);
        end;
        rhoThresh = rho/256^3;
    else
        z = mulAt( s ) ;
        if hasK
            z = z + Kd;
        end;
        r = r - a*z;
    end;
    rhoOld = rho;
    rho = r'*r;
    gam = rho/rhoOld;
    d = r + gam*d;
%     pause(.5);
    crit = [log([gApprx;x_updnrm/xnrm;rho]);-k];
    if dotrace
        if traceInd>size(trace,1)
            trace = [trace; zeros(size(trace,1),size(trace,2))];
        end;
        trace(traceInd,:) = crit;
        traceInd = traceInd+1;
    end;
    if verbose>0
        progressbar(crit);
    end;
% %     plot(1:length(x),[x ]);
%     if ~exist('xPrep','var')
%         xPrep = ones(size(x));
%     end;
%     plot(1:length(x),[x_base-xPrep+x]); %pause(1);
end;
if verbose>0
    progressbar('ready');
end;
x = x+ x_base;
if dotrace
	trace(traceInd:end,:) = [];
end;
