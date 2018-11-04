function [ x , fval, final_grad ] = conjgrad_nonlin( funcs , x , varargin)
% [xout, fval, final_grad] = conjgrad_nonlin( funcs , xinit , [option, value]);
%
% Non linear conjugate gradient method.
% Best if hessian is always positive definite.
%
% Implementation based on section B4 of 
% "An Introduction to the Conjugate Gradient Method Without the Agonizing Pain.pdf"
%
% INPUTS:
% funcs : a 2 element structure that specifies two functions, the sum of
%         which is minimized.
%         The structure should contain the fields:
%    fun     : the actual functions; format as for fminunc:
%               [fval , grad, hess_info] = fun( x );
%              third output is hessian information used for:
%    hessmul : an function that multiplies with the hessian.
%               H*x = hessmul( hess_info , x);
%              (same format as for fminunc)
%     funcs(1) : should specify a quadratic function 
%     funcs(2) : should specify a non-quadratic function. 
%
% option-value pairs or scalar structure with options; valid options:
%  epsilon : scalar convergence tolerance. When gradient norm reduces by at least 
%            this factor from it's initial value, optimization stops. 
%  max_iter  : maximum number of (outer) iterations
%  revaluate_quadratic_it : Default = 50; number of iterations after which the 
%            quadratic function is revaluated. Can be used to avoid accumulation 
%            of roundoff error (default, although higher might be acceptable as
%            well). When the 'quadratic' function is not perfectly quadratic 
%            this should be set to a lower value. 
%  revaluate_quadratic_hessian : boolean; default=true; Update hessian
%            information when revaluating the quadratic function.
%  searchHistory : default 0; 
%            Number of previous search directions that are cached and used
%            in the non linear search part. Set to 5 to 20 when the
%            nonlinear function evaluates much quicker (CPU time) than a
%            hessian multiplication of the quadratic function.
%  othogonalizationmode : 0 = Fletcher-Reeves search vector update
%                        -1 = Polak-Ribiere search vector update
%                        >0 = explicit orthogonalization w.r.t. previous
%                             search vectors, requires searchHistory > 0.
%                             This typically reduces the number of
%                             iterations at the expense of extra hessian
%                             multiplication of the non linear part.
%  collect_statistics : default false => fval is 2 element vector with
%            final function value of quadratic and non linear parts.
%            When true: fval is a structure with information at each iteration.
%  
% OUTPUTS:
%  xout : best x that I found.
%  fval : function cost at xout
%  final_grad : gradient of cost function at xout.
%
% Created by Dirk Poot, Erasmus MC, 31-10-2012
%       May 2013: Added searchHistory and orthogonalization mode

% Constants/initialisation:
options.epsilon = .0001; % convergence tolerance.
options.max_iter  = 200; % maximum number of iterations.
options.max_NR = 3;     % maximum number of newton-raphson iterations.
options.epsilon_NewtRaphs = .1; % Line search tolerance.
options.revaluate_quadratic_it = 50; % number of iterations after which the function value and gradient of the quadratic part are re-evaluated.
options.revaluate_quadratic_hessian = true;
options.searchHistory = 0; % search history cache. Requires extra memory and quadratic 
                           % part to be much more computationally expensive than non linear 
                           % part.
options.othogonalizationmode = 0; 
options.collect_statistics = false;
options.display = 3;  % 0 : no display of progress, 
                      % 1 : final summary (not implemented yet)
                      % 2 : display each iteration.

if nargin>=3
    options = parse_defaults_optionvaluepairs( options, varargin{:});
end;

szx = size(x); 
n = numel(x);

% parse inputs:
fun_quad = funcs(1).fun;
hessmul_q = funcs(1).hessmul;
fun_nonl = funcs(2).fun;
hessmul_n = funcs(2).hessmul;

if options.searchHistory>0
    opt = fmin_fast();
    opt.HessMult = @(h,x) subspace_fun_hessmul( h, x , hessmul_n);
    opt.Preconditioner = @(h,x) [];
    opt.Preconditioner_Multiply = @(dum, x) x ;

end;

% start of algorithm
it =0 ;
k = 0;
[f_q, g_q, h_q_info] = fun_quad( x );
g_q = g_q(:);
[f_n, g_n, h_n_info] = fun_nonl( x );
old_f_q = [];

r = -(g_q + g_n(:));
d = r;
delta_0 = r'*r;
delta_new = delta_0;
if options.othogonalizationmode<0
    r_old =r ;
end;
if options.collect_statistics
    % define statistics structure; update fields in each loop
    numstats =1;
    stats.iter = it;         % ok
    stats.f_q = f_q;         % ok
    stats.f_q_approx = f_q;  % ok
    stats.f_n = f_n;         % ok
    stats.num_f_q = 1;       % ok
    stats.num_f_n = 1;       % ok
    stats.num_hessmul_q = 0; % ok
    stats.num_hessmul_n = 0; % ok
end;
if options.display>1    
    fprintf(' Iteration   f(x)    |update x|      |d(f(x))|   [predicted f(x)]\n');
    fprintf('%5d   %12.9g           %15.3g\n',it, f_q+f_n , sqrt(delta_0) );
end;
while it < options.max_iter && delta_new > options.epsilon^2*delta_0
    it = it+1;
    
    delta_d = d'*d;
    
    Hd_q =  hessmul_q( h_q_info, d);
    dHd_q = d'*Hd_q;
    
    % first 1D search step: (this is the standard conjugate gradient step)
    Hd_n = hessmul_n( h_n_info, d);
    dHd = d'*Hd_n + dHd_q;
    alpha = -((g_q+g_n(:))'*d)/( dHd);
    if options.collect_statistics
        numstats = numstats + 1;
        stats(numstats).iter = it;
        stats(numstats).num_f_q = 0;
        stats(numstats).num_hessmul_q = 1;
        stats(numstats).num_hessmul_n = 1;
    end;

    if options.searchHistory==0 
        % no history
        % check if this step reduced the gradient enough:
        gd_q = g_q'*d;
        j=0;numhessmul = 0;
        while true 
            % f = .5 * x H x + b x + c
            % df/dx = H x + b
            % df/dx | (x -> x + alpha*d) - df/dx | (x -> x ) = H * (alpha*d)
            clear g_n h_n_info % save memory
            x_n = x + alpha * reshape(d, szx);
            [f_n, g_n, h_n_info] = fun_nonl( x_n );
            j = j + 1;
            doNewtRaphs = (j<options.max_NR) &&  (alpha^2*delta_d > options.epsilon_NewtRaphs);
            if ~doNewtRaphs 
                break;
            end;
            Hd_n = hessmul_n( h_n_info, d);numhessmul = numhessmul +1;
            alpha = alpha -((g_q+alpha * Hd_q+g_n(:))'*d)/( d'*Hd_n + dHd_q);
        end;
        x = x_n;
        f_q = f_q + gd_q * alpha + .5*alpha^2 * dHd_q;
        g_q = g_q + alpha * Hd_q; % update gradient of quadratic part.
        x_step_norm = sqrt(delta_d)*alpha ; 
        if options.collect_statistics
            stats(numstats).f_q_approx = f_q;
            stats(numstats).f_n = f_n;
            stats(numstats).num_f_n       =  j;
            stats(numstats).num_hessmul_n = stats(numstats).num_hessmul_n + numhessmul;
        end;
    else
        % use history and check if we need to proceed the subspace search.
        % 
        % x_new = x + D * beta
        % beta^0(hid) = alpha
        % Gradient in subspace:
        % f( x + D * beta) = f_q( x + D * beta) + f_n( x + D * beta)
        %                  = f_q(x) + g_q' * D * beta + .5* beta'* D' * H_q * D * beta + f_n( x + D * beta)
        % D[ f , beta] = g_q' * D +  DHD_q * beta + (f_n'[x+D*beta])^T * D
        % H[ f , beta] = DHD_q +  D^T * (f_n''[x+D*beta]) * D
        hid = mod(it-1, options.searchHistory)+1;
        dscale = 1./sqrt(dHd); % scale search direction so subspace hessian becomes identity (deviates only due to nonlinearities)
        if it==1
            % initialize variables:
            D = d* dscale;
            nD = 1;
            HD_q = Hd_q *dscale;
            DHD_q = (d'*HD_q) * dscale;
        else
            D(:,hid) = d * dscale;
            nD = size(D,2);
            HD_q(:,hid) = Hd_q *dscale;
            DHD_q(hid,1:nD) = (d'*HD_q) * dscale;
            DHD_q(1:nD,hid) = DHD_q(hid,:)';
        end;
        gD_q = g_q'*D;
        beta_in = zeros(nD,1);
        beta_in(hid)= alpha/dscale;
        fun = @(bta) subspace_fun( bta, DHD_q, gD_q, fun_nonl , D, x);
        [beta, fsub_opt, exflg, gsub_opt, hsub_ifo, fminfastoutp] = fmin_fast( fun, beta_in, opt );
        g_n = hsub_ifo.gn;
        h_n_info = hsub_ifo.hn;
        x_step = reshape( D * beta, szx);
        x = x + x_step;
        x_step_norm = sqrt(x_step(:)'*x_step(:));
        f_q_prev = f_q;
        f_q = f_q + gD_q * beta + .5*(beta'*DHD_q*beta);
        g_q = g_q + HD_q * beta; % update gradient of quadratic part.
        f_n = fsub_opt+(f_q_prev-f_q);

        if options.collect_statistics
            stats(numstats).f_q_approx = f_q;
            stats(numstats).f_n = f_n;
            stats(numstats).num_f_n       =  fminfastoutp.funcCount;
            stats(numstats).num_hessmul_n = stats(numstats).num_hessmul_n + fminfastoutp.cgiterations ;
        end;

    end;

    if mod(it, options.revaluate_quadratic_it ) == 0
        old_f_q = f_q;
        if options.revaluate_quadratic_hessian
            clear g_q h_q_info % save memory
            [f_q, g_q, h_q_info] = fun_quad( x );
            g_q = g_q(:);
        else
            clear g_q % save memory
            [f_q, g_q] = fun_quad( x );
            g_q = g_q(:);
        end;
        if options.collect_statistics
            stats(numstats).f_q = f_q; % don't set f_q if we dont explicitly compute it.
            stats(numstats).num_f_q = 1;
        end;
    end;
    
    r = -(g_q+g_n(:)); % new gradient
    delta_old = delta_new;
    delta_new = r'*r;
    if options.othogonalizationmode <= 0 
        if options.othogonalizationmode==0
            % Fletcher-Reeves search vector update
            beta = delta_new/delta_old;
        else
            % Polak-Ribiere search vector update:
            beta = max( (r'*(r- r_old))/delta_old ,0 );
            r_old = r;
        end;
        d = r + beta * d;
        k = k+1;
        restartCG = (k >= n ) || (r'*d<=0);
        if restartCG
            d = r;
            k=0;
        end;
    else
        % make search direction H-orthogonal to current search directions:
        % =>  d_new * H * D = 0 
        %   == (r + D * beta)' * H * D = 0
        %  => r' * H * D + beta' * D' * H * D = 0 
        %  => beta' = -(r'*H*D)* inv( D'*H*D )
        selcols = mod( hid-1:-1:hid-min(nD,options.othogonalizationmode) , nD)+1;
        HD_n = zeros( numel(d), numel(selcols));
        for kH = 1:numel(selcols)
            HD_n(:,kH) = hessmul_n( h_n_info, D(:,selcols(kH) ) );
        end;
        if options.collect_statistics
            stats(numstats).num_hessmul_n = stats(numstats).num_hessmul_n + numel(selcols);
        end;
        HD = HD_n + HD_q(:,selcols);
        clear HD_n
        beta =  ( -(r'*HD) / ( D(:,selcols)'*HD ) )';
        d = r + D(:,selcols)*beta;
    end;
    
    if options.display>1    
        if ~isempty(old_f_q)
            fprintf('%5d   %12.9g %9.4g %15.3g  %12.9g\n',it, f_q+f_n, x_step_norm, sqrt(delta_new) , old_f_q+f_n);
            old_f_q = [];
        else
            fprintf('%5d   %12.9g %9.4g %15.3g\n',it, f_q+f_n, x_step_norm , sqrt(delta_new) );
        end;
    end;
end;
if options.collect_statistics
    fval = stats;
else
%     fval = [f_q , f_n , fun_quad( x ) , fun_nonl( x )]; % DEBUG, explicit evaluation as well.
    fval = [f_q , f_n]; 
end;
if nargout>2
    final_grad = g_q+g_n(:);
end;

function [f, g, h] = subspace_fun( x, A, b, f_n , D, x0);
% [f,g,h] = quadratic_test_fun( x, A, b, f_n , D, x0);
%
% function that evaluates
% f = 0.5 * x'*A*x + b'*x + f_n( x0 + D * x )
% 
% and is able to evaluate the derivative and hessian as well. 
%
% Created By Dirk Poot, Erasmus MC, 28-5-2013

Ax = A*x;
f = 0.5*( x'*Ax ) + b*x;
sx = size(x0);
if nargout>1
    if nargout>2
        [fn ,gn, hn] = f_n( x0 + reshape(D* x,sx) );
        h.A = A;
        h.hn = hn;
        h.gn = gn;
        h.D = D;
    else
        [fn ,gn] = f_n( x0 + reshape(D* x,sx));
    end;
    g = Ax + (b + (gn(:)'*D))';
else
    fn = f_n( x0 + reshape(D* x,sx));
end;
f = f +fn;

function [hx] = subspace_fun_hessmul( h, x , f_n_hessmul)
hx = h.D' * f_n_hessmul( h.hn, h.D*x)  + h.A*x;

function test_simple_quadratic_function()
%%
A1 = randn(90,100);A1 = A1'*A1; 
A2 = randn(90,100);A2=A2'*A2; 
b1 =randn(100,1);b2=randn(100,1);
funcs(1).fun= @(x) quadratic_test_fun( x, A1,b1);
funcs(2).fun= @(x) quadratic_test_fun( x, A2,b2);
funcs(1).hessmul = @(H,x) H*x; 
funcs(2).hessmul = @(H,x) H*x;
xin = zeros(100,1);
[ x_h , fval ] = conjgrad_nonlin_v2( funcs , xin ,'searchHistory',5);
[ x_1 , fval ] = conjgrad_nonlin_v2( funcs , xin ,'searchHistory',0);
[ x_0 , fval ] = conjgrad_nonlin( funcs , xin ,'searchHistory',0);
plot([x_h x_0 x_1])
function test_TV_regul_function()
%%
ncon = 3000;
x_GT = phantom( 'Modified Shepp-Logan' , 100);
sz = size(x_GT);
A1 = randn(ncon,prod(sz))*20;A1 = A1'*A1; 
b1 = -A1 * x_GT(:);
funcs(1).fun= @(x) quadratic_test_fun( x(:), A1,b1);
funcs(2).fun= @(x) totalVariationVecRegularizer( x, [],5,[],[],1);
funcs(1).hessmul = @(H,x) H*x; 
funcs(2).hessmul = totalVariationVecRegularizer( [], [],[],5,[],2);
xin = randn([1,sz]);
%% 
% profile on
tic;
[ x , fval ] = conjgrad_nonlin_v2( funcs , xin ,'searchHistory',5,'max_iter',50,'epsilon', eps);
toc
tic;
[ x2 , fval2 ] = conjgrad_nonlin_v2( funcs , x ,'searchHistory',5,'max_iter',50);
toc
tic;
[ xL , fvalL ] = conjgrad_nonlin_v2( funcs , xin ,'searchHistory',5,'max_iter',100,'epsilon', eps);
toc
tic;
[ xL2 , fvalL2 ] = conjgrad_nonlin_v2( funcs , xin ,'searchHistory',5,'max_iter',1000,'epsilon', eps);
toc
tic;
[ xL3 , fvalL3 ] = conjgrad_nonlin_v2( funcs , xL2 ,'searchHistory',5,'max_iter',1000,'epsilon', eps);
toc
tic
[ x_1 , fval_1 ] = conjgrad_nonlin_v2( funcs , xin ,'searchHistory',0, 'max_iter',51,'epsilon', eps);
toc
tic
[ x_0 , fval_0 ] = conjgrad_nonlin( funcs , xin ,'searchHistory',0,'max_iter',51,'epsilon', eps);
toc
tic
[ x_0L , fval_0L ] = conjgrad_nonlin( funcs , xin ,'searchHistory',0,'max_iter',1100,'epsilon', eps);
toc
% profile viewer
%%
imagebrowse(cat(4,x_0, x_0L,x,x2,xL,xL2 ,xL3, reshape(x_GT,[1 sz])),[-.2 1.2])
q=[sum( fval_0), sum( fval_0L) , sum( fval) ,sum(fval2),sum(fval3),sum(fvalL),sum(fvalL2),sum(fvalL3)];
q-min(q)
% [sum(fval_0),sum(fval)]*[eye(2) [1;-1]]