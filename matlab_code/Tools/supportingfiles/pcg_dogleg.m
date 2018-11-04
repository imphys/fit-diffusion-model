function [x, f, exflag, k] = pcg_dogleg( H_fun , H_info, precon_fun, precon_info, r, opts)
% [update_x, update_f, exflg, numits] = pcg_dogleg( H_fun , H_info, precon_fun, precon_info, G, opts );
% Computes:
%   update_x = arg min_{ norm(update) < maxR }  .5 * update'* H * update +  update' * G
%   update_f = .5 * update_x' * H  * update_x + update_x'*G
%
% with the preconditioned conjugate gradient method.
% The analytical solution is given by 
%   update_x = inv(H) * y
% constrained by 
%   norm( update_x ) < maxR
%
% INPUTS:
%  H_fun  : hessian multiply function (H * x = H_fun( H_info, x) )
%  H_info : hessian information passed to H_fun
%  precon_fun : preconditioner function, approximates inv(H)
%  precon_info: preconditioner information passed to precon_fun
%  G      : gradient
%  opts   : scalar structure with options
%      stop_tol_rel   : stop tolerance on residue norm, 
%      kmax           : maximum number of iterations
%      maxR           : trust region radius
%      maxR_frac_continue : stop conjugate gradient iterations when 
%                           norm(x) >= maxR * maxR_frac_continue 
%   Call otps = pcg_dogleg()  to obtain a default option structure.
%
% NOTE: H * x === H_fun( H_info, x) 
%       inv(H) * x  ==approx==  precon_fun( precon_info, x )
%
% OUTPUTS:
%   update_x : best value for update within trust region.
%   update_f : improvement in cost function value by update_x
%   exflag   : exit condition
%              0  : exhausted kmax iterations
%              1  : converged to within stop_tol(_rel)
%              2  : CG iterations stepped out of trust region, returned
%                   update at trust region boundary (maxR)
%              3  : stopped iterating since update close to trust region
%                   boundary. (maxR_frac_continue)
%              4  : negative curvature found, stepped to trust region
%                   boundary (maxR)
%  k         : number of iterations used.
% Created by Dirk Poot, Erasmus MC, 27-11-2012

if nargin<1
    opts.stop_tol_rel = .01;
    opts.kmax     = [];
    opts.maxR     = inf;
    opts.maxR_frac_continue = .95;
    opts.precon_first_step = true;
    x=opts;
    return;
end;

%% preconditioned conjugate gradient algorithm based on
% 'An introduction to the conjugate gradient method without the agonizing
% pain' by J.R. Shewchuk
% trust region dogleg extension based on 
% http://www.numerical.rl.ac.uk/nimg/course/lectures/raphael/lectures/lec7slides.pdf

x = zeros(size(r));
xx = 0;
f =0;
maxR2 = opts.maxR^2;
maxR2_continue = maxR2 * opts.maxR_frac_continue.^2;

% r = G;
if opts.precon_first_step
    d = precon_fun( precon_info, r );
else
    % dogleg method 'prefers' fist step in steepest descend direction.
    d = r;
end;
Mr = d;
rMr = r'*Mr;
stop_tol_r  = rMr * opts.stop_tol_rel.^2;
if isempty(opts.kmax)
    kmax = numel(r);
else
    kmax = opts.kmax;
    if kmax <= 0
        x = d; % assume perfect preconditioner used. 
    end;
end;
exflag = 0;
k=0;   
for k = 1 : kmax 
    Ad = H_fun( H_info, d );
    dAd = d'*Ad;
    if dAd<0 % negative curvature direction found:
        % dogleg step to maxR
        % solve  norm(x + alpha * d) = maxR  for alpha<0
        % => (x + alpha * d )'*(x + alpha * d ) == maxR^2
        % => x'*x + 2 * alpha * x'*d + alpha^2 * d'*d == maxR^2
        % => alpha = (- (x'*d) +/- sqrt( (d'*d) * maxR^2 + (x'*d)^2 -(d'*d)*(x'*x) ) )/ (d'*d)
        dd = d'*d;
        if k == 1 % first iteration x is scalar 0.
            xd = 0;
        else
            xd = x'*d;
        end;
        if isfinite(maxR2)
            det = max(0, dd * (maxR2-xx) + xd^2 ); 
            alpha = (-xd - sqrt(det))/dd;
        else
            % We dont have any (finite) trust region radius, 
            % So make a guess that negating the step required to get to 
            % the maximum is a good step.
            alpha = rMr/dAd; % note sign-reversed!
        end;
        if alpha>0 % DEBUG check
            error('alpha not in expected range');
        end;
        x = x + alpha * d;
        f = f + alpha * (r'*d) + .5*alpha^2*dAd ;
        exflag = 4;
        break;
    end;
    alpha = rMr/dAd;
    x_prev = x;
    x = x - alpha * d;
    xx = x'*x;
    if xx > maxR2
        % we stepped out of trust region
        %    solve  norm(x_prev - alpha_con * d) = maxR  for alpha > alpha_con > 0
        % (initially we solved a modified version:  solve  norm(x - alpha_mod * d) = maxR  for 0 > alpha_mod > -alpha
        %  but that may have numerical problems)
        %
        %  x_prev'*x_prev - 2* alpha_con * x_prev'*d + alpha_con^2 * d'*d == maxR2
        %  det = d'*d (maxR2 - xprev'*xprev) + (xprev*d)^2
        %  alpha = (xprev*d +/- sqrt(det))/(d'*d)
        if 0 
            % previous version
            dd = d'*d;
            xd = x'*d;
            det = max(0, dd * (maxR2-xx) + xd^2 ); 
            alpha_mod = (xd + sqrt(det))/dd;
            if alpha_mod>0 || alpha_mod<-alpha  % DEBUG check
                error('alpha_mod not in expected range');
            end;
            x = x - alpha_mod * d;
            f = f - (alpha_mod+alpha) * (r'*d) + .5*(alpha_mod+alpha)^2*dAd ;
        else
            clear x
            dd = d'*d;
            xpd = x_prev'*d;
            xp2 = x_prev'*x_prev;
            det = max(0, dd * (maxR2-xp2) + xpd^2 ); 
            alpha_mod = (xpd + sqrt(det))/dd;
            if alpha_mod<0 || alpha_mod>alpha  % DEBUG check
                error('alpha_mod not in expected range');
            end;
            x = x_prev - alpha_mod * d;
            f = f - (alpha_mod) * (r'*d) + .5*(alpha_mod)^2*dAd ;
        end;
        exflag = 2;
        break;
    else
        clear x_prev;
    end;
    % f(x) = .5 x' H x + x' G
    % f( x + alpha * d)  = .5 * x'* H * x + x'* G + alpha * x'* H *d + .5*alpha^2 * d'*H*d + alpha * d*G;
    %  x'*H *d =0 (due to construction of d)
    % d*G == d*r
    % => f += alpha^2 * dHd + alpha * dr
    f = f - alpha * (r'*d) + .5*alpha^2*dAd ;
    
    if xx >= maxR2_continue
        % we stepped so close to trust region boundary that it's not worth
        % continuing. 
        exflag = 3;
        break;
    end;
    r = r - alpha * Ad;
    Mr = precon_fun( precon_info, r );
    rMr_prev = rMr;
    rMr = r'*Mr;
    if rMr < stop_tol_r 
        exflag = 1;
        break;
    end;
    beta = rMr/rMr_prev;
    d = Mr + beta * d;
end;


%%
function tests
%% Test with random posdef matrix
opts = pcg_dogleg();
hmul = @(h,x) h*x;
a = randn(10,5);
% A = a'*a;
% b = randn(5,1);
A = [  31.981061182750981  -3.025380703368191   8.666283746422723 -17.066243985917485  -6.845750641896006
  -3.025380703368191   3.253257090829214  -4.723788830211969  -1.014117873648843   0.282712758740409
   8.666283746422723  -4.723788830211969  15.415235730342490  -2.604283635100411   4.090892855053212
 -17.066243985917485  -1.014117873648843  -2.604283635100411  21.023265497994799   1.922816777234540
  -6.845750641896006   0.282712758740409   4.090892855053212   1.922816777234540  15.641653266790021];
b =[-1.068848667544271
  -1.604978265317102
   2.959097104323258
   0.880916721900641
  -0.721948759481716];
precon = eye(size(A));
cfun = @(x) .5*x'*A*x+x'*b;
%% exflag = 0
opts.kmax = 1;
[x0, f0, exflag0, k0] = pcg_dogleg( hmul , A, hmul, precon, b, opts);
cfun(x0)-f0
%% exflag = 1
opts.kmax = [];
[x1, f1, exflag1, k1] = pcg_dogleg( hmul , A, hmul, precon, b, opts);
cfun(x1)-f1
%% exflag = 2
opts.maxR = norm(x1)*.99;
[x2, f2, exflag2, k2] = pcg_dogleg( hmul , A, hmul, precon, b, opts);
cfun(x2)-f2
norm(x2)-opts.maxR
%% exflag = 3
opts.maxR = norm(x1)*.90;
[x3, f3, exflag3, k3] = pcg_dogleg( hmul , A, hmul, precon, b, opts);
cfun(x3)-f3
norm(x3)-opts.maxR
%% exflag = 4;
B = A;B(2,2)=-1;
opts.maxR = 5*norm(x1);
% b = [0 0 0 -1 0]';
[x4, f4, exflag4, k4] = pcg_dogleg( hmul , B, hmul, precon, b, opts);
cfunB = @(x) .5*x'*B*x+x'*b;
cfunB(x4)-f4
norm(x4)-opts.maxR
grfun = gradest(cfunB, x4);
%% perfect preconditioning, step internal
iA = inv(A);
[x, f, exflag, k] = pcg_dogleg( hmul , A, hmul, iA, b, opts);
cfun(x)-f
gradest(cfun, x)'
%% good preconditioning, step external
A = diag((1:100));
B = diag(1./((2:101)/2));
b= randn(100,1);
%%
[x, f, exflag, k] = pcg_dogleg( hmul , A, hmul, B, b, opts);
cfunC = @(x) .5*x'*A*x+x'*b;
grest = gradest(cfunC, x)';
plot(grest)

