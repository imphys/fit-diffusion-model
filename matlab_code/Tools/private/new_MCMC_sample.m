function [xnew, fnew ] = new_MCMC_sample( fun, x , options)
% [xnew, fvalnew ] = new_MCMC_sample( fun, x , options)
% Creates a new sample from the log likelihood function specified by fun
% with the Metropolis–Hastings algorithm.
% 
% As proposal distribution (which is a free parameter) in this algorithm
% I use a multivariate normal distribution of which the log(pdf) is
% specified by the gradient and hessian of fun. Actually, to be more
% precise by the proposal distribution is derived from the preconditioner
% of the hessian; which might be different when the hessian is not given
% explicitly or not positive definite. 
%
% With this proposal distribution, very fast mixing should be obtained.
% (assuming the preconditioner is good)
%
% General info:
%  http://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm
%  http://en.wikipedia.org/wiki/Gibbs_sampling
%  http://en.wikipedia.org/wiki/Markov_Chain
%
% Created by Dirk Poot,
% Erasmus MC, 6-2-2013

if nargin<1
    % create anonymous functions with as empty as possible workspace:
    default_HessMult       = @(H, x) H*x;
    default_MultiplyRTinv  = @(R, x) (R' \ x);
    default_MultiplyRinv   = @(R, x) (R \ x);
    default_MultiplyR      = @(R, x) (R * x);
    default_logdet         = @(R)   (-2)*sum(log(diag(R)));
    
    %return default options:
    options.HessMult                     = default_HessMult;
    options.Preconditioner               = @make_advanced_preconditioner; % make preconditioner
    options.Preconditioner_Multiply      = @mul_default_preconditioner;
    options.Preconditioner_MultiplyR     = default_MultiplyR;
    options.Preconditioner_MultiplyRinv  = default_MultiplyRinv;
    options.Preconditioner_MultiplyRTinv = default_MultiplyRTinv;
    options.Preconditioner_logdet        = default_logdet;
    options.MCMCstdscale = 1; % scale factor for standard deviation of proposal distribution.
    options.maxIter = 1; % default only 1 iteration
    xnew=options;
    return;
end;

% evaluate log likelihood function:
[f, G, H_info ] = fun( x );

% proposal distribution is 
%  newtonstep = -H\G;
%  Normal( x + newtonstep, inv(H) ) 
%    = (2*pi)^(-k/2) * 1/sqrt( det( inv(H) ) ) * exp( -.5* (x-mu)' * inv( inv(H) ) * (x-mu) )
% 
% R' * R  = H
% then a sample from the proposal distribution is given by 
%  xnew = x+newtonstep + R\ r;
% with r = randn(size(x)); (that is : independed standard normal)

precon_info = options.Preconditioner( H_info , x);

for iter = 1 :options.maxIter
    % propose a random step:
    % proposal distribution = normal distribution whose log is given by
    % gradient and hessian at current position
    r = randn(numel(x),1);
    xstep0 = options.Preconditioner_MultiplyRTinv( precon_info, -G(:) );
    xnew = x + reshape( options.Preconditioner_MultiplyRinv( precon_info, xstep0 + options.MCMCstdscale * r ) , size(x));

    [fnew, Gnew, H_infonew ] = fun( xnew );
    precon_infonew = options.Preconditioner( H_infonew , xnew);

    % Metropolis–Hastings algorithm
    % a = P( xnew) / P( x ) * Q( x | xnew) / Q( xnew | x)
    % with 
    %    P is actual likelihood function
    %    Q = proposal distribution.

    Q_xnew_x = - sum(r.^2)/2 - .5 * options.Preconditioner_logdet( precon_info ) ;

    xstep_new = options.Preconditioner_Multiply( precon_infonew, -Gnew(:) );
    r_new = options.Preconditioner_MultiplyR( precon_infonew, xnew(:) + xstep_new - x(:) );
    Q_x_xnew = - sum(r_new.^2)/2 - .5 * options.Preconditioner_logdet( precon_infonew ) ;

    % fun specifies the -log likelihood, so negate sign of f and fnew:
    log_a = -(fnew - f) + (Q_x_xnew - Q_xnew_x);

    if log_a<0
        reject = rand(1) > exp(log_a);
        if reject
            xnew = x;
            fnew = f;
        end;
    else 
        reject = false;
    end;
    if ~reject && iter < options.maxIter
        % set new as current position and update all variables dependent on
        % the current position:
        f = fnew;
        x = xnew;
        G = Gnew;
        H_info = H_infonew;
        precon_info = precon_infonew;
    end;
end;