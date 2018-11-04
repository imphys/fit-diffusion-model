function [x, f, exflag, G, H_info, outp] = fmin_fast( fun, x, options , constraints)
% [x, f, exflag, G, H, outp] = fmin_fast( fun, xinit, options , constraints)
% Fast non linear optimization routine.
% 
% INPUTS:
%  fun : function that should be minimized. Called as:
%         [fx, Gx, Hx ] = fun( x )
%       where fx is the function value at the location x 
%             Gx is the derivative of fun at x
%             Hx is the hessian at x, or hessian information when HessMult is
%                 provided
%  xinit: vector, matrix or ND array with initial position.
%  options: scalar structure with options; obtain the default values with 
%           options = fmin_fast();
%    Some of the valid options:
%   HessMult : hessian multiplication function
%               HessMult( Hx, y ) should multipy the vector y with the
%               hessian evaluated at x
%   Preconditioner : preconditioner creation function
%               precon_info = fun_makeprecon( H_info , x);
%   Preconditioner_Multiply : function that multiplies with the
%               preconditioner.
%               d = precon_fun( precon_info, y );
%               approximates  inv(Hx) * y
%                
% OUTPUTS:
%  x      : Optimal location of fun that I found
%  f      : optimal value of fun at x
%  exflag : exit flag; 
%            0 : run out of iterations; not converged.
%            1 : norm of gradient satisfies convergence tolerance (norm(G) < abs_tol_G)
%            2 : norm(step in x) <  abs_tol_step_x
%            3 : converged. Function was reasonably approximated by a
%                quadratic function and the improvement in f was small.
%  G      : Gradient at x (should be close to 0)
%  H_info : Hessian (information) at x
%  outp   : extra information about the optimization; 
%            iterations  : maximum number of iterations or function calls
%                          reached.
%            funcCount   : number of function evaluations
%            cgiterations: Total number of conjugate gradient iterations,
%                          which equals the number of hessian
%                          multiplications that are performed.
%            stepstaken  : number of steps that are accepted. (Any
%                          difference with iterations indicates
%                          'wasted' computations.)
%            preconditionerPrepare : number of times the preconditioner is
%                          constructed.
% 
% Created by Dirk Poot, Erasmus MC, 21-11-2012

% error('not finished yet')
if nargin<1
    %return default options:
    options.HessMult = @(H, x) H*x;
    options.Preconditioner = @make_default_preconditioner; % make preconditioner
    options.Preconditioner_Multiply = @mul_default_preconditioner;
    options.TR_scale_invalid_f = .1;
    options.TR_scale_good_step = 2;
    options.TR_scale_bad_step  = .25;
    options.TR_scale_verybad_step = .1;
    %options.Preconditioner_Multiply = ...;%
    options.InitialTrustRegionRadius = [];
    options.pcg_options = pcg_dogleg();
    options.maxIter = 100;
    options.abs_tol_G = .0001;
    options.abs_tol_x = .001;
    options.MaxFunEvals = 100;
    x=options;
    return;
end;

size_x = size(x);

[f, G, H_info ] = fun( x );numFunEvals = 1;
H_updated = true;
fun_hessmul    = options.HessMult;
fun_makeprecon = options.Preconditioner;
fun_mulprecon  = options.Preconditioner_Multiply;
pcg_opts = options.pcg_options;
pcg_opts.maxR = options.InitialTrustRegionRadius;
abs_tol_G          = options.abs_tol_G;
abs_tol_step_x     = options.abs_tol_x;
max_numFunEvals    = options.MaxFunEvals;
preconditionerPrepare =0;
cgiterations = 0;
numstepsaccepted = 0;
posdef = 1;
exflag = 0;

for iter = 1 : options.maxIter
    % checks value and gradient:
    if ~isfinite(f) || any(~isfinite(G(:)))
        error('Function value and gradient should be finite.');
    end
    
    % Check convergence:
    G_norm = norm( G(:) );
    if posdef && (G_norm < abs_tol_G )
        % norm of gradient satisfies convergence tolerance.
        exflag = 1;
        break;
    elseif iter>1
        if ( pcg_exflag <=1 ) && ... %(x_step_norm < pcg_opts.maxR_frac_continue * pcg_opts.maxR) && ...
           (quadratic_approx_ratio > .25) && ...
           (-delta_f < abs_tol_G * (1+abs( f )) )
            exflag = 3;
            break;
        elseif (x_step_norm < abs_tol_step_x )
            exflag = 2;
            break;
        elseif (numFunEvals > max_numFunEvals)
            exflag = 0;
            break;
        end;
       
    end
    
    % Compute Newton step, limitted to within trust region:
    if H_updated 
        precon_info = fun_makeprecon( H_info , x);
        H_updated = false; % preconditioner and hessian are now consistent
        preconditionerPrepare = preconditionerPrepare +1;
    end;
    [x_step, delta_f_predicted, pcg_exflag, k] = pcg_dogleg( fun_hessmul , H_info, fun_mulprecon, precon_info, G(:), pcg_opts);
    cgiterations = cgiterations + k;
    posdef = (pcg_exflag ~= 4);
    
    x_step_norm = norm(x_step);
    
    x_new = x + reshape(x_step, size_x);
    
    [f_new,G_new ,H_info_new ] = fun( x_new );numFunEvals = numFunEvals + 1;
    % ss = linspace(0,1,1000);ft=ss;fp=ft;g = x_step(:)'*G(:);H=x_step(:)'*fun_hessmul(H_info,x_step);for k=1:numel(ss);ft(k) = fun(x + ss(k)*reshape(x_step, size_x));fp(k) = ft(1)+ss(k)*g+.5*ss(k).^2*H;end;plot(ss,[ft' fp'] )
    delta_f = f_new - f;
    if ~isfinite(f_new)
        if isempty( pcg_opts.maxR) % trust region radius not specified. 
            pcg_opts.maxR = inf; 
        end;
        pcg_opts.maxR = min( x_step_norm, pcg_opts.maxR) * options.TR_scale_invalid_f;
        quadratic_approx_ratio = 0;
    else
        quadratic_approx_ratio = delta_f/delta_f_predicted;
        if isempty(pcg_opts.maxR) % first iteration when trust region radius not specified
            if quadratic_approx_ratio<.75
                % quadratic approximation not good
                
                % fit f + alpha * g + alhpa^2 * h + |g|*(exp( b* alpha) -1 -b*a -.5*b^2*a^2)
                % for fun( x + alpha * x_step )
                % on function value, gradient and hessian of fun( x ) 
                % and function value of fun( x + 1*x_step )
                g_f = G(:)'*x_step(:);
                h_f = -g_f;

                % solve b from f_new-(f+g_f+h_f/2) = |g|*(exp(b)-1-b- .5* b^2 )
                lhs = (f_new - (f + g_f + .5*h_f) ) /(-g_f);
                
                if 0 
%%
                    b_0 = log( lhs+3.5 );
                    lhs_1 = lhs + b_0 + .5* b_0^2;
                    b_1 = log( lhs_1+1 );
                    alpha_t = linspace(0,1,1000)';
                    f_t = zeros(size(alpha_t));
                    for k=1:numel(alpha_t); 
                        f_t(k) = fun( x+ alpha_t(k)*reshape(x_step, size_x));
                    end;
                    pred_q = f+alpha_t*g_f+.5*h_f*alpha_t.^2;
                    pred_e = pred_q + (-g_f)*(exp(alpha_t*b_1)-1-b_1*alpha_t-.5*b_1.^2.*alpha_t.^2);
                    plot(alpha_t,[f_t pred_q pred_e])
                end;
                % find trust region radius for which ratio is approx .9
                % (based on approximation function)
                % => find point on which (exp( b* alpha) -1 -b*a -.5*b^2*a^2)
                %    = -.1 * delta_f_predicted
                
                % approximated alpha. Approximation formula obtained with
                % Mathematica (approx_find_initial_trust_region_radius.nb)
                h = log(lhs);
                alpha = .79/(h+1.6*exp(-.42*h));
                pcg_opts.maxR = x_step_norm * alpha;
            elseif quadratic_approx_ratio < 2
                % initial step was good or a bit better than expected
                if quadratic_approx_ratio<1.1 && quadratic_approx_ratio>.9
                    % very good approximation 
                    pcg_opts.maxR = 4 * x_step_norm;
                else
                    pcg_opts.maxR = 1.5 * x_step_norm;
                end;
            else % improvement much larger than expected: quadratic approximation not good
                % set trust region smaller than current step:
                pcg_opts.maxR = .7 * x_step_norm;
            end;
        elseif (quadratic_approx_ratio >= .75) && (x_step_norm>=.7*pcg_opts.maxR)
            % Good step: increase trust region radius
            pcg_opts.maxR = pcg_opts.maxR * options.TR_scale_good_step;
        elseif (quadratic_approx_ratio < -1)
            pcg_opts.maxR = min( pcg_opts.maxR , x_step_norm) * options.TR_scale_verybad_step;
        elseif (quadratic_approx_ratio < .25)
            pcg_opts.maxR = pcg_opts.maxR * options.TR_scale_bad_step;
        end;
    end;
    if f_new < f
        % accept step 
        x = x_new;
        f = f_new;
        G = G_new;
        H_info = H_info_new;
        H_updated =true;
        clear precon_info f_new G_new H_info_new;
        numstepsaccepted = numstepsaccepted +1;
    end;
end;
    

if nargout>=6
    outp.iterations   = iter;
    outp.funcCount    = numFunEvals;
    outp.cgiterations = cgiterations;
    outp.stepstaken   = numstepsaccepted;
    outp.preconditionerPrepare = preconditionerPrepare;
end;