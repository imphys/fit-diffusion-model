function outp = validateDerivativeAndHessian(fun, x, hessmulfun, preconditioner, bandw)
% validateDerivativeAndHessian(fun, x [, hessmulfun [, preconditioner, bandwidth]])
% This function validates the function 'fun', which typically is function you want to pass
% to an optimization routine. 
% fun should compute 3 outputs:
%   [f, g, H] = fun( x ) 
% where 
%   f : function value in x; should be a scalar value.
%   g : df/dx ;  gradient in x
%   H : d2f/dxdx ; hessian in x, or hessian info that can be passed to hessmulfun
%
% This function compares the analytic gradient and hessian to a numerical estimate. 
% The numerical estimates are computed by the 'gradest' and 'jacobianest' function from 
% the 'DERIVESTsuite' by John D'Errico.
% http://www.mathworks.com/matlabcentral/fileexchange/13490-automatic-numerical-differentiation
% 
% NOTES:
%  - Don't use on large problems, since this function is quite slow and the hessian is always computed explicitly. 
%  - In order to fully trust the validation, pick several x (at random) that cover all (special) cases of your function (if it has them). 
%
% [Hy] = hessmulfun(H, y) : 
%         optional function that can multiply one or more vectors with the hessian.
%         Allows you to avoid computing and/or storing the full hessian of the function.
% [R, mulfun/pvec] = preconditioner(H, bandwidth, DM, DG);
%         optional function that computes a preconditioner R for M = DM*H*DM + DG
%         Note: this function only tests DM = eye and DG = zero.
%         if the second output is a function, it should be able to multiply 1 vector (with numel(x) elements) with R
%         (This is for using my hack into the MATLAB optimization toolbox that allows such function.)
%         The condition number of H and H preconditioned with R is compared. The latter should be lower (obviously). 
%
% Created by Dirk Poot, Erasmus MC

% set defaults:
if nargin<3
    hessmulfun = [];
end;
if nargin<4
    preconditioner = [];
end;

% call function:
[f1] = fun(x);
[f2,g2] = fun(x);
[f,g, H] = fun(x);
outp.f_1outp = f1;
outp.f_2outp = f2;
outp.f_3outp = f;
outp.g_2outp = g2;
outp.g_3outp = g;
outp.H = H;
goodresults = true;
if ~isequal(f1,f2,f) || ~isequal(g2,g)
    warning('Function value or gradient not equal if the number of requested outputs changes. This should not happen.'); % (When the difference is just roundoff errors it might occaisionally be acceptable).
    goodresults = false;
end;
if ~isempty( preconditioner )
    if nargin<5
        % fmin_fast preconditioner creation function
        [R, mulfun] = preconditioner( H , x);
    else
        % fminunc type preconditioner creation function
        DM =  speye(numel(x));
        DG = sparse(numel(x),numel(x));
        [R, mulfun] = preconditioner(H, bandw, DM, DG );
    end;
    eR = zeros(numel(x));
    for k=1:numel(x); 
        xk = zeros(numel(x),1);
        xk(k)=1; 
        eR(:,k) = mulfun(xk, R);
    end;
    sqrteR = sqrtm(eR);
    outp.R = eR;
    outp.sqrtmR = sqrteR;
end;
% Using gradest from the 'DERIVESTsuite' by John D'Errico
% http://www.mathworks.com/matlabcentral/fileexchange/13490-automatic-numerical-differentiation
maxnumx_grad = 100;
reducedtest_grad = 50;
maxnumx_hess = 100;
reducedtest_hess = 50;
if numel(x) > maxnumx_grad 
    redproj_grad = randn(numel(x), reducedtest_grad);
    g_red = (g(:)'*redproj_grad)';
    [gest , err] = gradest(@(xr) fun( x + reshape(redproj_grad*xr,size(x))), zeros(reducedtest_grad,1) );
    outp.g_project = redproj_grad;
else
    % full gradient
    g_red = g(:);
    [gest , err] = gradest(fun, x);
end;
gerr = (g_red-gest(:));
outp.gest = gest;
outp.gest_numericalerror = err;
if any( abs(gerr)./err(:) >2  &  abs(gerr)./abs(g_red+gest(:))>1e-5 )
    disp('[ gradient,   numerical_gradient,  difference_in_gradient    absolute_difference_in_gradient/numerical_error]: ')
    fprintf( '   %15f  ,  %15f  , %15f  ,    %15f\n',[g_red gest(:)  gerr  abs(gerr)./err(:)]');
    warning('error in gradient computation seems to be too large: ');
    goodresults = false;
end;
if numel(x)>maxnumx_hess 
    redproj_hess = randn(numel(x), reducedtest_hess);
    outp.H_project = redproj_hess;
end;
if ~isempty(hessmulfun)
    outp.Hinfo = H;
    if numel(x)>maxnumx_hess 
        H = hessmulfun(H, redproj_hess);
    else
        H = hessmulfun(H, eye(numel(x)));
    end;
    Hred = H;
    outp.H = H;
else
    if numel(x)>maxnumx_hess 
        Hred = H*redproj_hess;
    else
        Hred = H;
    end;
end;
% Using jacobianest from the 'DERIVESTsuite' by John D'Errico
% http://www.mathworks.com/matlabcentral/fileexchange/13490-automatic-numerical-differentiation
if numel(x)>maxnumx_hess 
    % evaluate reduced hessian:
    [Hest , errHest] = jacobianest( @(xr) reshape(shuffleoutputs(fun, 2, 2, {x+reshape(redproj_hess*xr,size(x))}),[],1) , zeros(reducedtest_hess,1));
else
    % evaluate full hessian
    [Hest , errHest] = jacobianest( @(x) reshape(shuffleoutputs(fun, 2, 2, {x}),[],1) , x);
end;
outp.Hest = Hest;
outp.Hest_numericalerror = errHest;

if ~isempty( preconditioner )
    if numel(x)>maxnumx_hess && ~isempty(hessmulfun)
        warning('error probably follows; Cannot evaluate preconditioner performance with reduced hessian.');
    end;
    Hresid_a = sqrteR * H * sqrteR; 
    Hresid_n = sqrteR * Hest * sqrteR; 
    outp.condH = cond(H);
    outp.condHest = cond(Hest);
    outp.condRH = cond(Hresid_a);
    outp.condRHest = cond(Hresid_n);
    if outp.condH<outp.condRH
        warning('Condition number of preconditioned hessian is larger than original hessian => Bad preconditioner');
        goodresults = false;
    end;
    if outp.condRH>50
        warning(['Large condition number of preconditioned hessian. Condition number : ' num2str(outp.condRH) '\nIf the eigenvalues are evenly spread (not tested), many (>30) PCG iterations are needed for decent convergence (1%)']);
        goodresults = false;
    end;
end;

Herr = (Hred-Hest);
if any( abs(Herr(:))./(errHest(:)+1e-4*max(errHest(:))) >2 & abs(Herr(:))./abs(Hred(:)+Hest(:)+1e-6*max(Hred(:)))>1e-4)
    warning('error in hessian computation seems to be too large.');
    goodresults = false;
end;
outp.functionCorrect = goodresults;
if goodresults
    disp('The gradient and hessian appear to be computed correctly');
end;