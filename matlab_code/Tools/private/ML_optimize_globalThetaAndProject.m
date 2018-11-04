function [theta , opts] = ML_optimize_globalThetaAndProject(theta, data, opts)
% theta = ML_optimize_global(theta, data, opts)
% Optimizes all of theta in one global optimization.
% The assumed distribution is specified by opts.logPDFfun
%
% Created by Dirk Poot, Erasmus MC, 22-3-2011
if opts.optimizer==3
    opt = fmin_fast;
    opt.maxIter = opts.maxIter;
    opt.abs_tol_G = opts.tolFun;
    opt.InitialTrustRegionRadius = 10;
else
    opt = optimset('fminunc');
    opt = optimset(opt,'gradObj','on','Hessian','on', 'Display','iter', 'MaxIter',opts.maxIter,'TolFun',opts.tolFun,'largescale','on','PrecondBandWidth',inf,'OutputFcn',@progbupd);
end;
opts.explicitHessian = false;
% Never compute hessian explicitly, it might be very large (and if it's small, thats probably for testing).
% if opts.doRegularize
%     if opts.optimizer==3
%     else
%         regularizer.fun = opts.spatialRegularizer{1};
%         regularizer.explicitHessian = opts.explicitHessian;
%         regularizer.hessMulFun = opts.spatialRegularizer{2};
%     end;
% end;

if ~opts.explicitHessian
    % if not explicit hessian, I use the following convention for hessinfo:
    %
    % hessinfo  = the third output of a function to be optimized
    %
    % [Hy] = hessinfo.hessMulFun( hessinfo, Y )  multiplies the (possibly multi column vector) Y with the hessian
    %  => Hy = H * Y
    % [R, fun] = hessinfo.makePreconditioner( hessinfo , upperbandw, DM, DG )
    % prepare R as preconditioner of
    %    M = DM*H*DM + DG
    % so
    %    fun(x, R)   =approx=   inv(M) * x
    if isempty(opts.maxPCGiters )
        maxPCGiter = min(prod(opts.blockSize)*size(theta,1),60+20*ndims(theta));
    else
        maxPCGiter = opts.maxPCGiters;
    end;
    
    hessmulfun = @(hessinfo, Y) fullThetaAndProjectParCostfun_HessMul(hessinfo, Y, opts);
    if opts.optimizer==3
        %         makePreconditioner = @(hessinfo , x) hessinfo.makePreconditioner( hessinfo );
        opt.HessMult = hessmulfun;
        
        opt.Preconditioner = @(Hinfo, x) [];  % creates identity  preconditioner
        opt.Preconditioner_Multiply = @(precon_info, x) x;  % applies identity preconditioner
        opt.pcg_options.kmax = maxPCGiter;
    else
        opt = optimset(opt,'HessMult',hessmulfun, 'MaxPCGIter',maxPCGiter,'PrecondBandWidth',0);
        % set precondbandwidth to less than inf, to display cgiter correctly.
    end;
else
    hessmulfun = [];
end;
if ~isempty(opts.constraints)
    %opts.constraints_A,  opts.constraints_b
    npar = size(theta,1);
    constrFun = @(x) constraintBlockFun(x, opts.constraints , npar );
end;
% % select 'fields' (i.e. additional parameters per voxel)
% if ~isempty(opts.fields)
%     opts.function = make1arg_anonfun( opts.function, opts.fields(:,:) );
% end;
init = cell(1, 1+numel(opts.projectParameters) );
init{1} = theta(:);
for k=1:numel(opts.projectParameters)
    init{k+1} = opts.projectParameters{k}(:);
end;
init = vertcat(init{:});

progressbar('start',[0 maxPCGiter;0 opts.maxIter]);
if opts.optimizer == 2
    opts.doRegularize = false;
    funcs(1).fun = make1arg_anonfun( @fullThetaAndProjectParCostfun, opts );
    funcs(1).hessmul = hessmulfun;
    funcs(2).fun = opts.spatialRegularizer{1};
    funcs(2).hessmul = opts.spatialRegularizer{2};
    [theta_opt, fval] = conjgrad_nonlin(funcs, init);
    hess.LLfun = [];
    hess.LLregularization = [];
else
    optfun = make1arg_anonfun( @fullThetaAndProjectParCostfun, opts );
    % [f,g,h] = optfun( theta) ;HF = @(x) opt.HessMult( h, x); [a,b,c,d,e] = inverse_free_condest( HF, numel(theta), 1, 1000 );
    % imagebrowse(cat(5,reshape(d,size(theta)), reshape(e,size(theta))),[-1 1]*.01)
    if isempty(opts.constraints)
        if opts.optimizer == 3
            [theta_opt, fval, exflag, grad, hess, outp] = fmin_fast(optfun, init, opt);
        else
            [theta_opt, fval, exflag, outp, grad, hess] = fminunc(optfun, init, opt);
        end;
        
    else
        [theta_opt, fval, exflag, outp, lambd, grad, hess] = fmincon(optfun, theta, [],[],[],[],[],[], constrFun ,opt);
    end;
end;

procdpars = numel(theta);
theta = reshape( theta_opt(1:numel(theta)), size(theta));
for k=1:numel(opts.projectParameters)
    lastparuse = procdpars + numel(opts.projectParameters{k});
    
    opts.projectParameters{k} = reshape( theta_opt(procdpars+1:lastparuse), size(opts.projectParameters{k}));
    procdpars = lastparuse;
end;

progressbar('ready');
opts.curpos.LLdata = hess.LLfun;
opts.curpos.LLregularization = hess.LLregularization;
if opts.doComputeDerivativeRegularizationScale
    error('cannot do doComputeDerivativeRegularizationScale');
%     % Compute dThetadLamba = - H^(-1) * dRegularizationdTheta
%     %
%     opts.curpos.dRegularizationdTheta = hess.dRegularizationdTheta;
%     %     opts.curpos.dThetadLamba = pcg( @(Y) fullThetaCostfun_HessMul(hess, Y,opts), hess.dRegularizationdTheta, 1e-4, 100);
%     opts.curpos.dThetadLamba =  - .5 * cgiterLS( @(x) 0, @(x) x, hess.dRegularizationdTheta,[],@(Y) fullThetaCostfun_HessMul(hess, Y,opts), 100);
end;

function [stop] =progbupd(xOutputfcn,optimValues,state)
progressbar( [optimValues.cgiterations ; optimValues.funccount ] );
stop = false;


