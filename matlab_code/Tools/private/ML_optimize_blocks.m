function [theta , opts] = ML_optimize_blocks(theta, data, opts)
% theta = ML_optimize_blocks(theta, data, opts)
% Optimizes all of theta, a block at a time.
% The assumed distribution is specified by opts.logPDFfun
% the used blocks are specified by blockSize, blockOverlap, startBlockPos and computeRange
%
%
% Created by Dirk Poot, Erasmus MC, 22-3-2011
% Adapted by Gwendolyn Van Steenkiste, Vision Lab University of Antwerp,
% 5-12-2013 : extra blocks method that does not need/use the PSF

opts.explicitHessian = (~opts.doRegularize || isempty(opts.spatialRegularizer{2})) && isempty(opts.project) ;
% can't compute explicit hessian if regularizer uses hessmul function.
% can't compute explicit hessian if projecting.
doconstrainedoptimization = ~isempty(opts.constraints);
    
        
if doconstrainedoptimization
    opt = optimset('fmincon');
    opt = optimset(opt,'gradObj','on','Hessian','user-supplied', 'Display','off', 'MaxIter',opts.maxIter,'TolFun',opts.tolFun,'largescale','on','PrecondBandWidth',inf,'GradConstr','on');
else
    if opts.optimizer==3
        opt = fmin_fast;
        opt.maxIter = opts.maxIter;
        opt.TolFun = opts.tolFun;
        opt.InitialTrustRegionRadius = 10;
    elseif opts.optimizer==4
        opt = new_MCMC_sample;
        opt.maxIter = opts.maxIter;
        opts.explicitHessian = true; % really need explicit hessian (for now)
    else
        opt = optimset('fminunc');
        opt = optimset(opt,'gradObj','on','Hessian','on', 'Display','off', 'MaxIter',opts.maxIter,'TolFun',opts.tolFun,'largescale','on','PrecondBandWidth',inf);
    end;
end;

if opts.doRegularize
    regularizer.fun = opts.spatialRegularizer{1};
    regularizer.explicitHessian = opts.explicitHessian;
    regularizer.hessMulFun = opts.spatialRegularizer{2};
    if numel(opts.spatialRegularizer)>=3
        regularizer_prepare = opts.spatialRegularizer{3};
    else
        regularizer_prepare = [];
    end;
end;

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
        maxPCGiter = min(prod(opts.blockSize),8+3*ndims(theta))*size(theta,1);
    else
        maxPCGiter = opts.maxPCGiters;
    end;
    hessmulfun = @(hessinfo, Y) hessinfo.hessMulFun(hessinfo, full(Y)); % The input Y might be sparse if optimization includes just 1 parameter.
    if opts.optimizer==3  
       
        makePreconditioner = @(hessinfo , x) hessinfo.makePreconditioner(hessinfo );
        opt.HessMult = hessmulfun;
        opt.Preconditioner =  makePreconditioner;
        opt.Preconditioner_Multiply = compactedHessianMakePreconditioner();
        opt.pcg_options.kmax = maxPCGiter;
        

    elseif opts.optimizer==4
        opt.HessMult = hessmulfun;
        makePreconditioner = @(hessinfo , x) hessinfo.makePreconditioner( hessinfo );
    else
        makePreconditioner = @(hessinfo , upperbandw, DM, DG) hessinfo.makePreconditioner( hessinfo , upperbandw, DM, DG );
        PrecondBandWidth = prod(opts.blockSize)*size(theta,1);
        opt = optimset(opt,'HessMult',hessmulfun, 'Preconditioner', makePreconditioner, 'MaxPCGIter',maxPCGiter,'PrecondBandWidth',PrecondBandWidth);
    end;
    % set precondbandwidth to less than inf, to display cgiter correctly.
else
    hessmulfun = [];
    PrecondBandWidth = [];
    makePreconditioner = [];
end;
if doconstrainedoptimization
    %opts.constraints_A,  opts.constraints_b
    npar = size(theta,1);
    if isempty( opts.constraints.Aineq ) && isempty( opts.constraints.Aeq ) && isempty( opts.constraints.nonlcon )
        constrFun = [];    
    else
        constrFun = @(x) constraintMulBlockFun(x, opts.constraints , npar );
%         opt = optimset(opt,'Algorithm','trust-region-reflective','HessMult',@(H,x) H*x);
        opt = optimset(opt,'Algorithm','active-set','HessMult',@(H,x) H*x,'Hessian',[]);
        % Unfortunately non of the algorithms of fmincon can be used with the hessian that we know:
        % active-set and sqp  : cannot handle hessians (does not accept the
        %                hessian of the cost function, nor a hessian multiplicaiton
        %                function )
        % trust-region-reflective : cannot handle nonlinear constraints.
        % interior-point : Could in principle do it. However only very inefficiently and cumbersome
        %                  It only accepts a 'hessFcn' which should return
        %                  an explicit hessian when provided x and lambda.
        %                  So this function should be recreated for each
        %                  block (that's just programming annoyance) and
        %                  return the third output of the cost function
        %                  (assuming the hessian is explicitly calculated).
        %                  This last thing means computing the value and gradient twice.
        %                  Alternatively a hessmult function can be
        %                  provided, however it cannot be given hessian
        %                  information (so it should recompute the hessian
        %                  at each iteration !!!! How stupid can it be)
        %                  and has a different interface (ARGH!!)
        % 
    end;
end;
spatialNoiseLevel = false;
if isempty(opts.project) || ~opts.MakePSFApproxOfProj
    if opts.numPDFoptpar==0
        NLsize = size(opts.noiseLevel);
        spatialNoiseLevel = ~all( NLsize(2:end)==1 );
        if ~spatialNoiseLevel
            noiseLevel = opts.noiseLevel;
            LLfun = @(data, A) opts.logPDFfun(data, A, noiseLevel, [false true false]);  % derivative w.r.t. A
        end;
    else
        LLfun = opts.logPDFfun; % exclusively passed to voxelLLfun_sigma_m
    end;
else
    dLLdf_full = cell(1,numel(opts.projectGrad));
    for k= 1 : numel(opts.projectGrad)
        dLLdf_full{k} = project1( opts.projectGrad{k}.logPDFgrad , opts.projectGrad{k}.projectGrad , opts.project(k,:), 3 );
        dLLdf_full{k} = reshape(dLLdf_full{k}, [1, opts.spatialSize] );
    end;
    dLLdf_full = vertcat(dLLdf_full{:});
    hLLifo = opts.projectScaledPartPSF.prepareGlobal( dLLdf_full, opts);
    predict_in = data;
    theta_in = theta;
    dLLdf_full_in = dLLdf_full;
end;

fastExtraInitialValues = opts.initialValueSpecifierVect(4) && isempty(opts.project) && ~opts.doRegularize && isempty(opts.fields) && (opts.numPDFoptpar==0);
% To use fastExtraInitialValues we need
%  - specific (fixed) test vectors (extraInitialValues) that are the same for all voxels
%  - no projection, since that mixes the signals of voxels
%  - no regularization, since that couples the optimization of the voxels
%  - no fields since a field (almost by definition) is not equal for all
%    voxels (so it would be required to construct the dot product matrices
%    that are constructed right below for each voxel separately).
%  - no PDF optimization parameters since we (implicitly) assume that all
%    parameters change the predicted signals and we implicitly assume
%    gaussian noise.
if fastExtraInitialValues
    extraInitialValuesMat = [opts.extraInitialValues{:}];
    if isempty( opts.linearParameters )
        fastInitialPred = opts.function( extraInitialValuesMat );
        if isempty(opts.parameterPrior.fun)
            extraInitialValuePriorLL = zeros(1,size(extraInitialValuesMat,2));
        else
            extraInitialValuePriorLL = opts.parameterPrior.fun( extraInitialValuesMat );
        end;
    else
        % (implicitly) compute least squares estimate of scaling parameters
        % => assumes gaussian PDF.
        if ~isempty(opts.parameterPrior.fun)
            error('cannot handle parameter prior when scaling parameters and initial values are present.');
            % we would need to know how parameterprior factors over the
            % parameters so we can split its contribution as well.
        end;
        % What is computed here?
        % assume fun = Xi * theta_scaling
        % with Xi(:, k) = fun( extraInitialValues{i}_with_all_scaling_parameters_0_except_k )
        % This requires fun to be of the form 
        %  fun = theta_scaling(1) * fun1(theta_red) + theta_scaling(2)*fun2(theta_red) ....
        % where theta_scaling = theta( opts.linearParameters , ..)
        % and theta_red = extraInitialValues{i}( ~opts.linearParameters , ..)    (i.e. theta with the scaling parameters removed)
        %
        % In this case the best fitting initial value is given by
        %  i = arg min_i   min_{theta_scaling} || Xi * theta_scaling - Y||^2
        % where Y is the data in the 'current' voxel
        %  The internal minimization can be solved in closed form.
        %  Let [Q,Ri] = qr( Xi )
        %  and Pi = Xi / Ri
        %  then Pi'*Pi = eye( size(Pi,2)) 
        %  and theta_scaling = Ri \ Pi'*Y
        %  The residu norm is then given by 
        %    || Pi * Pi'*Y - Y ||^2 = Y' * ( Pi*Pi' - eye ) * ( Pi * Pi' -eye) * Y
        %    = Y' * (eye - Pi * Pi') * Y
        %  since Y is fixed, arg min residue is searching for the largest
        %  || Pi'*Y ||^2
        % Directly below we compute Ri and Pi
        % For each block we compute Pi'*Y and search for the (one or more)
        % largest norms of that.
        % Only afterwards theta_scaling is computed explicitly.
        fastInitialPred = cell(1,numel(  opts.linearParameters ));
        for k = 1 : numel(  opts.linearParameters )
            tmp = extraInitialValuesMat; tmp(opts.linearParameters,:)=0;
            tmp( opts.linearParameters(k) ,:) = 1;
            fastInitialPred{k} = opts.function( tmp ); clear tmp;
        end;
        fastInitialPred = cat(3,fastInitialPred{:});
        fastInitialPred_chol = cell( 1, size(fastInitialPred ,2) );
        for k=1:size(fastInitialPred ,2)
            tmp = reshape( fastInitialPred(:,k,:), size(fastInitialPred,1), size(fastInitialPred,3));
            [dummyQ, fastInitialPred_chol{k}] = qr(tmp,0);
            fastInitialPred(:,k,:) = tmp/fastInitialPred_chol{k}; % automatically reshaped on store.
        end;
        clear tmp dummyQ;
    end;
end;

szData = size(data);
ndim = numel(szData)-1;
nblocksDim = zeros(1,ndim);

st = cell(1,ndim);
ed = cell(1,ndim);
allsel = cell(1,ndim);
allselLrg = cell(1,ndim);
for dim = 1:ndim
    if opts.periodicDomain(dim)
        endcmptdim = opts.computeRange(2,dim) + opts.blockOverlap(dim);
    else
        endcmptdim = opts.computeRange(2,dim);
    end;
    st{dim} = opts.startBlockPos(dim):opts.blockSize(dim)-opts.blockOverlap(dim):endcmptdim;
    ed{dim} = min(st{dim}+opts.blockSize(dim)-1, endcmptdim);
    if ed{dim}(1)<opts.computeRange(1,dim)
        % if first block ends before start of computeRange, remove all those blocks:
        selectblksdim = ed{dim}>=opts.computeRange(1,dim);
        ed{dim} = ed{dim}(selectblksdim);
        st{dim} = st{dim}(selectblksdim);
    end;
    st{dim} = max(opts.computeRange(1,dim),st{dim});
    if opts.blocksOptimizeDirection(dim)<0
        % flip ordering of blocks in dim
        st{dim} = st{dim}(end:-1:1);
        ed{dim} = ed{dim}(end:-1:1);
    end;
    nblocksDim(dim) = numel(st{dim});
    
    allsel{dim} = cell(1,nblocksDim(dim));
    allselLrg{dim} = cell(1,nblocksDim(dim));
    for kdim=1:nblocksDim(dim)
        allsel{dim}{kdim} = st{dim}(kdim):ed{dim}(kdim);
        if opts.periodicDomain(dim)
            allsel{dim}{kdim} = mod( allsel{dim}{kdim} -1, szData(dim+1))+1;
            allselLrg{dim}{kdim} = mod( st{dim}(kdim)-opts.spatialRegularizeBorder(dim)-1: ed{dim}(kdim)+opts.spatialRegularizeBorder(dim)-1, szData(dim+1))+1;
        else
            allselLrg{dim}{kdim} = max(1,st{dim}(kdim)-opts.spatialRegularizeBorder(dim)): min(ed{dim}(kdim)+opts.spatialRegularizeBorder(dim), szData(dim+1));
        end;
    end;
end;
nblocks = prod(nblocksDim);
sel = cell(1,ndim+1); sel{1} = ':';
selLrg = sel;
indexmap = zeros([szData(2:end) 1]);
linindex = reshape(1:prod(szData(2:end)),[szData(2:end) 1]);
linindexDimstep = cumprod([1 szData(2:end-1)]);
nerrors = 0;
if isempty(opts.mask)
    nprocesstraces = prod(szData(2:end) + (nblocksDim-1).*opts.blockOverlap);
else
    if all(opts.blockOverlap==0)
        nprocesstraces = nnz(opts.mask);
    else
        fullsel = cell(1,ndim);
        for dim = 1: ndim
            %             s = cell(1,numel(st{dim}));
            %             for k=1:numel(s)
            %                 s{k} = st{dim}(k):ed{dim}(k);
            %             end;
            %             fullsel{dim} = [s{:}];
            fullsel{dim} = [allsel{dim}{:}];
        end;
        nprocesstraces = nnz(opts.mask(linindex(fullsel{:})));
    end;
end;
if opts.save_initializationorder
    opts.saved_initializationorders = cell(1,nblocks);
end;
nprocessedtraces =[0 0 0];
repeatPSFscaleIterations = true;
showprogress = (nblocks>1) && (opts.progressbarUpdateStep>0) && (nprocesstraces>opts.progressbarUpdateStep);
while repeatPSFscaleIterations % Only repeat when we need to rescale the PSF, since we can't find an acceptable point; probably due to non linearity of fun.
    repeatPSFscaleIterations = false;
    if showprogress
        progressbar('start',nprocesstraces,[],'EstTimeLeft','on','minTimeInterval',1);
        nextshow = sum(nprocessedtraces) + opts.progressbarUpdateStep;
    end;
    for blockIdx=0:nblocks-1
        krem = blockIdx;
        % Locate block index and create linear index of the (within mask) voxels:
        szSelLrg = ones(1,max(ndim,2));
        blkidxdim = zeros(1,ndim);
        for dim = opts.blocksOptimizeDimOrder
            kdim = mod(krem, nblocksDim(dim))+1;
            blkidxdim(dim) = kdim;
            krem = floor(krem/ nblocksDim(dim));
            sel(dim+1) = allsel{dim}(kdim); %st{dim}(kdim):ed{dim}(kdim);
            selLrg(dim+1) = allselLrg{dim}(kdim);
            szSelLrg(dim)  = numel(selLrg{dim+1});
        end;
        sellin = linindex(sel{2:end});
        if ~isempty(opts.mask)
            sellin = sellin( opts.mask( sellin ) );
            if isempty(sellin)
                continue;
            end;
        end;
        
        % select data:
        datasel = double(data(:,sellin));
        
        % get initial values.
        if opts.initialValueSpecifierVect(1) % provided initial value
            thetasel = {theta(:,sellin)};
        else
            thetasel = {}; % Why would you not want to use the provided value??? But I allow you to.
        end;
        if opts.initialValueSpecifierVect(2) % previously computed neighboring blocks.
            for dim = find(blkidxdim>1)
                % if not first block in dimension dim, shift sellin to previous block.
                sellin_adj = sellin + linindexDimstep(dim)*( st{dim}(blkidxdim(dim)+[0 -1]) *[-1;1]);
                thetasel{end+1} = theta(:,sellin_adj);
            end;
        end;
        if opts.initialValueSpecifierVect(3) % not yet computed neighboring blocks.
            for dim = find(blkidxdim<nblocksDim)
                % if not last block in dimension dim, shift sellin to previous block.
                sellin_adj = sellin + linindexDimstep(dim)*( st{dim}(blkidxdim(dim)+[0 1]) *[-1;1]);
                thetasel{end+1} = theta(:,sellin_adj);
            end;
        end;
        if opts.initialValueSpecifierVect(4) && ~fastExtraInitialValues% Extra global initial values
            thetasel = [thetasel opts.extraInitialValues];
            if numel(sellin)>1
                for idx = numel(thetasel)-(0:numel(opts.extraInitialValues)-1)
                    thetasel{idx} = repmat(thetasel{idx}, 1, numel(sellin));
                end;
            end;
        end;
        
        % select 'fields' (i.e. additional parameters per voxel)
        if isempty(opts.fields)
            fun = opts.function;
        else
            fieldsel = opts.fields(:,sellin);
            %         fun = @(par) opts.function(par, fieldsel);
            fun = make1arg_anonfun( opts.function, fieldsel );
        end;
        if spatialNoiseLevel
            noiseLevel = opts.noiseLevel(:,sellin);
            LLfun = @(data, A) opts.logPDFfun(data, A, noiseLevel, [false true false]);  % derivative w.r.t. A
        end;
        
        % Build final function to optimize in this iteration:
        
        if opts.numPDFoptpar==0
            if ~isempty(opts.project) && opts.MakePSFApproxOfProj
                
                % projecting, so dont use LLfun explicitly, but rather its derivative and hessian:
                dLLsel   = dLLdf_full(:,sellin);
                [selLrgPSF , hLLifoBlk]= opts.projectScaledPartPSF.prepareRegion( sellin(:).', hLLifo );
                selLrgPSF = filterStillToProcess( selLrgPSF, st, ed, blkidxdim, opts);
                selLrgPSF = ((selLrgPSF-1)*cumprod([1 opts.spatialSize(1:end-1)])')'+1;
                [hLLifoBlk , psfdiag]  = opts.projectScaledPartPSF.prepareBlock( sellin(:).', selLrgPSF, hLLifo, hLLifoBlk);
                hLLselfun = make1arg_anonfun( opts.projectScaledPartPSF.mulsel, hLLifoBlk );
                voxLLfun = make1arg_anonfun( @voxelLLfun_proj_m, fun, datasel, dLLsel, hLLselfun, psfdiag, opts.maxfunArgsOut, opts.explicitHessian , opts.parameterPrior, opts.project_PSFscaleBlock, opts.funHessian_I, opts.funHessian_J );
                
            elseif ~opts.MakePSFApproxOfProj
                [opts.project,opts.projectParameters] = opts.MakeProjBlockFun(theta,fun,sellin); % projplan en transformation zijn al meegegeven
                voxLLfun = @(p) voxelLLfun_proj(      p, fun, LLfun, datasel , opts.maxfunArgsOut, opts.explicitHessian , opts.parameterPrior, opts.funHessian_I, opts.funHessian_J,opts);
%                 voxLLfun = @(p) voxelLLfun_m(      p, fun, LLfun, datasel , opts.maxfunArgsOut, opts.explicitHessian , opts.parameterPrior,                    opts.funHessian_I, opts.funHessian_J);

            else
                voxLLfun = @(p) voxelLLfun_m(      p, fun, LLfun, datasel , opts.maxfunArgsOut, opts.explicitHessian , opts.parameterPrior,                    opts.funHessian_I, opts.funHessian_J);
            end;
        else
            voxLLfun = @(p) voxelLLfun_sigma_m(    p, fun, LLfun, datasel , opts.maxfunArgsOut, opts.explicitHessian , opts.parameterPrior, opts.numPDFoptpar, opts.funHessian_I, opts.funHessian_J);
        end;
        
        if opts.doRegularize
            indexmap(selLrg{2:end}) = reshape(1:prod(szSelLrg),szSelLrg);
            subindexsel = indexmap(sellin);
            thetaselLrg = theta(selLrg{:});
            if isempty(regularizer_prepare)
                optfun = @(p) regularizedLLFun(p, voxLLfun, regularizer, thetaselLrg, subindexsel);
            else
                [regulfun , regulHessmul] = regularizer_prepare( opts.spatialRegularizer, regularizer.explicitHessian, thetaselLrg, subindexsel );
                if isempty(regulfun)
                    optfun = @(p) regularizedLLFun(p, voxLLfun, regularizer, thetaselLrg, subindexsel);
                else
                    regularizer.fun =regulfun;
                    regularizer.hessMulFun = regulHessmul;
                    optfun = @(p) regularizedLLFun(p, voxLLfun, regularizer);
                    %                 optfun = makeregularizedLLFun(voxLLfun, regularizer);
                end;
            end;
        else
            optfun = voxLLfun;
        end;
        if (islogical(opts.validateFullFunction) &&  opts.validateFullFunction) || (~islogical(opts.validateFullFunction) && ismember(blockIdx, opts.validateFullFunction))
            if exist('PrecondBandWidth','var')
                dummy = validateDerivativeAndHessian(optfun, thetasel{1}, hessmulfun, makePreconditioner, PrecondBandWidth);
            else
                dummy = validateDerivativeAndHessian(optfun, thetasel{1}, hessmulfun, makePreconditioner);
            end;
        end;
        
        % Reduce number of initial values with which the optimization is started:
        if numel(thetasel)>opts.onlyKeepnBestStartPoints || fastExtraInitialValues || opts.maxIter==0
            % general case, can't really do much else than evaluate full
            % function:
            fval = zeros(1,numel(thetasel));
            for initidx = 1:numel(thetasel)
                fval(initidx) = optfun(thetasel{initidx});
            end;
            if fastExtraInitialValues
                if isempty( opts.linearParameters )
                    % use initially predicted signal intensities and only evaluate
                    % data likelihood
                    if numel(sellin)~=1
                        [felI,felJ]=  ndgrid(1:numel(sellin), 1:size(fastInitialPred,2 ));
                        ll = -reshape( sum( LLfun( datasel(:,felI), fastInitialPred(:,felJ) ) ,1 ) + extraInitialValuePriorLL(felJ(:)) , numel(sellin), size(fastInitialPred,2 ));
                        [ lls, fastInitialReord ] = sort(ll,2);
                        ll = sum(lls,1);
                    else
                        ll = -sum( LLfun( datasel(:,ones(1,size(fastInitialPred,2))), fastInitialPred) ,1 ) + extraInitialValuePriorLL;
                    end;
                    fval =[fval ll];
                else
                    % Find the intial values that (after optimizing the
                    % scaling parameters) leave the lowest residue:
                    % See above at construction of fastInitialPred for more
                    % explanation.
                    dotprod = reshape( datasel' * fastInitialPred(:,:) , [size(datasel,2) size(fastInitialPred,2), size( fastInitialPred,3)]);
                    sumnorm2 = -sum( dotprod.^2 ,3) ;
                    [ lls, fastInitialReord ] = sort( sumnorm2 , 2); % sort over initial vectors first indices have largest norm
                    nkeep = min( opts.onlyKeepnBestStartPoints,size(fastInitialReord ,2) );
                    % Now actually compute the optimal scaling parameters
                    % and construct the initial parameter vectors:
                    thetaseladd = cell(1,nkeep);
                    ll = zeros(1,nkeep);
                    for k = 1 : nkeep
                        newparam = zeros( size(theta,1), numel(sellin) );
                        for k2 = 1 : size(datasel,2)
                            scaleparamvalue =  fastInitialPred_chol{ fastInitialReord(k2,k) } \ permute( dotprod(k2, fastInitialReord(k2,k), : ) ,[3 1 2]);
                            newparam(:,k2) = opts.extraInitialValues{ fastInitialReord(k2,k) };
                            newparam(opts.linearParameters ,k2) = scaleparamvalue;
                        end;
                        thetaseladd{k} = newparam;
                        ll(k) = optfun( thetaseladd{k} );
                    end;
                    fval = [fval ll];
                    thetasel = [thetasel thetaseladd];
                end;
            end;
            [dummy, bestinitidxs ] = sort(fval);
            if opts.save_initializationorder
                opts.saved_initializationorders{blockIdx+1} = bestinitidxs;
            end;
            bestinitidxs = bestinitidxs(1:min(opts.onlyKeepnBestStartPoints,numel(bestinitidxs)));
            if fastExtraInitialValues && isempty( opts.linearParameters )
                bestinitidxlrg = bestinitidxs(bestinitidxs>numel(thetasel) )- numel(thetasel);
                if numel(sellin)~=1 && ~isempty(bestinitidxlrg)
                    bestinitidxlrgOrig = fastInitialReord(:,bestinitidxlrg);
                    thetaseladd = cell(1,numel(bestinitidxlrg));
                    for sel_vox_id = 1 : numel(bestinitidxlrg)
                        thetaseladd{sel_vox_id} = [opts.extraInitialValues{bestinitidxlrgOrig(:,sel_vox_id)}];
                    end;
                else
                    thetaseladd = opts.extraInitialValues(bestinitidxlrg);
                end;
                thetasel = [thetasel(bestinitidxs(bestinitidxs<=numel(thetasel))) thetaseladd];
            else
                thetasel = thetasel(bestinitidxs);
                fval = fval(bestinitidxs);
            end;
            
        end;
        
        % Call optimizer for all initial values.
        if opts.maxIter==0
            thetaselopt = thetasel;
            % fval = fval ; % fval is already set. 
            nprocessedtraces(3) = nprocessedtraces(3) + numel(sellin);
        else
        fval = zeros(1,numel(thetasel));
        thetaselopt = cell(size(thetasel));
        for initidx = 1:numel(thetasel)
            try
                if doconstrainedoptimization
                    if ~isempty(opts.constraints.lb) || ~isempty(opts.constraints.ub)
                        lb = opts.constraints.lb*ones(1, size(thetasel{initidx},2));
                        ub = opts.constraints.ub*ones(1, size(thetasel{initidx},2));
                    else
                        lb = [];ub=[];
                    end;
                    [thetaselopt{initidx}, fval(initidx), exflag(initidx), outp] = fmincon(optfun, thetasel{initidx}, [],[],[],[],lb,ub, constrFun ,opt);
                else
                    % unconstrained optimization methods:
                    if opts.optimizer==3
                        [thetaselopt{initidx},fval(initidx), exflag, G, H_info, outp] = fmin_fast(optfun, thetasel{initidx}, opt);
                    elseif opts.optimizer==4
                        [thetaselopt{initidx},fval(initidx)] = new_MCMC_sample(optfun, thetasel{initidx}, opt);
                    else
                        [thetaselopt{initidx}, fval(initidx), exflag(initidx), outp] = fminunc(optfun, thetasel{initidx}, opt);
                    end;
                end;
                %         if isequal(thetaselopt{initidx},thetasel{initidx})
                %             error('Optimization returned initial position, probably incorrect');
                %         end;
                nprocessedtraces(3) = nprocessedtraces(3) + numel(sellin);
            catch ME  % REMOVE 'ME' for old MATLAB (and adjust below as indicated)
                fval(initidx) = inf;
                nerrors = nerrors +1;
                if nerrors <= opts.maxErrorNotifications
                    str = [];
                    for dim = 1:ndim
                        if numel(sel{dim+1})==1
                            str = [str ', ' num2str(sel{dim+1})];
                        else
                            str = [str ', ' num2str(sel{dim+1}(1)) ':' num2str(sel{dim+1}(end))];
                        end;
                    end;
                    if nerrors >= opts.maxErrorNotifications;
                        extra = '\nSuppressing further error notifications.';
                    else
                        extra = '';
                    end;
                    % FOR OLD MATLAB replace:
                    warning('Fit_MRI:FailedBlock',['Error during optimization of block ' num2str(blockIdx+1) '/' num2str(nblocks) ';  theta(:' str ') \n   identifier : ''' ME.identifier '''\n   message    : ' ME.message '\n   file       : ' ME.stack(1).name ' at line ' num2str(ME.stack(1).line) extra]);
                    % by:
                    %   warning('Fit_MRI:FailedBlock',['Error during optimization of block ' num2str(blockIdx+1) '/' num2str(nblocks) ';  theta(:' str ')' ]);
                    % (to avoid using 'ME' at the expense of less detailed information)
                end;
                nprocessedtraces(2) = nprocessedtraces(2) + numel(sellin);
            end;
        end;
        end;
        % Select the best value found and update theta with it.
        [bestfval, bestoptidx] = min(fval);
        if isfinite(bestfval) || bestfval==-inf
            % only update when best function value is finite (except for -inf)
            if ~isempty(opts.project) && opts.MakePSFApproxOfProj
                newfval = fun( thetaselopt{bestoptidx} );
                deltafval =  newfval - data(:,sellin);
                % Update dLLdf_full with update in function value (=data)
                % since
                % dLL / df_i =approx= dLL/df_i + d2LL/ df_i df_j * ( new_f_j - old_f_j )
                dLLdf_full(:, selLrgPSF)  = dLLdf_full(:, selLrgPSF) + opts.projectScaledPartPSF.mulfull( opts.projectScaleDeltaPSFupdate * deltafval , hLLifoBlk);
                data(:,sellin) = newfval;
            elseif ~isempty(opts.project) && ~opts.MakePSFApproxOfProj
                newfval = fun( thetaselopt{bestoptidx} );
                data(:,sellin) = newfval; % update data, required for overlapping blocks and evaluate new total likelihood afterwards.
            end;
            % Store update
            theta(:,sellin) = thetaselopt{bestoptidx};
        end;
        if showprogress && (nextshow <= sum(nprocessedtraces) )
            progressbar(nprocessedtraces);
            nextshow = sum(nprocessedtraces) + opts.progressbarUpdateStep;
        end;
    end;
    if showprogress
        progressbar('ready');
    end;
    if ~isempty(opts.project) && opts.MakePSFApproxOfProj
        % Test if PSF (and PDF hessian) was modelled accurately enough for the update:
        LLpredict_in = projectAllAndLL( predict_in, opts );
        delta_predict = data-predict_in;
        gradLL0 = dLLdf_full_in(:)'*delta_predict(:);
        [LLnew] = projectAllAndLL( data, opts );
        
        deltaLL = sum(LLnew - LLpredict_in ); % first subtract, then sum for (slightly) improved numerical accuracy.
        % fit polynome p(x) through  (0,0) , (1,deltaLL), with d p / d x |_{x=0} = gradLL0
        % => 2nd order polynome;  a x^2 + b x + c
        %      => c =0 , b = gradLL0, a = deltaLL - gradLL0
        % Bad step if a<0 or extremum not close to x=1.
        if deltaLL <= gradLL0
            % negative curvature
            disp('Negative curvature detected => PSF not accurate. Step taken, since likelihood improved more than expected. However, since this was unexpected, the step probably is rather poor.');
        else
            % compare gradient at step. Also take gradient of regularization into account!
            delta_theta = theta-theta_in;
            if opts.doRegularize
                %             if isempty(opts.previousRegularizationGradient)
                %                 [opts.previousRegularizationValue, opts.previousRegularizationGradient] = opts.spatialRegularizer{1}(theta_in);
                %             end;
                [LLreg0, gLLreg0] = opts.spatialRegularizer{1}(theta_in);
                %             gR0 = opts.previousRegularizationGradient(:)'*delta_theta(:);
                gR0 = gLLreg0(:)'*delta_theta(:);
                opts.previousRegularizationGradient = [];
                %             [opts.previousRegularizationValue, opts.previousRegularizationGradient ] = opts.spatialRegularizer{1}(theta);
                [LLreg1, gLLreg1] = opts.spatialRegularizer{1}(theta);
                %             gR1 = opts.previousRegularizationGradient(:)'*delta_theta(:);
                gR1 = gLLreg1(:)'*delta_theta(:);
            else
                LLreg0 = 0;
                gR0 = 0;
                gR1 = 0;
            end;
            % gradient of p(x) |_x==1  = 2 a + b
            gAt_x_1 = gR1 + 2 * deltaLL - gradLL0; % gradient at step end, assuming quadratic polynomial for 'PSF', no assumptions about regularization function.
            gAt_x_0 = gR0 + gradLL0;
            if abs(gAt_x_1) < abs(gAt_x_0) * .2
                % accept step, gradient in step direction is expected to be close enough to 0
            else
                % don't accept step. (Making this choise is compuationally expensive)
                % Do a line search to locate the best step fraction (using all knowledge available):
                linesearchcritfun = @(alpha) lineSearchThetaCostfun(theta_in, delta_theta, alpha, opts);
                %             [alpha , fn ] = linesearch( linesearchcritfun, 0, sum(opts.LLdata) + opts.previousRegularizationValue, gAt_x_0, 1, [2 2 .99 .1 5]);
                [f0,g0]=linesearchcritfun(0); % can largely be optimized away.
                [f1,g1]=linesearchcritfun(1); % can largely be optimized away.
                %             [alpha  ] = linesearch_1( linesearchcritfun, sum(LLpredict_in) + opts.previousRegularizationValue, gAt_x_0, f1, g1, 5);
                [alpha  ] = linesearch_1( linesearchcritfun, f0, g0, f1, g1, 5)
                
                theta = theta_in + alpha * delta_theta;
                
                % if function value is not sufficiently decreased, start optimization again with higher PSF scaling.
                %  (This essentially waists all computations performed for the current step)
                
            end;
        end;
        
        if 0
            adj = linspace(0,2,30);
            LLtst = zeros(6000,30);
            LLreg = zeros(1,30);
            for k=1:30;
                dtst = predict_in + delta_predict*adj(k);
                LLtst(:,k) = opts.logPDFfun(opts.data_in{1}, opts.project{1}( dtst ), 1);
                LLreg(k) = opts.spatialRegularizer{1}( dtst );
            end;
            sLL = -sum(LLtst);
            plot(adj,[(sLL-min(sLL))' LLreg'-min(LLreg)]*[1 0 1;0 1 1])
        end;
    end;
end; % repeatPSFscaleIterations
if nerrors > 1 %= opts.maxErrorNotifications
    warning('Fit_MRI:FailedBlocks', ['In total, optimization of ' num2str(nerrors) '/' num2str(nblocks) ' blocks failed with an error. The initial value is returned for these blocks.']);
end;
opts.encounteredErrors = opts.encounteredErrors + nerrors;




function points = filterStillToProcess( points, st, ed, blkidxdim, opts)
% Preserves the points that are still used in subsequent block optimizations.

% make sure points are inside reconstruction volume:
keepmask = all( points >= 1, 2) & all( bsxfun(@minus, points, opts.spatialSize)<=0 ,2);
points = points( keepmask , :);
% make sure points are inside mask:
steps = cumprod([1 opts.spatialSize(1:end-1)]);
if isempty(opts.mask)
    keepmask = true( size(points,1), 1);
else
    keepmask = opts.mask( (points-1)*steps' +1 );
end;
% check if we still need to process that point:
indetermined = find(keepmask);
for dim = opts.blocksOptimizeDimOrder(end:-1:1) % reverse order, since slowest dimension should make selection first.
    if opts.blocksOptimizeDirection(dim)<0
        % downwards direction
        if blkidxdim(dim)<numel(ed{dim})
            include = points(indetermined,dim)<=ed{dim}(blkidxdim(dim)+1); % include points that should be processed by next block
        else
            include = false(size(indetermined)); % there is no next block, so shouldnt include anything;
        end;
        exclude = points(indetermined,dim)>ed{dim}(blkidxdim(dim)); % exclude points that are higher than end of current block
    else
        % upward direction
        if blkidxdim(dim)<numel(st{dim})
            include = points(indetermined,dim)>=st{dim}(blkidxdim(dim)+1); % include points that should be processed by next block
        else
            include = false(size(indetermined));
        end;
        exclude = points(indetermined,dim)<st{dim}(blkidxdim(dim));
    end;
    keepmask( indetermined( exclude ) ) = false;
    indetermined = indetermined( ~(exclude | include) );
end;
keepmask( indetermined ) = false; % only points inside current block can be indetermined here, and they are not processed in future optimizations.
points = points( keepmask , :);

% function [optfun] = makeregularizedLLFun( voxLLfun, regularizer)
% optfun = @(p) regularizedLLFun(p, voxLLfun, regularizer);