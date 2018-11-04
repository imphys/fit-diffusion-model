function [theta, residue, opts] = fit_MRI( fun, data_in, theta, varargin)
% [par, residue, finalOptions]  =  fit_MRI( fun, data , par0 [,option, value , ...] )
% This function fits a parametric model to several measurements (MR images).
%
% When is this function applicable? :
%  - You have (series of) independent data (or images).
%     e.g.: images, potentially with different orientations& resolutions
%  - You know the statistical distribution of the data
%     e.g.: 'rician', for magnitude MR images
%  - You want to compute parametermaps (with arbitrary number of dimensions) from this data
%     e.g.: The density and T1 value in each voxel of a volume.
%  - Optionally: These parameter maps can be regularized.
%     e.g.: with Total Variation.
%  - You have a function that, given parameters in one voxel, can 
%    predict the 'image intensity' for a set of images in this voxel.
%     e.g.: MLE_T1_d
%  - Optionally: If the data is not a set of images that co-inside with the map-volume
%                you can provide a (approximately linear) function that can 'project'
%                each predicted image to the corresponding data. Extra transform specifying parameters
%                can be estimated as well, 
%     e.g.: affineTransform; where the affine transform parameters can be fine-tuned.
%
% The parametric model is given by fun. There are many options that allow you to specify initial 
% values, whether the rician nature of the data should be taken into account, etc. 
%
% Some basic checks on the provided options are performed, but they are not thoughroughly 
% validated (as this would be time consuming and difficult to do with the huge 
% flexibility provided by this function). Thus, when invalid values are provided 
% errors in unsuspected places might occur. 
% 
% INPUTS:
% fun  : Function that predicts MR intensities or the string 'eye' if you don't want to use a model.
%        fun should be a function handle with one input, the parameters that should be optimized.
%        (See option 'fields' below for a special case when fun should accept two inputs.)
%        Two or three outputs should be returned (if requested): 
%           1) a vector containing the predicted magnitude of the MR images, 
%           2) a matrix containing the derivatives of these intensities w.r.t. the parameter vector 
%           3) the hessian, if the option 'hessian' is 'on' 
%        Inputs and outputs of fun are:
%         [A, dAdpar, d2Adparidparj] = fun( par(:, index_range) )
%          A      = nMRI x numel(index_range) predicted magnitudes
%          dAdpar = nMRI x numel(index_range) x size(par,1) 
%                           derivative of each A with respect to (w.r.t.) each parameter.
%                           ( = Jacobian of function for each column of theta)
%          d2AdparIdparJ = nMRI x numel(index_range) x (numel(I)) 
%                           second derivative (hessian) of each A w.r.t. each parameter combination.
%                           d2 A / d par(I(k)) d par(J(k)) 
%                           with [I,J] = find(triu(ones(size(par,1))));
%        When the option 'fields' is provided, fun should be a function which expects 2 inputs:
%           [ ...] = fun( par(:, index_range) , fields(:, index_range) )
%           The resulting derivatives should (still) be w.r.t. par.
% data : A N-dimensional matrix in which the data of each voxel is stored in a single
%        column (first dimension, data(:,i)). Any datatype convertable to double is allowed.
%        Dimensions after the first one are treated as spatial dimensions.
% par0 : Initial values from which the local optimization starts. 
%        Can be a single column vector that is used as initial value for all voxels
%        or an (size(par,1) x spatialsize) array with an initial
%        pararameter for each voxel. 
%        If no initial value is provided (par0=[]), a value is requested by calling:
%        par0 = fun([]);
%
% Options - value pairs: (May also be in a scalar structure with the options as fieldnames)
% 'noiseLevel'         : default = [];  => estimated from residue.
%                        Scalar, or column vector with one element for each MR image, 
%                        or N-dimensional matrix with one element per voxel with the noise 
%                        standard deviation of the rician or normal MR intensity distribution. 
%                        If project is specified and noiseLevel is empty, one scalar noise level is estimated
%                        from the residue between each project and corresponding data.
%                        if noiseLevel is non-empty, you can specify a global noise level, a noise level for each 
%                        image, a noise level for each project, or for each data point individually. 
%                        (The last one cannot be 'refined')
%                        Specify the desired structure with a cell array. 
% 'mayRefineNoiseLevel': boolean. default: false if 'noiseLevel' is provided, true otherwise.
%                        if true, the noiseLevel is estimated after the initial fitting round
%                        and if it deviates much (see acceptNoiseLevelAdjustFactors) the fitting
%                        is redone
% 'numPDFoptpar'       : The main reason of this option is to specify if the noise level should
%                        be estimated globally or per voxel.
%                        0 (default): global noise level estimate, or using fixed prespecified noise level
%                        1 : for the default 'imageType', an initial noise level estimate
%                           should be included as last 'row' in par0. The final estimated noise level of
%                           each voxel is returned in the same position in par.
%                           (Since in this case the noise level estimate is part of par, regularization 
%                           can be applied on this estimate as well).
% 'imageType'          : allowed: {'magnitude' / 'rician'}, 'real'/'normal', 'complex', 'phase'
%                                 'StudentTx'
%                        Specifies the type of MR image that is provided. 
%                        - For magnitude images (the default), the Rice distribution
%                          of the image data is taken into account. 
%                        - Real indicates that a single channel (either real or imaginary) is selected 
%                          and thus the noise distribution is gaussian. Since 'Real' is a bit faster it 
%                          might also be usefull for magnitude images when the SNR>3 for all voxels.
%                        - Complex indicates complex images, with complex valued gaussian noise
%                        - phase indicates a phase image. The fitting is done in the complex
%                          domain to avoid wrap around problems. When a complex image is available,
%                          typically using that is to be prefered.
%                        - StudentTx , where x is a printed real number (num2str). Default: x=4
%                          Assumes that the data have a studentT distribution with x degrees of freedom.
%                          This is more robust to outliers in data than the normal distribution.
%                          Typically it is slower in the optimization.
% 'gradObj'            : {'on'}, 'off'
%                        specifies if fun can compute the gradient of the function value.
% 'hessian'            : 'on',{'off'}
%                        specifies if fun can compute the hessian of the function value.
% 'maxTracesInBlock'   : default: [inf inf inf]
%                        The maximum number of elements in index_range that is supported by 
%                        fun (when computing function value (el. 1), gradient (el. 2) and hessian (el.3) ). 
%                        Set to 1 if fun can predict the MR values of a single trace only.
% 'blockSize'          : N element vector with (maximum) size of block that should be computed 
%                        at once. Default: 2 in each dimension.
% 'spatialRegularizer' : default [] : no spatial regularization
%                        Specifies the spatial regularization that should be performed
%                        'TV' or 'total variation' : total variation penalty
%                        'laplacian'  : regularize with the norm of the laplacian of the
%                                       estimated parameters.
%                        Also specify: 'voxelSpacing' and 'spatialRegularizerWeights'
%                        Alternatively spatialRegularizer might be a function or 2 element cell 
%                        array with functions such that the regularization penalty is 
%                           [f,g,h] = spatialRegularizer(x) 
%                        or [f,g,hessinfo] = spatialRegularizer{1}(x) with Hessian*x = spatialRegularizer{2}(hessinfo, x)
%                        with x being a 'box' round the voxels that are currently being optimized. 
%                        so size(x,1)=size(par,1) , ndims(x)=ndims(par)
%                        Specify 'spatialRegularizeBorder' to indicate how many voxels around the
%                        ones being optimized should be put inside the 'box' to avoid erroneous border effects. 
% 'voxelSpacing'       : N-1 element vector specifying the voxel spacing of the spatial dimensions.
%                        Only used for the default regularization options ('TV' or 'laplacian'). 
% 'spatialRegularizerWeights' : size(theta,1) vector or size(theta,1) x size(theta,1) matrix 
%                               with the weights of each parameter in the regularization
%                               Only used for the default regularization options ('TV' or 'laplacian'). 
%                               The distance d between parameters of neighboring voxels i and j that is 
%                               used for the regularization is computed by
%                                  vecd = theta(:,i)-theta(:,j)
%                                  d = sqrt( vecd' * spatialRegularizerWeights * vecd );
% 'validateFullFunction' : boolean, default = false
%                          Perform extensive validation of the full function, including regularization
%                          if applicable, to check the gradient/jacobian and hessian.
%                          Can be used to check the derivative/hessian computations.
%                          Only use on small problems as the full hessian or jacobian is
%                          evaluated and many function calls are used. 
%
% Some other options. These are typically documented quite well in the code of this function (where the default values are provided).
%  Some occasionally important ones are:
%    maxIter, tolFun, validateFullFunction, logPDFfun, parameterPrior,  project, projectParameters, 
%    startBlockPos, blockSize, blockOverlap, blocksOptimizeDirection, blocksOptimizeDimOrder,
%    computeRange, mask, fields, constraints_A, constraints_b, 
%    onlyKeepnBestStartPoints, initialValueSpecifierVect, extraInitialValues, save_initializationorder
%  Still other parameters, trying to exhaustively list options that (may) make sense to specify:
%    spatialSize, maxIter_align, tolFun_align, projectOptimization.outerloopTolerance, progressbarUpdateStep,
%    acceptNoiseLevelAdjustFactors, maxNoiseLevelAdjustRepeats, spatialRegularizerLS,
%    maxErrorNotifications, parameterPrior_mu, parameterPrior_Sigma
%
% Some known 'fun' functions:
%  verysimpletestfitfun, MLE_Haytonfit_*, MLE_T1_d,  MLE_T12_w_TE_TR, DifusionTensor_Apred_m
%
% OUTPUTS:
% par     : N dimensional matrix with the Maximum Likelihood fitted parameters. 
%           par(:, i) are the parameters of voxel i.
% residue : data - predicted intensities.
% finalOptions : Scalar structure with the final options, including noise level.
%
% NOTE: This function uses the optimization toolbox. Unfortunately there 
%       are a few bugs in the version of this toolbox that I have that need fixing.
%       See comments just below this help text in this file for the fixes.
%
% created by Dirk Poot (d.poot@erasmusmc.nl), Erasmus Medical center, Rotterdam

% Original version: 22-3-2011
% Based on ComputeDT which was originally created 14-8-2007
%  4-10-2011 : add support for estimating noise level for each voxel (in the main ML optimization)
% 13-12-2013 : adapted by Gwendolyn Van Steenkiste, Introduced a block
% metho that doesn`t use the PSF.
% ************************************************************************************************
% ************************  PATCHES FOR THE OPTIMIZATION TOOLBOX *********************************
% ************************************************************************************************
%  NOTE: On windows >= vista you probably need administrator rights to modify these files. 
%
%  file: [matlabroot       '\toolbox\optim\optim\private\snls.m']
%     or:[matlabroot 'R2011b\toolbox\shared\optimlib\snls.m']
%  line 443 Original  : JACOB = sparse(A);     % A is the Jacobian, not the gradient.
%           Corrected : JACOB = A;	           % A is the Jacobian, not the gradient.
%  Reason/impact: This is the parsing of the jacobian information to the jacobian output argument. 
%                      It is not used by the optimization routines. The original fails when the jacobmult function 
%                      uses anything else than a double 2D matrix as information about the jacobian (as we do).
%                      The adjustment might potentialy break other users code that relies on automatic conversion
%                      of their own jacobinfo output from double to sparse. However, this is not likely to occur
%                      since when you needed sparse jacobinfo, you probably already made it sparse to begin with.
%
%  file: [matlabroot       '\toolbox\optim\optim\private\pcgr.m']
%     or:[matlabroot 'R2011b\toolbox\shared\optimlib\pcgr.m']
%  line 108,109 Original: wbar = RPCMTX'\r(ppvec);
%                         w(ppvec,1) = RPCMTX\wbar;
%               Hacked  : if isa(ppvec,'function_handle')
%                            w = ppvec(r, RPCMTX);
%                         else 
%                            wbar = RPCMTX'\r(ppvec);
%                            w(ppvec,1) = RPCMTX\wbar;
%                         end;
%  Reason/impact: This is a hack to allow a truly custom preconditioner. I need this since a normal preconditioner 
%                 would/might consume too much memory (even as sparse matrix) and would be much slower.
%                 ppvec is (/should) not be used explicitly anywhere else, so this hack should not cause any compatibility problems.
%                 You can use it 'from' the preconditioner function; the preconditioner function should return a function_handle 
%                 in the variable ppvec and RPCMTX should contain any options that are required.
%
% ************************************************************************************************
% ************************  END OF PATCHES FOR THE OPTIMIZATION TOOLBOX **************************
% ************************************************************************************************

residue = [];
% Create structure with default settings. Each of these can be overridden by a option-value pair.

opts.function = fun;          % function predicting intensities. Since this is one of the 3 special inputs, it typically should not be overridden by a option-value pair.
opts.function_jacMul = @function_jacMul;
                              % Function that multiplies jacobian of function with x
                              % J = nimg  x  ntraces  x  nparameters
                              % x = nparameters x ntraces
opts.function_jacMulAdj = @function_jacMulAdj;
                              % Function that multiplies adjoint of jacobian of function with x
                              % J = nimg  x  ntraces  x  nparameters
                              % x = nimg  x  ntraces
opts.maxfunArgsOut = [];      % maximum number of output arguments; counting: [ predicted_magnitude, gradient, hessian ]
                              % typically use 'gradObj' and 'hessian' options to fill this.
                              % Specifying this option disables all usage of the options 'gradObj' and 'hessian'
opts.funHessian_I = [];       %  [I,J] = find(triu(ones(size(par,1))));
opts.funHessian_J = [];       %  
opts.spatialSize = [];        % Spatial size of images, 
opts.numParam  = [];          % number of parameters in theta. Always overwritten.
opts.numImages = [];          % Number of MR images. (=size(data,1))
opts.noiseLevel = [];         % noise level in the MR images
opts.tolFun_meanLL_NoiseLevel = .1;% If mean LL improvement due to noise level update is less than this value, iterations are stopped.
opts.imageType = 'magnitude'; % magnitude, real, complex, phase.
% opts.initialTheta = [];       % initial estimate
opts.gradObj = 'on';          % can function compute gradient ?    
opts.hessian = 'off';         % can function compute hessian ?
opts.maxIter = 20;            % maximum number of iterations in optimizations.
opts.tolFun  = 1e-4;          % convergence tolerance for function. Fixed since we're optimizing likelihood values.
opts.projectOptimization.maxIter = 20;    % maximum number of iterations in optimization of alignment .
opts.projectOptimization.tolFun  = 1e-4;    % convergence tolerance for function. Fixed since we're optimizing likelihood values.
opts.projectOptimization.outerloopTolerance = [];    % Tolerance for the log likelihood update by the alginment.
                              % If the alignment changes the log likelihood less than this value, repeat estimation iterations are stopped.
opts.projectOptimization.skip = false; % boolean that allows to skip optimization of the projection parameters.
                             
opts.maxTracesInBlock = [];   % default maximum number of traces that function can compute simultanuously.
                              % defaults to [inf inf inf]
                              % Set to the maximum number of elements in index_range that is supported 
                              % (and efficiently evaluated) by fun when computing:
                              %  function value (el. 1), gradient (el. 2) and hessian (el.3) . 
opts.progressbarUpdateStep = 32; % number of traces to compute before updating the progressbar (only for speed)
                                 % Zero indicates no progress information display.
opts.mayRefineNoiseLevel = [];   % allow refinement of noise level?
opts.acceptNoiseLevelAdjustFactors = [.8 1.25]; % min, max adjust factor of noise level that is acceptable (doesnt require recompute).
opts.maxOuterloopIterations = 5;           % maximum number of adjustments of the noise level. Each adjustment requires full recomputation(/adjustment) of all parameters.
opts.doRegularize = [];       % boolean automatically set, conditional on filling 'spatialRegularizer' 
opts.spatialRegularizer = {};   % see help text above.
opts.spatialRegularizerLS = {}; % LS version of regularization. format for lsqnonlin (regularization is sum of squares. Function output is non squared terms of the sum, second output is jacobian).
                                % second element in the cell array is jacobian multiply function.
opts.spatialRegularizeBorder = 0; % border size of spatial regularizer.                                
opts.doComputeDerivativeRegularizationScale = false; % if true computes derivative w.r.t. regularisation scaling. 
                                                     % Only applies when regularizing
                                                     % See problem as   theta = arg min f(theta; data) + lambda * R(theta)
                                                     %  Solved with lambda ==1 , compute derivative of theta w.r.t. lambda.
opts.voxelSpacing = [];         % voxel spacing in each of the spatial dimensions. 
opts.spatialRegularizerWeights = 1; % dummy weight. Value is very important when using spatial regularization!
opts.validateFullFunction = false;  % Validate gradient/jacobian and hessian. Only use on small problems!!

opts.logPDFfun = [];                % logarithm of the PDF of the data with the 'mean' values predicted by fun.
                                    % 4 argument function that evaluates the PDF for each predicted magnitude:
                                    % logPDFfun( data, A, sigma , derivative_selector) 
                                    % derivative_selector : 3 element logical vector; when true in 
                                    % 1st element: compute derivatives w.r.t. data
                                    % 2nd element: compute derivatives w.r.t. A
                                    % 3nd element: compute derivatives w.r.t. sigma
                                    % default: filled in by 'imageType'
                                    % Output: full or summed over images (but not over voxels!)
opts.numPDFoptpar = 0 ;         % Number of parameters in theta (last elements) that specify the PDF of the data.
                                % Usually 0 if noise level is estimated globally. Set to 1 if you want to estimate
                                % the noise level standard deviation for each voxel independently. 
                                % (logPDFfun should accept 3 inputs and also compute derivative (&hessian) with respect to sigma)
opts.project = [];              % [numel(data_in) x 5] cell array with projector functions, i.e. simulations of (MR) image(/data) acquisition.
                                % Each element of the cell array data corresponds to a single row of project
                                % and a single contrast image (specific element in output of fun)
                                % Each row 'k' consists of:
                                %  [  @projData, @gradientImgMul, @gradientImgAdjMul, @gradientParMul, @gradientParAdjMul ] 
                                % where 
                                %  [ im, g ] = projData( img{k} , projectParameters{k} )
                                %                projection function that is 'almost' linear in img{k} 
                                %                (second derivative of projData is ignored during optimization)
                                %  [ gI, g ] = gradientImgMul( I , g )
                                %                gI( a ) = sum_b  d im( a ) / d img{k}( b ) * I( b )
                                %  [ J, g ]  = gradientImgAdjMul(gJ, g ) 
                                %                J(  b ) = sum_a  d im( a ) / d img{k}( b ) * gJ( a )
                                %  [ gP, g ] = gradientParMul( P, g )
                                %                gP( a ) = sum_b  d im( a ) / d projectParameters{k}( b ) * P( b )
                                %  [ P_, g ] = gradientParAdjMul(gP_, g) 
                                %                P_( b ) = sum_a  d im( a ) / d projectParameters{k}( b ) * gP_( a )
                                % with:
                                %  img{k} : image in reconstruction volume formed from k-th element of ouput of 'fun'
                                %  projectParameters{k} : parameters that ajust the projection and which should be optimized as well.
                                %  im     : img{k} projected to data{k}; From the predicted image simulate the acquired data.
                                %  g      : fully user specified; should contain information to be able to multply with gradient afterwards.
                                %           Preferably 'projData' should not perform 'heavy' calculations to compute g, 
                                %           as it's output g might (occasionally) be unused.
                                %           If 'heavy' computations can be reused by several invocations of the gradient multiply functions, 
                                %           these may update g.
                                %  I,J    : An vector with the size of img{k}               that should be/was multiplied with the gradient
                                %  gI, gJ : An vector with the size of im                   that should be/was multiplied with the gradient
                                %  P, P_  : An vector with the size of projectParameters{k} that should be/was multiplied with the gradient
                                %  gP,gP_ : An vector with the size of im                   that should be/was multiplied with the gradient
                                % 
                                % Alternatively: when multiple 'projects' are performed for each image, 
                                % project{k,1} (=>'projData') might be a cell array with function handles, one for each projection. 
                                % 'projectParameters{k}' should be a cell array with corresponding size, specifying the parameters 
                                % for each project. The gradient multiplication functions may be a cell array of corresponding size,
                                % or a single function. Data_in{k} should also be a cell array with the data for each projection.
                                %
opts.projectSelectSource =[];   % [numel(data_in) x 1] cell array with the integer vector selecting which projected images are provided to which projection function (&parameters).
                                % default : projectSelectSource{i} = [i];
                                %           => predicted image i is provided to project i.
                                % if multiple predicted images are selected, they are concatenated such that 
                                % size(img{k}) = [spatialSize numel(projectSelectSource{i})]
                                % Let j = projectSelectSource{i}(k)
                                %  if j>0 : select image from row j of output of function
                                %  if 0 < (-j) <= size(par,1) : select image from par(-j)
                                %                              NOTE: this is mainly usefull to select the noise level parameter(s) (when numPDFoptpar>0)
                                %  if (-j) >= size(par,1)     : select field(-j-size(par,1))
                                %                              NOTE: this is provided for convenience only. 
                                %                                    Including the (static) field within the project function will always allow improved performance.
opts.projectIsLinear = false;   % if true: project should be exactly a linear operator. 
opts.projectScaledPartPSF =[];  % Scalar structure with the following fields (as function name), each containing a (anonymous) function with 
                                % the following signature:
                                % - prepareGlobal:  [ hLLifo ] = fun( dLLdf_full, opts )
                                %           Prepares the hLLifo structure. With this structure, for each processing block,
                                %           prepareRegion constructs a list of voxels in which the PSF is nonzero 
                                %           and prepareBlock constructs the hLLifoBlk structure. 
                                %           hLLifo is not used for anything else.
                                %           INPUTS:
                                %             opts.project : the projection functions as input to fit_MRI
                                %             opts.projectGrad : a cell array with the required gradient info as returned by 
                                %                                the third output (hinfo) of 'adjustProjectParametersCritFun'
                                % - prepareRegion: [ selLrgPSF , hLLifoBlk0 ] = fun( sellin, hLLifo )
                                %           Prepares a list of voxels in which the PSF is(/might be) nonzero.
                                %           INPUTS:
                                %             sellin : a row vector with linear indices (sub2ind) of the reconstruction volume
                                %                      of the voxels that are optimized within the current block.
                                %             hLLifo : output of PSFMulprepareGlobal_default (opts.projectScaledPartPSF.prepareGlobal)
                                %           OUTPUTS:
                                %             selLrgPSF : n x ndims   matrix with full indices of voxels that (migth) have non-zero PSF.
                                %                         These points should be unique, but may be outside the reconstruction 
                                %                         volume (they are clipped afterwards).
                                %             hLLifoBlk0 : may contain something; passed to 'prepareBlock' and not used otherwise.
                                % - prepareBlock:  [ hLLifoBlk , psfdiag ] = fun( sellin, selLrgPSF, hLLifo, hLLifoBlk0)
                                %           For each processing block, prepares hLLifoBlk, which is a structure with which mulsel
                                %           and mulfull can perform their multiplications with the psf 
                                %           INPUTS:
                                %             sellin : a row vector with linear indices (sub2ind) of the reconstruction volume
                                %                      of the voxels that are optimized within the current block.
                                %             selLrgPSF : a row vector with linear indices (sub2ind) of the reconstruction volume
                                %                      of the voxels of which the 'dLLdf_full' needs to be updated. 
                                %             hLLifo : output of PSFMulprepareGlobal_default (opts.projectScaledPartPSF.prepareGlobal)
                                %             hLLifoBlk0 : second output of prepareRegion.
                                %           OUTPUTS:
                                %             hLLifoBlk : Structure with which PSFMulSelectedBlock_default (opts.projectScaledPartPSF.mulsel)
                                %                         and PSFMulLargeBlock_default (opts.projectScaledPartPSF.mulfull) can do their 
                                %                         multiplications with the psf.
                                %             psfdiag   : scalar, (numimg x 1)  vector or (numimg x 1 x numvoxelsinblock) array
                                %                         with (approximate) diagonal of the PSF. Used for desgining the (voxelwise) preconditioner.
                                % - mulsel:  [ Hx ] = fun( x, hLLifoBlk )
                                %           Multiply x with the point spread function.
                                %           INPUTS:
                                %             x : numimg x numvoxelsinblock
                                %             hLLifoBlk : output of PSFMulprepareBlock_default
                                %           OUTPUTS:
                                %             Hx : numimg x numvoxelsinblock
                                %                  Hx(i,j) = sum_klm d2 L /d im_i(k) d im_i(l) * (d im_i(k)/d f_i(j) ) * (d im_i(l)/d f_i(m) ) * x(i,m)
                                %                  with L = total log likelihood, im_i = project{i,1}( f_i )
                                %                     with f = fun( theta ), f_i = f(i,:,:,:)
                                %                  and j = sellin
                                %                  Note that we ignore the non linearity of project, which typically is 
                                %                  small and thus can be ignored. However, if you can (easily) include it for your problem, 
                                %                  please do so.
                                % - mulfull
                                %            Same as mulsel, but now with j = selLrgPSF (as input to 'prepareBlock')
                                % 
                                % The default implementation computes the PSF numerically.
                                % This is slow (using many gradientImgMul and gradientImgAdjMul), 
                                % but exact (if we ignore the second derivative of project).
                                % If the PSF is (almost) constant (invariant) over the entire image, please specify this
                                % by setting 'projectPSFspatiallyInvariant'. 
opts.projectPSFnonzeroOffsets = [];    % numel(data_in) cell array in which each element specifies where 
                                % the PSF of the corresponding project is non zero with a
                                % n x numel(spatialSize) matrix specifying the offsets 
                                % of (potentially) non zero psf elements.
opts.projectPSFspatiallyInvariant = false; % default = 0 (false): cant assume PSF is constant (see projectScaledPartPSF)
                                % if 1 (=true): the PSF does not depend on spatial location. Note that this also requires that the 
                                %               second derivative of the likelihood function is spatially constant, which is only 
                                %               exactly true for 'imageType' = 'normal' 
                                % if .5       : The PSF is (approximately) constant within a block.
                                %               The PSF of the center of each block is used (except when this is too close to an edge of the reconstruction volume). 
                                %               This is especially more efficient than 0 (default) when the blocksize is large. 
                                %               However, for this option, the blocsize should not be larger than half the reconstruction volume size in any dimension. 
                                %               (which typically is not a good situation anyway, since optimization in larger blocks is more difficult).
opts.project_PSFscaleBlock  = 1; % factor with which the PSF is multiplied during opmization. If >1 => stronger second derivative => smaller step.
                                % Smaller steps mean less prone to fail due to oscilations, but less efficient.
                                % Should not be <1 (causes instability in optimization). 
                                % Note for writers of 'projectScaledPartPSF' functions: you dont need to take this option into account, it is handled for you.
opts.projectScaleDeltaPSFupdate = 1; % scale the update of the gradient due to the PSF
                                % should be 1 if PSF is exact. 
                                % Note for writers of 'projectScaledPartPSF' functions: you dont need to take this option into account, it is handled for you.
opts.projectEqualPSFgroups = [];% cell array with integer vectors. Each vector specifies the 'project' elements with equal PSF.
                                % Default: each project has a unique PSF.
opts.projectPSFrelativecutoff=.01; % default: Elements that in magnitude are below 1% of the maximum PSF magnitude are ignored.
opts.adjProject = [];           % cell column with adjProjection functions. 
                                %   [data , dData_dAlignpar ] = adjProject{i}( data_in{i}, projectParameters{i} )
                                %   function that warps the data_in to ND MR images, dependent on alignpar
                                % or a 3 column cell array with 
                                %  adjProject(i,:) = {@projData, @gradientParMul, @gradientParAdjMul } 
                                % where 
                                %  [ img, g ] = projData( data{k} , projectParameters{k} )
                                %                projection function that is 'almost' linear in img{k} 
                                %                (second derivative of projData is ignored during optimization)
                                %  [ gP, g ] = gradientParMul( P, g )
                                %                gP( a ) = sum_b  d img( a ) / d projectParameters{k}( b ) * P( b )
                                %  [ P_, g ] = gradientParAdjMul(gP_, g) 
                                %                P_( b ) = sum_a  d img( a ) / d projectParameters{k}( b ) * gP_( a )
                                % with:
                                %  data{k} : the k-th element in the source image data cell array.
                                %  projectParameters{k} : parameters that ajust the projection and which should be optimized as well.
                                %  img    : data{k} projected to reconstruction space.
                                %  g      : fully user specified; should contain information to be able to multply with gradient afterwards.
                                %           Typically 'projData' should not perform 'heavy' calculations to compute g, 
                                %           as it's g output might (occasionally) be unused.
                                %           If 'heavy' computations can be reused by several invocations of the gradient multiply functions, 
                                %           these may update g.
                                %  P, P_  : An vector with the size of projectParameters{k} that should be/was multiplied with the gradient
                                %  gP,gP_ : An vector with the size of img                  that should be/was multiplied with the gradient
                                % 
                                
opts.projectParameters = {};    % cell array of same size as project/adjProject with for each project the initial project adjustment parmeters.
                                %  These parameters are optimized further. 
opts.projectParameterPrior = {[]}; % cell array with functions that compute the log likelihood of the prior distribution of the project parameters
                                 % Same structure as project; 1 or 2 functions (columns);
                                 % { logprior_fun [, hessmulfun] }
                                 %  [f, g , H/Hinfo] = logprior_fun( projectParameters{k} )
                                 % if hesmulfun is provided:
                                 %      H * x == hessmulfun( Hinfo, x)
opts.projectGrad = [];         % internal parameter, should not be provided unless you are debugging AND know what you are doing.
opts.explicitHessian = false;  % internal parameter.
opts.optimizeBlocks  = true;  % true: optimise theta in blocks, 
                              % false: optimise all voxels at once (I think this can only be usefull when the PSF is non local)
opts.maxPCGiters     =  [];   % Maximum number of PCG iterations in each quadratic sub-problem.
                              % Default based on blocksize, #parameters, per voxel and number of spatial dimensions.
                              % The default is too low when spatial regularization 'causes' a badly conditioned
                              % hessian. (The default should be fine when the voxel function has a badly conditioned hessian)
                              % This typically happens when the regularisation has a strong (blurring) effect.
opts.computeRange    =   [];  % 2 x ndims matrix that specifies the [start;end] voxel which might be updated in each dim.
                              % Allows to select a box, especially usefull when you do not want to reduce the image size
                              % and/or specify a mask. Or when the nonzero part of the mask is much smaller than the image.
                              % By default the minimum range that still selects all voxels selected by the mask is used.
opts.startBlockPos   =   [];  % start position of first block in index coordinates, specify for each spatial dimension, might be negative (smaller first block)
opts.blockSize       = [];  % (maximum) block size, specify for each spatial dimension
opts.blockOverlap    =  []; % overlap in voxels between neighboring blocks. Should be >=0 && <blockSize. specify for each spatial dimension.
opts.blocksOptimizeDirection = 1; % direction in which the blocks are evaluated in each spatial dimension
                                   % positive : low to high
                                   % negative : high to low.
opts.blocksOptimizeDimOrder = []; % order of dimensions along which the blocks are optimized. 
                                  % default: 1:nspatialdims;      
                                  %     meaning: blocks adjacent in first spatial dimension are optimized right after each other. 
                                  % Alternatives should be a permutation of this default.
opts.periodicDomain = false;    % logical scalar or nspatialDims element vector.
                                % if true, the domain is treated as periodic.
                                % That is: when (e.g. the regularization) would access
                                % the voxel theta(:, ..., end+1, ... )
                                % the voxel theta(:, ...,     1, ... ) is used
                                % (instead of the default end case handling
                                % of the regularisation).
opts.maxErrorNotifications = 5; % Maximum number of errors that are shown during single optimization of all blocks. 
opts.encounteredErrors =  0 ;   % number of errors encounterd untill now.
opts.parameterPrior.fun = [];   % 'normal' or function that computes the prior log likelihood of a parameter vector. 
                                % To correct the noise level estimate when numPDFoptpar==1 use param_prior_debias_sigmaest
                                % fun should be a 1 input function 
                                % [f,g,h] = fun( par ) ;   
                                %      par : nparameters x n
                                %      f   : scalar (cost)
                                %      g   : numel(grad_linindex) x n
                                %      h   : numel(hess_linindex) x n 

opts.parameterPrior.mu = [];    % if parameterPrior= 'normal' : mean location of the prior
opts.parameterPrior.sigma = []; % if parameterPrior= 'normal' : covariance matrix of the normal prior distribution. Should be positive definite.
opts.parameterPrior.gradient_linindex = []; % linear index of the parameters for which the gradient is computed by parameterPrior.fun
                                        % default: 1:nparameters
opts.parameterPrior.hessian_I = []; % row index of the parameters for which the hessian is computed by parameterPrior.fun
opts.parameterPrior.hessian_J = []; % column index of the parameters for which the hessian is computed by parameterPrior.fun
                                        % default: [I,J]=find( triu(ones(nparameters) ) )
opts.mask = [];                 % Update only computed within this mask. Default= all voxels.
opts.fields = [];               % Extra parameter fields for function. Same format as resulting x
                                % first dimensions contains all diferent fields, remaining dimensions are spatial dimensions of size(x)
                                % fun should accept 2 inputs.
opts.constraints = [];          % scalar constraints specifying structure with fields:
                                %  lb, ub, Aineq, bineq, Aeq, beq, nonlcon
                                %  which specify constraints for each voxel, e.g.:  Aineq * theta(:,i) <= bineq
                                %  see fmincon for more details on the constraints.
opts.onlyKeepnBestStartPoints = 1; % Number of start points to keep for full optimization of block. 
                                   % if >1 : optimize for (at most) this many start point and select best after optimization.
                                   % If more than this many start points are given, the function that is optimized
                                   % is evaluated for all initial values. The resulting function values are sorted 
                                   % and only the best onlyKeepnBestStartPoints start points are selected for actual optimization. 
                                   % After optimization the very best is selected and put in the output.
opts.initialValueSpecifierVect = [1 0 0 0]; % vector that specifies which initial values should be selected.
                                   % Element i of this vector means, when true 
                                   % 1 : select current value of theta for point block. (This should almost always be on)
                                   % 2 : select the N previously computed neighboring blocks (i.e. takes blocksOptimizeDirection into account)
                                   % 3 : select the N neighboring blocks not yet computed (i.e. takes blocksOptimizeDirection into account)
                                   % 4 : select all 'extraInitialValues'
                                   % 5 : select elements from 'initialValueRange'
                                   %        Specifically, the number specified here is the number of
                                   %        elements taken from a Hayton sequence over the
                                   %        parameterrange specified by 'initialValueRange'
                                   % For 2 and 3 the mask is not taken into account. if blockOverlap>0 then 
                                   % the used blocks overlap as well, which they also might at the edge cases (depends on startBlockPos)
opts.initialValueRange = [];      % nparameters x 2 matrix with lower and upper bound used for initial values
opts.linearParameters = [];       % (vector with) indices of the scaling parameter(s) in fun
                                   %  E.g. theta = [A R1 B R2 C1 C2]'
                                   %       fun(theta) = A * exp( R1 * x ) + B * sin( R2*Y + C1) + C2
                                   % then linearParameters = [1 3 6]  (or [])
                                   % Used (only) in combination with
                                   % initial values. (Both extraInitialValues and initialValueRange )
                                   % (Implicitly) performs a least squares estimation for these parameters, which
                                   % allows a substantial reduction in the number of initial paramater vector
                                   % required to cover the parameterspace. 
opts.extraInitialValues =[];     % Provide a cell array with in each element an initial 1-column-vector. 
                                 % Enable using of this setting with: initialValueSpecifierVect(4)=true
                                 
opts.save_initializationorder = false; % If true: Save order of initializations of blocks for which the 
                                       % number of initializations exceeds onlyKeepnBestStartPoints.
                                       % Saved to fields 'saved_initializationorders' in the opts output.
                                       % Can be used to study if all provided intializations do get selected occasionally.
% opts.saved_initializationorders   % Not an input.
opts.optimizer = 1;     % The optimizer that is (prefably) used (in global optimization)
                        % 1 fminunc/fmincon
                        % 2 congrad_nonlin
                        % 3 fmin_fast
                        % 4 new_MCMC_sample, => one run of Gibbs sampling a markov
                        %   chain with the Metropolis-Hastings algorithm to
                        %   update each point.
opts.computeLeverage = false; % boolean , if true : compute leverage factors instead of 
                              % performing any other computation.
                              % leverage factors indicate how much a change in a 
                              % data sample is absorbed by theta.
                              % Idealy this should be
                              %     size(theta,1)/size(data,1)  (<1)
                              % Indicating (in a way) that each sample is
                              % equally important. 
                        
opts.MakePSFApproxOfProj = true; % boolean, if true : uses the PSF of the projection
                                 % if false: the complete projection is
                                 % used.

opts.MakeProjBlockFun = []; % Function that addapts the projections so only the ROI is projected.


if isequal(opts.function,'eye')
    % override defaults for special (dummy) function. 
    opts.function = @functionEye;
    opts.function_jacMul = @functionEye_jacMul; 
    opts.function_jacMulAdj = opts.function_jacMul;
end;

% parse input arguments (set fields <option> of opts with corresponding <value>)
opts = parse_defaults_optionvaluepairs( opts, varargin{:});
optsOutputAdjust = struct; % empty structure.

% Conditional settings:
% Get size of reconstruction volume:
if isempty(opts.spatialSize) 
    if iscell(data_in)
        dummy = size(theta);
        if ~all(dummy(2:end)==1)
            % read size from theta, if not provided for just 1 column.
            opts.spatialSize = dummy(2:end);
        else
            if isempty(opts.fields)
                if isempty(opts.project) 
                    if isempty(opts.adjProject)
                        error('A cell array as data input is only allowed when projecting. (So provide project, adjProject, or an ''normal'' data array.)');
                    end;
                    warning('fit_MRI:SpatialSizeNotProvided', 'Efficiency: Please provide ''spatialSize''. Now an extra adjProject is used to obtain this value.');
                    dummy = opts.adjProject{1}( data_in{1}, opts.projectParameters{1}); % also crashes if adjProject is not specified, which here means that the provided options are invalid.
                else
                    error('fit_MRI needs to know the spatial size; Please provide "spatialSize", or an initial value that fills the desired space.');
                end;
                opts.spatialSize = size(dummy);
            else
                % read size from fields.
                dummysz = size(opts.fields);
                opts.spatialSize = dummysz(2:end);
            end;
        end;
    else
        dummy = size(data_in);
        opts.spatialSize = dummy(2:end);
    end;
end;
if isempty(opts.numImages)
    if iscell(data_in)
        if isempty(opts.projectSelectSource)
            opts.numImages = numel(data_in);
        else
            % Assume the maximum selected source image is output by fun and all 
            % output images are projected:
            opts.numImages = max([opts.projectSelectSource{:}]);
        end;
    else
        opts.numImages = size(data_in,1);
    end;
end;
szData = [opts.numImages opts.spatialSize];
nvoxels = prod(opts.spatialSize);
nspatialdims = numel(opts.spatialSize);
if isempty(opts.voxelSpacing)
    % set default voxelSpacing
    opts.voxelSpacing = ones(1,nspatialdims);
end;
% Check and possibly create computeRange, i.e. a box containing the voxels that get updated:
if ~isempty(opts.computeRange) && (~isequal(size(opts.computeRange),[2,nspatialdims]) || any(opts.computeRange(1,:)<1) || any(opts.computeRange(1,:)>opts.spatialSize))
    warning('fit_MRI:InvalidInput:ComputeRange','Invalid computeRange specified, using default');
    opts.computeRange = [];
end;
if isempty(opts.computeRange)  % no default, or invalid default provided.
    % set default:
    opts.computeRange = [ones(1,nspatialdims);opts.spatialSize];
    if ~isempty(opts.mask)
        % compact computeRange based on mask:
        tmpmask = opts.mask;
        if size(tmpmask,1)==1 && opts.spatialSize(1)~=1
            szmask = size(tmpmask);
            tmpmask = reshape(tmpmask,szmask(2:end));
        end;
        opts.computeRange = box_mask( tmpmask );
        if opts.computeRange(2,1)==0
            warning('fit_MRI:EmptyMask', 'Mask is empty, initial parameters returned.');
            residue = nan;
            return;
        end;
        clear tmpmask 
    end;
end;
% fill defaults for block processing:
if isempty(opts.startBlockPos)
    opts.startBlockPos = ones(1,nspatialdims);
end;
if isempty(opts.blockSize)
    opts.blockSize = 2*ones(1,nspatialdims);
end;
if isempty(opts.blockOverlap)
    opts.blockOverlap = zeros(1,nspatialdims);
end;
if isempty(opts.maxTracesInBlock)
    opts.maxTracesInBlock = nvoxels;
end;
if numel(opts.maxTracesInBlock)==1
    opts.maxTracesInBlock = repmat(opts.maxTracesInBlock,1,3);
end;
if numel(opts.blocksOptimizeDirection)==1
    opts.blocksOptimizeDirection = repmat(opts.blocksOptimizeDirection,1,nspatialdims);
end;
if isempty(opts.blocksOptimizeDimOrder) || ~isequal(sort(opts.blocksOptimizeDimOrder(:)),(1:nspatialdims)')
    if ~isempty(opts.blocksOptimizeDimOrder)
        warning('fit_MRI:InvalidOptionValue','The provided ''blocksOptimizeDimOrder'' is invalid; using the default order (1:nspatialdims).');
    end;
    opts.blocksOptimizeDimOrder = 1:nspatialdims;
end;
if numel(opts.periodicDomain)==1
    opts.periodicDomain = opts.periodicDomain(ones(1,nspatialdims));
end;
if isempty(opts.mayRefineNoiseLevel) 
    % if noiselevel provided : default do not adjust, otherwise adjust (first guess might be very wrong)
    opts.mayRefineNoiseLevel = isempty(opts.noiseLevel);
end;
if isempty(opts.projectOptimization.outerloopTolerance)
    opts.projectOptimization.outerloopTolerance = max(10, .01*prod(szData)); 
end;
if isempty(opts.doRegularize)
    opts.doRegularize = ~isempty(opts.spatialRegularizer) || ~isempty(opts.spatialRegularizerLS);
end;
if isempty(opts.maxfunArgsOut)
    if strcmpi(opts.hessian,'on')
        opts.maxfunArgsOut = 3;
    elseif strcmpi(opts.hessian,'off')
        if strcmpi(opts.gradObj,'on')
            opts.maxfunArgsOut = 2;
        elseif strcmpi(opts.gradObj,'off')
            opts.maxfunArgsOut = 1;
        else
            error('invalid value for option "gradObj", select "on" or "off", alternatively you may specify "maxfunArgsOut".');
        end;
    else
        error('invalid value for option "hessian", select "on" or "off", alternatively you may specify "maxfunArgsOut".');
    end;
end;
    
if isempty(theta) 
    % get initial theta, assuming fun provides default values:
    theta = opts.function([]);
end;
szTheta = size(theta);
if any(szTheta(2:end)==1 & szData(2:end)>1)
    theta = repmat( theta(:), [1 szData(2:end)./szTheta(2:end)]); % errors if incorrect size.
%     szTheta = size(theta);
end;
opts.numParam = szTheta(1);
if isempty(opts.funHessian_I)
    [opts.funHessian_I , opts.funHessian_J ]= find( triu( ones(opts.numParam ) ) );
elseif numel(opts.funHessian_I)~=numel(opts.funHessian_J)
    error('funHessian_I and funHessian_J should have the same size');
end;
% parse &check settings:
% theta = opts.initialTheta;

if opts.doRegularize && any(opts.maxTracesInBlock<=1)
    % if we are asked to regularize the solution, but the function is not vectorized, we should vectorize it.
    if isempty(opts.fields)
        opts.function = @(x) vectorizefun(opts.function, x);
    else
        opts.function = @(x, fields) vectorizefun(opts.function, x, fields);
    end;
    opts.maxTracesInBlock = max([1000 300 100],prod(opts.blockSize));
end;

% do we need to adjust the data to get real values:
dataAdjFun = @(x) x;
if strcmpi(opts.imageType , 'phase');
    opts.imageType = 'complex';
    if isempty(opts.fields)
        opts.function = @(x) phase2complex(opts.function, x, [], opts.funHessian_I, opts.funHessian_J);
    else
        opts.function = @(x, fields) phase2complex(opts.function, x, fields, opts.funHessian_I, opts.funHessian_J);
    end;
    dataAdjFun = @(x) exp(1i*dataAdjFun(x));
end;
if strcmpi(opts.imageType , 'complex');
    opts.imageType = 'real';
    if isempty(opts.fields)
        opts.function = @(x) splitsComplex(opts.function, x);
    else
        opts.function = @(x, fields) splitsComplex(opts.function, x, fields);
    end;
    % use splitsComplex to split data exactly the same way as the function output.
    dataAdjFun = @(x) splitsComplex(dataAdjFun , x);
end;
if isempty(opts.adjProject)
    if isempty(opts.project)
        % Conversion of data only required once if not projecting
        data = dataAdjFun(data_in);
        clear dataAdjFun;
    else
        % no conversion required when projecting.
        clear dataAdjFun;
    end;
end;
if strcmpi(opts.imageType , 'magnitude') && exist('data','var')
    wrongval = data<=0 | ~isfinite(data);
    if any(wrongval(:))
        if all(data(wrongval)==0)
            warning('fit_MRI:zeroDataValue','data cannot be zero for magnitude ML fit, zero values adjusted to 1');
        else
            warning('fit_MRI:InvalidDataValue',['zero, negative, infinite and/or NAN data found in ' num2str(nnz(wrongval)) ' elements, invalid for magnitude MR; adjusted to 1']);
        end;
        data(wrongval) = 1;
    end;
end;    

% Set logPDFfun based on image type:
if any( strcmpi(opts.imageType , {'magnitude', 'rician'}) )
    opts.logPDFfun =  @logricepdf;
    optsOutputAdjust.logPDFfun = 'logricepdf'; % overwrite functions in output, as otherwise the entire workspace gets copied as well.
elseif  strcmpi(opts.imageType , 'rician_logsigma')   
    opts.logPDFfun =  @logricepdf_logsigma;
    optsOutputAdjust.logPDFfun = 'logricepdf_logsigma'; % overwrite functions in output, as otherwise the entire workspace gets copied as well.
elseif  strcmpi(opts.imageType , 'rician_clamp')
    opts.logPDFfun = @(data, A, sigma, derivsel) clampargCallFun( @(data2) logricepdf( data2, A, sigma , derivsel) , data, 1);
    optsOutputAdjust.logPDFfun = 'logricepdf'; % overwrite functions in output, as otherwise the entire workspace gets copied as well.
elseif any( strcmpi(opts.imageType , {'real' ,'normal'}) ) 
    opts.logPDFfun = @lognormpdf;
    optsOutputAdjust.logPDFfun = 'lognormpdf'; % overwrite functions in output, as otherwise the entire workspace gets copied as well.
    opts.tolFun_meanLL_NoiseLevel = inf; % Noise level does not influence fit location; so we don't need to iterate. 
elseif numel(opts.imageType)>=8 && strcmpi(opts.imageType(1:8) , 'StudentT')
    nu = str2double(opts.imageType(9:end));
    if isempty(nu) || ~isfinite(nu)
        nu = 4;
    end;
    opts.logPDFfun = @(data, A, sigma, derivsel) logstudentTpdf(data, A, sigma, nu , derivsel);
    optsOutputAdjust.logPDFfun = ['logstudentTpdf_' num2str(nu)];
elseif isempty(opts.logPDFfun) || ~isa(opts.logPDFfun, 'function_handle')
    error('Unsuported imageType and logPDFfun does not contain an function.');
elseif ~isempty(opts.imageType)
    warning('fit_MRI:invalidImageType', 'option imageType is non empty but not recognised. The provided logPDFfun is used and imageType ignored.');
end;

% parse constraints: 
if ~isempty( opts.constraints )
    defaultconstraints = struct( 'lb',[],'ub',[],'Aineq',[],'bineq',[],'Aeq',[],'beq',[],'nonlcon',[]);
    opts.constraints = parse_defaults_optionvaluepairs( defaultconstraints, opts.constraints );
    if ~isempty(opts.constraints.lb)
        if numel(opts.constraints.lb)>opts.numParam || numel(opts.constraints.ub)>opts.numParam
            error('to many lower or upper bounds specified (opts.constraints.lb, or opts.constraints.ub)');
        else
            opts.constraints.lb = opts.constraints.lb(:);
            opts.constraints.lb(end+1:opts.numParam)=-inf; % due to replication for each voxel lb should be column vector and specified for each parameter, while specification for all parameters is not needed for fmincon
            opts.constraints.ub = opts.constraints.ub(:);
            opts.constraints.ub(end+1:opts.numParam)=-inf; % due to replication for each voxel ub should be column vector and specified for each parameter, while specification for all parameters is not needed for fmincon
        end;
    end;
end;

if ischar(opts.parameterPrior.fun)
    if strcmpi(opts.parameterPrior.fun,'normal')
        optsOutputAdjust.parameterPrior = opts.parameterPrior; % overwrite functions in output, as otherwise the entire workspace gets copied as well.
        neginvSigma = -inv( opts.parameterPrior.sigma );
        parameterPrior_mu = opts.parameterPrior.mu;
        hessIJlin = opts.funHessian_I + (opts.funHessian_J-1) * size(neginvSigma,1);
        opts.parameterPrior.fun = make1arg_anonfun(@lognormalVec, parameterPrior_mu, neginvSigma, hessIJlin);
    else
        error('unknown/unsupported prior distribution for the parameter vectors');
    end;
end;
if ~isempty( opts.parameterPrior.fun )
    if isempty(opts.parameterPrior.gradient_linindex)
        opts.parameterPrior.gradient_linindex = (1:opts.numParam)'; 
    end;
    if ~isequal(size(opts.parameterPrior.hessian_I), size(opts.parameterPrior.hessian_J)) || any(opts.parameterPrior.hessian_I(:)<0) || any(opts.parameterPrior.hessian_I(:)>opts.numParam) ||  any(opts.parameterPrior.hessian_J(:)<0) || any(opts.parameterPrior.hessian_J(:)>opts.numParam)
        error('parameterPrior.hessian_I and parameterPrior.hessian_J should have equal size and be integer vectors between 1 and numParam');
    end;
    if isempty(opts.parameterPrior.hessian_I)
        opts.parameterPrior.hessian_I = opts.funHessian_I;
        opts.parameterPrior.hessian_J = opts.funHessian_J;
    end;
    [test , opts.parameterPrior.hessian_write_funhessindex] = ismember( ...
            sub2ind(opts.numParam*[1 1], opts.parameterPrior.hessian_I, opts.parameterPrior.hessian_J) , ...
            sub2ind(opts.numParam*[1 1], opts.funHessian_I            , opts.funHessian_J )) ;    
    if any(~test)
        error('All hessian elements computed by the parameterPrior function should be contained in the function hessian.');
    end;
end;

if opts.doRegularize
    optsOutputAdjust.spatialRegularizer = opts.spatialRegularizer; % overwrite functions in output, as otherwize the entire workspace gets copied as well.
    [opts.spatialRegularizer, opts.spatialRegularizeBorder] = parse_spatialRegularizer( opts );
    
%     if isa(opts.spatialRegularizerLS,'function_handle')
%         opts.spatialRegularizerLS = {opts.spatialRegularizerLS ,[]};
%     end;
%     optsOutputAdjust.spatialRegularizer = opts.spatialRegularizer; % overwrite functions in output, as otherwize the entire workspace gets copied as well.
%     switch opts.spatialRegularizer
%         case {'TV','total variation'}
%             sp = [inf;opts.voxelSpacing(:)];
%             spatialRegularizerWeights = opts.spatialRegularizerWeights;
%             opts.spatialRegularizer = { make1arg_anonfun(@totalVariationVecRegularizer, sp , spatialRegularizerWeights,[],[],1) , totalVariationVecRegularizer([] , [], [], [], [], 2) };
%             opts.spatialRegularizeBorder = 1;
%         case {'laplacian'}
%             sp = [inf;opts.voxelSpacing(:)];
%             % TODO: make laplacian regularizer work with preconditioner and laplaceRegularizerHessMul
%             spatialRegularizerWeights = opts.spatialRegularizerWeights;
%             opts.spatialRegularizer = { make1arg_anonfun(@laplaceRegularizer, sp , spatialRegularizerWeights, 1 ) , make2arg_anonfun(@laplaceRegularizerHessMul, sp , spatialRegularizerWeights) };
%             opts.spatialRegularizeBorder = 2;
%         otherwise
%             if isempty(opts.spatialRegularizer)
%                 % assume opts.spatialRegularizerLS is provided
%                 spatialRegularizerLS = opts.spatialRegularizerLS{1};
%                 spatialRegularizerLS_JacMul = opts.spatialRegularizerLS{2};
%                 spatialRegularizerFun = make1arg_anonfun(@LSfun2normalfun, spatialRegularizerLS, spatialRegularizerLS_JacMul);
%                 spatialRegularizerHessMul = LSfun2normalfun('hessmulfun',spatialRegularizerLS, spatialRegularizerLS_JacMul);
%                 opts.spatialRegularizer = { spatialRegularizerFun , spatialRegularizerHessMul};
%             end;
%             if isa(opts.spatialRegularizer,'function_handle')
%                 opts.spatialRegularizer = {opts.spatialRegularizer , []};
%             elseif ~iscell(opts.spatialRegularizer) || ~isa(opts.spatialRegularizer{1},'function_handle')
%                 error('unsuported spatialRegularizer value');
%             end;
%     end;
end;
% set spatialRegularizeBorder for all spatial dimensions:
if numel(opts.spatialRegularizeBorder)==1
    opts.spatialRegularizeBorder = opts.spatialRegularizeBorder*ones(1,nspatialdims);
end;
    
opts.data_in = data_in;
if ~isempty(opts.project) || ~isempty(opts.adjProject)
    if isempty( opts.projectSelectSource )
        opts.projectSelectSource = cell(numel(opts.data_in),1);
        for k=1:numel(opts.projectSelectSource)
            opts.projectSelectSource{k} = k;
        end;
    end;
    if opts.numPDFoptpar~=0 && ~all( ismember(-size(theta,1)+(0:opts.numPDFoptpar-1), [opts.projectSelectSource{:}] ) )
        error('All noise level parameters should be included in projectSelectSource when projecting.');
    end;
end;
if ~isempty(opts.project) 
    if ~isfield(opts.projectScaledPartPSF,'prepareGlobal')
        opts.projectScaledPartPSF.prepareGlobal = @PSFMulprepareGlobal_default;
    end;
    if ~isfield(opts.projectScaledPartPSF,'prepareRegion')
        opts.projectScaledPartPSF.prepareRegion = @PSFMulprepareRegion_default;
    end;
    if ~isfield(opts.projectScaledPartPSF,'prepareBlock')
        opts.projectScaledPartPSF.prepareBlock  = @PSFMulprepareBlock_default;
    end;
    if ~isfield(opts.projectScaledPartPSF,'mulsel')
        opts.projectScaledPartPSF.mulsel        = @PSFMulSelectedBlock_default;
    end;
    if ~isfield(opts.projectScaledPartPSF,'mulfull ')
        opts.projectScaledPartPSF.mulfull       = @PSFMulLargeBlock_default;
    end;
    if isempty(opts.projectPSFnonzeroOffsets)
        opts.projectPSFnonzeroOffsets = cell(opts.numImages,1);
    end;
    if ~iscell(opts.projectPSFnonzeroOffsets)
        if ndims(opts.projectPSFnonzeroOffsets)>2 || size(opts.projectPSFnonzeroOffsets,2)~=numel(opts.spatialSize)
            error('Invalid value for projectPSFnonzeroOffsets. It should be a n x nSpatialDims matrix or a cell array with such matrices');
        end;
        opts.projectPSFnonzeroOffsets = {opts.projectPSFnonzeroOffsets};
        opts.projectPSFnonzeroOffsets(2:opts.numImages) = {opts.projectPSFnonzeroOffsets{1}(1,:)}; 
                    % fill, since otherwise an automatic search is launced. 
                    % Dont copy entirely to improve performance of 'unique', since currently, first the union 
                    %  of all non-zero voxels is taken and only afterwards the selection is made.
    end;
    if isempty(opts.projectEqualPSFgroups)
        opts.projectEqualPSFgroups = mat2cell(1:opts.numImages,1,ones(1,opts.numImages));
    end;
end;

% initialvalue parsing:
if numel(opts.initialValueSpecifierVect)>4 && (opts.initialValueSpecifierVect(5)>0) && ~isempty(opts.initialValueRange)
    if (opts.initialValueSpecifierVect(4)==0)
        opts.extraInitialValues = {};
    end;
    nin = opts.initialValueSpecifierVect(5);
    p = cell(size(opts.initialValueRange,1),1);
    primes = [2 3 5 7 11 13 17 19 23 29 31 37 41 43 47]; % if we run out of primes, there are too many parameters to do a search on...
    
    for k=1:size(opts.initialValueRange,1)
        p{k} = HaltonSequence( nin, primes(k), .5) * (opts.initialValueRange(k,2)-opts.initialValueRange(k,1))+opts.initialValueRange(k,1);
    end;
    opts.extraInitialValues = [opts.extraInitialValues  mat2cell(vertcat(p{:}),size(opts.initialValueRange,1),ones(1,nin) )];
    opts.initialValueSpecifierVect(4) = true;
end;

if opts.computeLeverage
    % Compute leverage factors, that is compute the sensitivity of the
    % estimation to each data sample. (0 = not affected, 1 = change in data
    % fully absorbed by change in theta)
    if opts.optimizeBlocks
        h = compute_leverage_factors(theta, data, opts);
    else
        error('cant compute leverage factors for global optimization');
    end;
    theta = h; % parse output arguments.
    return;
end;
% **************************    Finished setting opts; opts is fixed below  **************************    
% **************************    Define currentPoint   **************************    

currentPoint.theta = theta;
currentPoint.LLdata             = []; % total current log likelihood of data term.
currentPoint.LLregularization   = []; % total current log likelihood of regularization term.
currentPoint.LLparameterPrior   = []; % total current log likelihood of prior LL term.
currentPoint.noiseLevel = opts.noiseLevel;
currentPoint.predictedImages = [];    % if non empty, always has current predicted images.
currentPoint.projectParameters = opts.projectParameters;  
currentPoint.projectGrad = opts.projectGrad;


% **************************    Start processing Loop   **************************    
repeatGlobalEstimationlp = 1;
repeatGlobalEstimation = true; % repeatGlobalEstimation  ==  do we need another outer iteration?
while repeatGlobalEstimation % repeatGlobalEstimationlp <= opts.maxOuterloopIterations  % 
    repeatGlobalEstimation = false; % don't repeat, unless we find a reason to do so. 
    
    % Optimize projection parameters (if applicable)
    if isempty(opts.adjProject)
        if isempty(opts.project)
% %             data does not need to be projected
        else
            currentPoint.theta = theta;
            [currentPoint] = optimizeAllProjectParameters(currentPoint, opts);
            currentPoint.theta = []; % save memory; Need to update code that currentPoint is used everywhere.
            
%             if ~opts.projectOptimization.skip || (isempty(opts.projectGrad) && opts.optimizeBlocks) % if projectGrad is provided or not needed and project optimization should be skipped, we trust the projectGrad content and dont optimize.
%             if isempty(currentPoint.predictedImages)
%                 currentPoint.predictedImages = predictAllImagesfun( theta, opts );
%             end;
%             predictedImagesp = permute( currentPoint.predictedImages, [2:numel(szData) 1]); % TODO: can we avoid to permute this large(st) matrix?
%             projectGrad = cell(numel(data_in),1);
%             sel = repmat({':'},1,numel(szData));
%             % determine number of projection steps, for progressbar:
%             cnt = 0;
%             for k = 1 : size(opts.project,1)
%                 if iscell( opts.project{k,1} )
%                     cnt = cnt + numel(opts.project{k,1});
%                 else
%                     cnt = cnt + 1;
%                 end;
%             end;
%             progressbar('start', cnt);
%             done = 0;
%             % optimize all projection's:
%             for k=1:size(opts.project,1)
%                 sel{end} = k;
%                 if size(opts.noiseLevel,1)==1
%                     noiseLevel = opts.noiseLevel;
%                 else
%                     sznoiselvl = size(opts.noiseLevel);
%                     noiseLevel = reshape(opts.noiseLevel(k,:), [1 sznoiselvl(2:end)]);
%                 end;
%                 if size(opts.projectParameterPrior,1)==1
%                     projectParameterPriork = opts.projectParameterPrior;
%                 else
%                     projectParameterPriork = opts.projectParameterPrior(k,:);
%                 end;
%                 if numel(noiseLevel)==1 && iscell(noiseLevel)
%                     noiseLevel = noiseLevel{1};
%                 end;
%                 if iscell( opts.project{k,1} )
%                     nprojk = numel(opts.project{k,1});
%                     projectGrad{k} = struct('projectGrad',{cell(nprojk,1)}, 'logPDFhess',{cell(nprojk,1)},'logPDFgrad',{cell(nprojk,1)});
%                     predictedImagesps = predictedImagesp(sel{:});
%                     for projidx = 1: nprojk
%                         if iscell( opts.logPDFfun )
%                             logPDFfunk = opts.logPDFfun{projidx};
%                         else
%                             logPDFfunk = opts.logPDFfun;
%                         end;
%                         if iscell(noiseLevel)
%                             noiseLevelk = noiseLevel{1}(projidx);
%                         else
%                             noiseLevelk = noiseLevel;
%                         end;
%                         projectk = cell(1,size(opts.project,2));
%                         for idx = 1 : numel(projectk)
%                             if iscell(opts.project{k,idx})
%                                 projectk{idx} = opts.project{k,idx}{projidx};
%                             else
%                                 projectk{idx} = opts.project{k,idx};
%                             end;
%                         end;
%                         if iscell(projectParameterPriork{1})
%                             projectParameterPriorkp = projectParameterPriork{1}(projidx,:);
%                         else
%                             projectParameterPriorkp = projectParameterPriork;
%                         end;
% %                         par_init = opts.projectParameters{k}{projidx};
%                         [tmp1,tmp2,tmp3] = optimizeProjectParameters( data_in{k}{projidx}, predictedImagesps, projectk  , opts.projectParameters{k}{projidx} , logPDFfunk, noiseLevelk, opts.projectOptimization, projectParameterPriorkp);
%                         opts.projectParameters{k}{projidx} = tmp1;
%                         fval{k}(projidx) = tmp2;
%                         projectGrad{k}.projectGrad{projidx} = tmp3.projectGrad;
%                         projectGrad{k}.logPDFgrad{projidx}  = tmp3.logPDFgrad;
%                         projectGrad{k}.logPDFhess{projidx}  = tmp3.logPDFhess;
% %                         projected{k}{projidx}               = tmp3.projected;
% %                         par_upd = opts.projectParameters{k}{projidx} - par_init; % DEBUG
% %                         disp(par_upd) % DEBUG
%                         done = done+1;progressbar(done);
%                     end;
%                 else
%                     [opts.projectParameters{k}, fval(k), projectGrad{k}] = optimizeProjectParameters( data_in{k}, predictedImagesp(sel{:}), opts.project(k,:) , opts.projectParameters{k} , opts.logPDFfun, noiseLevel, opts.projectOptimization, projectParameterPriork);
%                     done = done+1;progressbar(done);
%                 end;
%             end;
%             progressbar('ready');
%             % DEBUG compactify projectGrad:
%             % for k1=1:numel(projectGrad);for k2 = 1:numel(projectGrad{k1}.projectGrad);projectGrad{k1}.projectGrad{k2}.img_out = [];projectGrad{k1}.projectGrad{k2}.pargradient = [];end;end;
%             % save dummy projectGrad
%             opts.projectGrad = projectGrad;
%             clear predictedImagesp
%             end;
            data = currentPoint.predictedImages; % NOTE: different meaning (but same size) for data when 'project' is active. 
                              %       Now it contains the predictedImages values, and projectGrad encodes the distance to the data_in.
        end;
    else
        [currentPoint, data] = optimizeAllAdjProjectParameters(currentPoint, opts, repeatGlobalEstimationlp);
        data = dataAdjFun(data); % conversion for phase or complex data.
    end;
    
    if repeatGlobalEstimationlp==1 || ~isempty(opts.adjProject)
            
    end;

    if (isempty(opts.noiseLevel) || opts.mayRefineNoiseLevel ) && (opts.numPDFoptpar==0)
        if isempty(currentPoint.predictedImages)
            currentPoint.predictedImages = predictAllImagesfun( theta, opts );
        end;
        residue = double(data) - currentPoint.predictedImages;
        % compute noise level from residue:
        resnrm = sum(residue.^2);
        if ~isempty(opts.mask)
            resnrm = resnrm(opts.mask);
        else
        end
        noiselvlMean = sqrt(sum(resnrm(:))/((size(data,1)-(size(theta,1)-opts.numPDFoptpar))*numel(resnrm)+1));
        noiselvl = sqrt(median(resnrm(:))/max(1,size(data,1)-(size(theta,1)-opts.numPDFoptpar))); % use max to avoid division by 0 if number of parameters equal number of measurements. Noise level inaccurate!!!
        opts.noiseLevel = noiselvl;
        disp(['Estimated (median) noise level from residu: ' num2str(noiselvl) ' (not used mean noise level = ' num2str(noiselvlMean) ')']);
    else
        noiselvl = opts.noiseLevel; 
    end;
    
    % now stop if needed:
%     if doStop( curpoint, opts)
%         break;
%     end;
    if repeatGlobalEstimationlp > 1 
        continueIterations = false; 
        if opts.mayRefineNoiseLevel 
            if isempty(oldNoiseEstimate)
                oldNoiseEstimate = noiselvl; % only happens in scalar gaussian, which does not need noise level. 
            end;
            if oldNoiseEstimate >= opts.acceptNoiseLevelAdjustFactors(1) * noiselvl && oldNoiseEstimate <= opts.acceptNoiseLevelAdjustFactors(2) * noiselvl
                % noise level of previous iterate was almost OK, so don't refine results.
            else
                disp('Noise level estimate changed substantially. Therefore the fit is refined.');
                continueIterations = true;
            end;
        end;
        if ~continueIterations && ~isempty(opts.adjProject)
            continueIterations = logLikAfterAlign - logLikBeforeAlign > opts.projectOptimization.outerloopTolerance;
        end;
        if ~continueIterations || repeatGlobalEstimationlp > opts.maxOuterloopIterations
            break;
        end;
    end;        

    if opts.optimizeBlocks
        [theta ,opts] = ML_optimize_blocks(theta, data, opts);
    elseif opts.projectOptimization.skip==-1
        [theta , opts] = ML_optimize_globalThetaAndProject(theta, data, opts);
    else
        [theta ,opts] =  ML_optimize_global(theta, data, opts);
    end;
    currentPoint.predictedImages = []; % theta is updated, so predictedImages is no longer current.
    
%     % now stop if needed:
%     if doStop( curpoint, opts)
%         break;
%     end;

    repeatGlobalEstimationlp = repeatGlobalEstimationlp + 1;
end;
% **************************    end processing loop   **************************    

if nargout>=3
    opts.projectParameters = currentPoint.projectParameters;
    % also output opts.
    % In opts: overwrite the functions that are created in fit_MRI.
    % Reason: when creating these function the workspace get's locked into them
    %         thus when you save opts with such function, the resulting file is huge
    %         and extra memory might be occupied as well.
    opts = parse_defaults_optionvaluepairs( opts, optsOutputAdjust);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Usage types  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Some options interact. Here important interactions are listed:
%
% 1) Project : 
%    - used : then noise level is estimated per data 
% 
% 
% 
% 

% function [stop] = doStop( curpoint, opts)
%     if repeatGlobalEstimationlp > 1 
%         continueIterations = false; 
%         if opts.mayRefineNoiseLevel 
%             if isempty(oldNoiseEstimate)
%                 oldNoiseEstimate = noiselvl; % only happens in scalar gaussian, which does not need noise level. 
%             end;
%             if oldNoiseEstimate >= opts.acceptNoiseLevelAdjustFactors(1) * noiselvl && oldNoiseEstimate <= opts.acceptNoiseLevelAdjustFactors(2) * noiselvl
%                 % noise level of previous iterate was almost OK, so don't refine results.
%             else
%                 disp('Noise level estimate changed substantially. Therefore the fit is refined.');
%                 continueIterations = true;
%             end;
%         end;
%         if ~continueIterations && ~isempty(opts.adjProject)
%             continueIterations = logLikAfterAlign - logLikBeforeAlign > opts.projectOptimization.outerloopTolerance;
%         end;
%         if ~continueIterations || repeatGlobalEstimationlp > opts.maxOuterloopIterations
%             break;
%         end;
%     end;


