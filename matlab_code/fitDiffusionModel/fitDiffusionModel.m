function fitDiffusionModel(varargin)
% fitDiffusionModel fits a diffusion model to diffusion-weighted data.
% A typical fit consists of the following steps:
% 1. Log-linear tensor fit
% 2. (Multiple) initializations for diffusion model (usually computed from log-linear tensor fit).
% 3. Non-linear least squares optimization for each initialization (Levenberg-Marquardt).
% 4. Maximum likelihood estimation with best initialization (fit_MRI).
% 5. Write results to disk.
%
% Note:
% 1. It expects niftitools from Jimmy Shen to be in the matlab-path.
% 2. It expects fit_mri tools from Dirk Poot to be in the matlab-path. It
% may be necessary to compile one or more mex-files in fit_mri.
% 3. fitDiffusionModel was implemented in matlab2013b and has not been
% thoroughly tested on other versions.
%
% Implemented models for diffusion-weighted images acquired at a single b-value:
% 'st_original':        The conventional single tensor model. (L[123]<opts.d_max)
% 'stv2':               A constrained single tensor model. (L2=L3, [AR]D<opts.d_max)
% 'bitensor_mdc':       A trace-constrained bi-tensor model (Trace[tensor] = 3*opts.MD, d_iso = opts.d_iso)
% 'bitensor_adc':       A axial-diffusivity constrained bi-tensor model (AD[tensor] = opts.AD, d_iso = opts.d_iso)
% 'dtv3':               A very constrained dual-tensor model. (AD1=AD2, RD1=RD2, fiso=0, [AR]D<opts.d_max)
% 'dtv4':               A very constrained dual-tensor model. (AD1=AD2, RD1=RD2, fiso=0, f1=0.5, f2=0.5, [AR]D<opts.d_max)
% 'balland1stick':      A very constrained dual-tensor model. (AD=diso, RD=0, [AR]D<opts.d_max)
% 'balland2sticks':     A very constrained dual-tensor model. (AD1=AD2=diso, RD1=RD2=0, [AR]D<opts.d_max)
%
% Implemented models for diffusion-weighted images acquired at multiple b-values:
% 'stv1':				A cigar-shaped bi-tensor model (L2=L3, d_iso = opts.d_iso)
% 'dtv1':				A loosely constrained dual-tensor model. (AD1=AD2, diso=opts.d_iso, [AR]D<opts.d_max)
% 'dtv2':				A loosely constrained dual-tensor model. (AD1=AD2, RD1=RD2, diso=opts.d_iso, [AR]D<opts.d_max)
% 'bitensor':           A bi-tensor model. (L[123]<opts.d_max)
%
% Implemented special models:
% 'balland2sticks_orientationprior':       A balland2sticks model that uses an orientation prior to stabilize the fit.
% 
% The input is an opts structure with at least the fields <data>, <mask>,
% <bvals>, <bvecs>, <out> and <model> specified.
%
% It can (quite easily) be compiled with the following matlab command line command: 
% mcc -R -nosplash -R -singleCompThread  -m -v fitDiffusionModel.m 
% And then it can be run from the command line (it may be required to set some library paths):
% fitDiffusionModel -data <fn_data> -mask <fn_mask> -bvals <fn_bvals> -bvecs <fn_bvecs> -out <dir_out> -model <model>
%
% created by Joor Arkesteijn (joorarkesteijn@gmail.com), Quantitative Imaging Group, Delft.
%
% 8-9-2014:       v1.0 Main changes:
%               - changed name to fitDiffusionModel (instead of more cumbersome fitDiffusionTensorToNifti)
%               - lsq estimation and ml estimation are now two different loops such the global lsq-residual can be used to compute a global noiselvl for the ml-estimation.
%               - added diffusion models 'stv1' and 'dtv2'
%               - added easy options to control memory use and show fitmri progressbar
%               - more options to estimate noiselevel including spatial regularization of the noiselvl
%               - the fitmri blocksize and blockoverlap can be set from the opts-structur
%               - all models now use two initializations to increase robustness
%               - added heuristic rules and weights to log-linear tensor fit to make it more robust against low-intensity outliers and multi-shell diffusion-weighted images
%               - decreased memory-use of all theta2nii functions
%               - no longer using the L11, L21, L31, L12, L22, L23 notation for the dual-tensor model, but switched to RD1, RD2, AD1 and AD2 to save disk space (and it makes more sense)
%
% 10-10-2014:     v1.1 Changes
%               - small bug fix in log file. .. slices out of X slices in total
%               - small bug fix in log file. .. No voxels in mask in slices
%               - added option to apodize dwis
%
% 30-03-2018:     v1.2 Changes
%               - added diffusion model 'stv2'
%               - fold periodical parameters (e.g. alphaX, thetaX, phiX) back to the interval 0 to 2*pi. Otherwise accuracy may be lost when saving the parameters as floats.
%               - added function theta2nii that can be used to recompute derived files
%               - added blurring over slice direction in optional apodization step (opts.apodization_sigmaz)
%               - added diffusion model 'bitensor'
%               - added a preliminary option to use a prior
%               - added diffusion model 'balland2sticks_orientationprior'
%
%               - TODO: strsplit is not defined in Matlab2011?
%               - TODO: Tools from Dirk doesn't work as expected in Matlab2014 (something about hessian that was set, but not specified)

%% Required options
opts.data = ''; %% Opts.data should contain a string with the file name to a nifti-archive (.nii) or compressed nifti-archive (nii.gz) containing the diffusion-weighted images.
opts.mask = ''; %% Opts.mask should contain a string with the file name to a nifti-archive (.nii) or compressed nifti-archive (nii.gz) containing a binary mask (1 = brain, 0 = no brain).
opts.bvals = ''; %% Opts.bvals should contain a string with the file name of a text-file with b-values.
opts.bvecs = ''; %% Opts.bvecs should contain a string with the file name of a text-file with b-vectors.
opts.out = ''; %% Opts.out should contain a string with the directory where the output should be written. If the location does not exist it will be created.
opts.model = ''; % Opts.model should contain a string specifying the model to be fit. Various options are possible.

%% Optional options
opts.lsq = 'true'; % The non-linear least squares optimization (step 3) can be skipped by setting opts.lsq to 'false'. 
opts.ml = 'true'; % Maximum-likelihood estimation (step 4) can be skipped by setting opts.ml to 'false'. 
opts.lsq_mapping = 'sin'; %% Select lsq-mapping for constrained parameters: 'none', 'erf' or 'sin'.
opts.ml_mapping = 'constrain'; %% Select ml-mapping for constrained parameters: 'none', 'erf', 'sin' or 'constrain'.
opts.crlb = 'false'; %% Computes the CRLB in each voxel and writes this to the disk.
opts.turnofwarning = 'true'; %% Turn of matrix nearly singular warning.
opts.fitsingleslice = 'false'; %% For testing purposes only a single slice in the center is fitted.
opts.d_max = '3e-3'; %% Maximum diffusivity of tensor compartment.
opts.d_iso = '3e-3'; %% Diffusivity of isotropic compartment.
opts.MD = '0.8e-3'; %% Mean diffusivity of the tensor.
opts.AD = '1.6e-3'; %% Axial diffusivity of the tensor.
opts.initializationfile = ''; %% Specify a file to use for the initialization, instead of initializing the model from a log-linear fit.
opts.fieldsfile = ''; %% Specify a file to use for the fields. This could be required for some diffusion-models.
opts.checkjacobian = 'false'; %% Output the analytically and numerically computed Jacobian for visual inspection.
opts.inspectoutliers = 'false'; %% Show a plot that can be used to detect outlier DWIs based on a log-linear tensor fit. Only works well for single b-values dwis.
opts.globalnoiselvl = 'lsqmodel'; %% Specify how the global noiselvl will be computed, choose between 'logtensor' or 'lsqmodel' or specify a value.
opts.refinenoiselvl = 'false'; %% Refine noiselvl with ML-estimation (or use the global noiselvl in each voxel).
opts.regularizationweight = '0'; %% Use spatial regularization parameters when refining the noiselvl by setting opt.regularizationweight to non-zero value, e.g. 1e6;
opts.memorylimit = '50'; %% Memorylimit determines the number of slices (in MegaBytes) loaded simultaneously. For effective spatial regularization this needs to be set to high value to prevent edges, e.g. 2000 (but this does increase memory consumption of course).
opts.fitmriprogressbar = 'false'; %% Show the fitmri progress bar by setting opts.fitmriprogressbar to 'true'.
opts.fitmriblocksize = '8 8 8'; %% The x, y, and z-dimensions seperated by white spaces of the blocksize that will be used for maximum likelihood estimation with fitmri.
opts.fitmriblockoverlap = '4 4 4'; %% The x, y, and z-dimensions seperated by white spaces of the blockoverlap that will used for maximum likelihood estimation with fitmri (but only when spatial regularization of the noiselvl is on)
opts.apodize_window = 'none'; %% Apply an apodization window to all DWIs to remove Gibbs-ringing. Choose between 'lukosz', 'boxcar', 'triangular', 'gauss' and 'staircase'
opts.apodize_kx = ''; %% Set width of apodization window for x-direction. 
opts.apodize_ky = ''; %% Set width of apodization window for y-direction.
opts.apodize_sigmaz = '0'; %% Set sigma for blurring in z-direction.
opts.number_priors = 0; %% Don't change this. Default estimation without prior.
opts.orientation_prior = '25'; %% This determines the width of the prior for the orientations

%% If deployed, remove the minus signs in front of the option fields (odd fields)
if isdeployed
    for i=1:2:numel(varargin)
        if strcmpi(varargin{i}(1),'-')
            varargin{i}=varargin{i}(2:end);
        end
    end
end

%% Parse optionvaluepairs
opts = parse_defaults_optionvaluepairs( opts, varargin{:});

%% Convert chars to doubles or booleans
opts.lsq = char2bool(opts.lsq);
opts.ml = char2bool(opts.ml);
opts.crlb = char2bool(opts.crlb);
opts.turnofwarning = char2bool(opts.turnofwarning);
opts.fitsingleslice = char2bool(opts.fitsingleslice);
opts.d_max = str2double(opts.d_max);
opts.d_iso = str2double(opts.d_iso);
opts.MD = str2double(opts.MD);
opts.AD = str2double(opts.AD);
opts.checkjacobian = char2bool(opts.checkjacobian);
opts.inspectoutliers = char2bool(opts.inspectoutliers);
opts.refinenoiselvl = char2bool(opts.refinenoiselvl);
opts.regularizationweight = str2double(opts.regularizationweight);
opts.fitmriprogressbar = char2bool(opts.fitmriprogressbar);
opts.memorylimit = str2double(opts.memorylimit);
opts.fitmriblocksize = str2double(strsplit(opts.fitmriblocksize));
opts.fitmriblockoverlap = str2double(strsplit(opts.fitmriblockoverlap));
opts.orientation_prior = str2double(opts.orientation_prior);
    
%% Assign useful names
fn_dwi = opts.data;
fn_mask = opts.mask;
fn_bval = opts.bvals;
fn_bvec = opts.bvecs;
dir_output = opts.out;

%% Create a log file
global fid;
if dir_output(end)~=filesep()
    dir_output = [dir_output filesep()];
end
if ~exist(dir_output,'dir')
    mkdir(dir_output)
end
fid = fopen([dir_output 'log.txt'],'w');

%% Print and log the estimation options
displayandlog('Using the following options for the optimization:')
displayandlog('==============================================================================')
print_opts(opts);
displayandlog('==============================================================================')

%% load b-values and gradients, and turn into column vectors
b = load(fn_bval);
if size(b,2)>size(b,1)
    b = b';
end
g = load(fn_bvec);
if size(g,2)>size(g,1)
    g = g';
end

%% Compute q and Q
q = g.*repmat(sqrt(b),1,3);
Q = zeros(length(q),6);
Q(:,1) = -q(:,1).*q(:,1);
Q(:,2) = -q(:,2).*q(:,2);
Q(:,3) = -q(:,3).*q(:,3);
Q(:,4) = -2*q(:,1).*q(:,2);
Q(:,5) = -2*q(:,1).*q(:,3);
Q(:,6) = -2*q(:,2).*q(:,3);

%% Unzip nifti-file if compressed, and delete it afterwards
delete_dwi = false;
if strcmpi('.gz',fn_dwi(end-2:end))
    displayandlog('Unzipping compressed nifti file...');
    fn_dwi = gunzip(fn_dwi,dir_output);
    fn_dwi = fn_dwi{1};
    delete_dwi = true;
end

%% Apply an apodization window if specified
if ~strcmpi(opts.apodize_window,'none')
    
    displayandlog('Applying apodization window to DWIs...')
    
    if delete_dwi % No need to preserve original dwis
        opts_apo.in = fn_dwi;
        opts_apo.out = fn_dwi;
    else % Place apodized dwis in output dir and delete afterwards
        opts_apo.in = fn_dwi;
        opts_apo.out = [dir_output 'apodized_dwis.nii.gz'];
        delete_dwi = true;
    end
    
    opts_apo.kx = opts.apodize_kx;
    opts_apo.ky = opts.apodize_ky;
    opts_apo.sigmaz = opts.apodize_sigmaz;
    opts_apo.w = opts.apodize_window;
    apodize(opts_apo);
    
end

%% Switch annoying warning off
if opts.turnofwarning
    warning('off','MATLAB:nearlySingularMatrix');
end

%% Load the header, and mask
img_header = load_untouch_header_only(fn_dwi);
img_dim = img_header.dime.dim(2:5);
mask_all = load_untouch_nii(fn_mask);
mask_all = mask_all.img>0;

%% Maybe for testing purposes you just want a subset of the mask?
if opts.fitsingleslice
    displayandlog('Only using subset of mask!')
    mask_all(:,:,[1:floor(end/2)-1 floor(end/2)+1:end]) = false;
end

%% Estimate (initial) noise level
[globalnoiselvl,globalsignallvl,noiselvl_perdwi,~,S0_all,DT_all] = checkDWI(fn_dwi,fn_mask,Q);
noiselvl_all = reshape(mask_all*globalnoiselvl,[1 size(mask_all)]);
displayandlog(['Voxels removed from mask because of insufficient SNR (S0<globalnoiselvl): ' num2str(sum(S0_all(:)<globalnoiselvl & mask_all(:))) '...']);
mask_all(S0_all<globalnoiselvl) = false;
displayandlog(['Total number of voxels in specified mask: ' num2str(sum(mask_all(:))) '...']);
displayandlog(['Global noise level computed from residuals logtensor fit (median all voxels): ' num2str(globalnoiselvl) '...']);
displayandlog(['Signal level (median all voxels): ' num2str(globalsignallvl) '...']);
displayandlog(['SNR (median all voxels): ' num2str(globalsignallvl/globalnoiselvl) '...']);
if opts.inspectoutliers
    lscatter(b(:)+100*rand(numel(b),1),noiselvl_perdwi(:),cellstr(num2str([1:numel(b)]')));
    xlabel('b-values of diffusion-weighted image');
    ylabel('residual of log linear single tensor fit averaged over all voxels in diffusion-weighted image');
    title('Please check this plot for outliers. Diffusion-weighted images acquired at similar b-values should have similar averaged residuals')
end

%% Do not load more than <opts.memorylimit> MB of uncompressed data simultaneously
memory_limit = opts.memorylimit*1024*1024*8; %MB -> bits
memory_singleslice = img_dim(1)*img_dim(2)*img_dim(4)*32; %bits
slices_per_iter = ceil(memory_limit/memory_singleslice);

%% Specify the model
if strcmpi(opts.model,'st_original')
    
    %% Define useful functions
    displayandlog('Fitting with diffusion model st_original...');
    fun = @(theta,~) st_original_fitmri(theta,Q,opts);
    model_theta0 = @(S0,DT,~,~,~,mask) st_original_theta0(S0,DT,mask,opts);
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) st_original_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
    %% Specify parameters that should be mapped or constrained between 0 and 1
    map_par = [false;false;false;true;true;true;false];
    
    %% Specify parameters that have period 2*pi
    period2pi_par = [true;true;true;false;false;false;false];
    npar = numel(map_par);
    pos_S0 = npar;
    
elseif strcmpi(opts.model,'stv2')
    
    %% Define useful functions
    displayandlog('Fitting with diffusion model stv2...');
    fun = @(theta,~) stv2_fitmri(theta,Q,opts);
    model_theta0 = @(S0,DT,~,~,~,mask) stv2_theta0(S0,DT,mask,opts);
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) stv2_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
    %% Specify parameters that should be mapped or constrained between 0 and 1
    map_par = [true;true;false;false;false];
    
    %% Specify parameters that have period 2*pi
    period2pi_par = [false;false;true;true;false];
    npar = numel(map_par);
    pos_S0 = npar;
    
elseif strcmpi(opts.model,'bitensor')
    
    %% Define useful functions
    displayandlog('Fitting with diffusion model bitensor...');
    fun = @(theta,~) bitensor_fitmri(theta,Q,opts);
    model_theta0 = @(S0,DT,~,~,~,mask) bitensor_theta0(S0,DT,mask,opts);
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) bitensor_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
    %% Specify parameters that should be mapped or constrained between 0 and 1
    map_par = [false;false;false;true;true;true;true;false];
    
    %% Specify parameters that have period 2*pi
    period2pi_par = [true;true;true;false;false;false;false;false];
    npar = numel(map_par);
    pos_S0 = npar;
    
elseif strcmpi(opts.model,'bitensor_mdc')
    
    %% Define useful functions
    displayandlog('Fitting with diffusion model bitensor_mdc...');
    fun = @(theta,~) bitensor_mdc_fitmri(theta,Q,opts);
    model_theta0 = @(S0,DT,~,~,~,mask) bitensor_mdc_theta0(S0,DT,mask,opts);
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) bitensor_mdc_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
    %% Specify parameters that should be mapped or constrained between 0 and 1
    map_par = [false;false;false;true;true;false;false];
    
    %% Specify parameters that have period 2*pi
    period2pi_par = [true;true;true;false;false;false;false];
    npar = numel(map_par);
    pos_S0 = npar;
    
elseif strcmpi(opts.model,'bitensor_adc')
    
    %% Define useful functions
    displayandlog('Fitting with diffusion model bitensor_adc...');
    fun = @(theta,~) bitensor_adc_fitmri(theta,Q,opts);
    model_theta0 = @(S0,DT,~,~,~,mask) bitensor_adc_theta0(S0,DT,mask,opts);
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) bitensor_adc_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
    %% Specify parameters that should be mapped or constrained between 0 and 1
    map_par = [false;false;false;true;true;false;false];
    
    %% Specify parameters that have period 2*pi
    period2pi_par = [true;true;true;false;false;false;false];
    npar = numel(map_par);
    pos_S0 = npar;

elseif strcmpi(opts.model,'stv1')
    
    %% Define useful functions
    displayandlog('Fitting with diffusion model stv1...');
    fun = @(theta,~) stv1_fitmri(theta,Q,opts);
    model_theta0 = @(S0,DT,~,~,~,mask) stv1_theta0(S0,DT,mask,opts);
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) stv1_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
    %% Specify parameters that should be mapped or constrained between 0 and 1
    map_par = [true;true;true;false;false;false];
    
    %% Specify parameters that have period 2*pi
    period2pi_par = [false;false;false;true;true;false];
    npar = numel(map_par);
    pos_S0 = npar;
    
elseif strcmpi(opts.model,'dtv1')
    
    %% Define useful functions
    displayandlog('Fitting with diffusion model dtv1...');
    fun = @(theta,~) dtv1_fitmri(theta,Q,opts);
    model_theta0 = @(S0,DT,~,~,~,mask) dtv1_theta0(S0,DT,mask,opts);
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) dtv1_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
    %% Specify parameters that should be mapped or constrained between 0 and 1
    map_par = [true;true;true;true;true;false;false;false;false;false];
    
    %% Specify parameters that have period 2*pi
    period2pi_par = [false;false;false;false;false;true;true;true;true;false];
    npar = numel(map_par);
    pos_S0 = npar;
    
elseif strcmpi(opts.model,'dtv2')
    
    %% Define useful functions
    displayandlog('Fitting with diffusion model dtv2...');
    fun = @(theta,~) dtv2_fitmri(theta,Q,opts);
    model_theta0 = @(S0,DT,~,~,~,mask) dtv2_theta0(S0,DT,mask,opts);
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) dtv2_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
    %% Specify parameters that should be mapped or constrained between 0 and 1
    map_par = [true;true;true;true;false;false;false;false;false];
    
    %% Specify parameters that have period 2*pi
    period2pi_par = [false;false;false;false;true;true;true;true;false];
    npar = numel(map_par);
    pos_S0 = npar;
    
elseif strcmpi(opts.model,'dtv3')
    
    %% Define useful functions
    displayandlog('Fitting with diffusion model dtv3...');
    fun = @(theta,~) dtv3_fitmri(theta,Q,opts);
    model_theta0 = @(S0,DT,~,~,~,mask) dtv3_theta0(S0,DT,mask,opts);
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) dtv3_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
    %% Specify parameters that should be mapped or constrained between 0 and 1
    map_par = [true;true;true;false;false;false;false;false];
    
    %% Specify parameters that have period 2*pi
    period2pi_par = [false;false;false;true;true;true;true;false];
    npar = numel(map_par);
    pos_S0 = npar;
    
elseif strcmpi(opts.model,'dtv4')
    
    %% Define useful functions
    displayandlog('Fitting with diffusion model dtv4...');
    fun = @(theta,~) dtv4_fitmri(theta,Q,opts);
    model_theta0 = @(S0,DT,~,~,~,mask) dtv4_theta0(S0,DT,mask,opts);
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) dtv4_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
    %% Specify parameters that should be mapped or constrained between 0 and 1
    map_par = [true;true;false;false;false;false;false];
    
    %% Specify parameters that have period 2*pi
    period2pi_par = [false;false;true;true;true;true;false];
    npar = numel(map_par);
    pos_S0 = npar;
    
elseif strcmpi(opts.model,'dtv3_constrainedangles')
    
    %% Define useful functions
    displayandlog('Fitting with diffusion model dtv3_constrainedangles...');
    fun = @(theta,fields) dtv3_constrainedangles_fitmri(theta,fields,Q,opts);
    model_theta0 = @(S0,DT,~,~,~,mask) dtv3_constrainedangles_theta0(S0,DT,mask,opts);
    model_theta2nii = @(theta_all,fields_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) dtv3_constrainedangles_theta2nii(theta_all,fields_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
    %% Specify parameters that should be mapped or constrained between 0 and 1
    map_par = [true;true;true;false];
    
    %% Specify parameters that have period 2*pi
    period2pi_par = [false;false;false;false];
    npar = numel(map_par);
    pos_S0 = npar;
    
elseif strcmpi(opts.model,'balland1stick')
    
    %% Define useful functions
    displayandlog('Fitting with diffusion model balland1stick...');
    fun = @(theta,~) balland1stick_fitmri(theta,Q,opts);
    model_theta0 = @(S0,DT,data,Q,~,mask) balland1stick_theta0(S0,DT,data,Q,mask,opts);
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) balland1stick_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
    %% Specify parameters that should be mapped or constrained between 0 and 1
    map_par = [true;true;false;false;false];
    
    %% Specify parameters that have period 2*pi
    period2pi_par = [false;false;true;true;false];
    npar = numel(map_par);
    pos_S0 = npar;
    
elseif strcmpi(opts.model,'balland2sticks')
    
    %% Define useful functions
    displayandlog('Fitting with diffusion model balland2sticks...');
    fun = @(theta,~) balland2sticks_fitmri(theta,Q,opts);
    model_theta0 = @(S0,DT,data,Q,~,mask) balland2sticks_theta0(S0,DT,data,Q,mask,opts);
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) balland2sticks_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
    %% Specify parameters that should be mapped or constrained between 0 and 1
    map_par = [true;true;true;false;false;false;false;false];
    
    %% Specify parameters that have period 2*pi
    period2pi_par = [false;false;false;true;true;true;true;false];
    npar = numel(map_par);
    pos_S0 = npar;
    
elseif strcmpi(opts.model,'balland2sticks_orientationprior')
    
    %% The opts structure is used to specify the Gaussian priors
    opts.number_priors = 2; %% Two priors, one for each orientation
    opts.prior_offset = [10; 10]; %% Use offset to prevent Rician bias
    
    %% Define useful functions
    displayandlog('Fitting with diffusion model balland2sticks...');
    fun = @(theta,fields) balland2sticks_orientationprior_fitmri(theta,fields,Q,opts);
    model_theta0 = @(S0,DT,data,Q,fields,mask) balland2sticks_orientationprior_theta0(S0,DT,data,Q,fields,mask,opts);
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) balland2sticks_orientationprior_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);

    %% Specify parameters that should be mapped or constrained between 0 and 1
    map_par = [true;true;true;false;false;false;false;false];
    
    %% Specify parameters that have period 2*pi
    period2pi_par = [false;false;false;true;true;true;true;false];
    npar = numel(map_par);
    pos_S0 = npar;
    
elseif strcmpi(opts.model,'balland2sticks_2timepoints')
    
    %% The opts structure is used to specify the Gaussian priors
    opts.number_priors = 2; %% Two priors, one for each orientation
    opts.prior_offset = [10; 10]; %% Use offset to prevent Rician bias
    
    %% Define useful functions
    displayandlog('Fitting with diffusion model balland2sticks_2timepoints...');
    fun = @(theta,fields) balland2sticks_2timepoints_fitmri(theta,fields,Q,opts);
    model_theta0 = @(S0,DT,data,Q,fields,mask) balland2sticks_2timepoints_theta0(S0,DT,data,Q,fields,mask,opts);
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) balland2sticks_2timepoints_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
    %% Specify parameters that should be mapped or constrained between 0 and 1
    map_par = [true;true;true;false;true;true;true;false;false;false;false;false];
    
    %% Specify parameters that have period 2*pi
    period2pi_par = [false;false;false;false;false;false;false;false;true;true;true;true];
    npar = numel(map_par);
    pos_S0 = [4 8];
    
elseif strcmpi(opts.model,'balland2sticks_3timepoints')
    
    %% The opts structure is used to specify the Gaussian priors
    opts.number_priors = 2; %% Two priors, one for each orientation
    opts.prior_offset = [10; 10]; %% Use offset to prevent Rician bias
    
    %% Define useful functions
    displayandlog('Fitting with diffusion model balland2sticks_3timepoints...');
    fun = @(theta,fields) balland2sticks_3timepoints_fitmri(theta,fields,Q,opts);
    model_theta0 = @(S0,DT,data,Q,fields,mask) balland2sticks_3timepoints_theta0(S0,DT,data,Q,fields,mask,opts);
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) balland2sticks_3timepoints_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
    %% Specify parameters that should be mapped or constrained between 0 and 1
    map_par = [true;true;true;false;true;true;true;false;true;true;true;false;false;false;false;false];
    
    %% Specify parameters that have period 2*pi
    period2pi_par = [false;false;false;false;false;false;false;false;false;false;false;true;true;true;true];
    npar = numel(map_par);
    pos_S0 = [4 8 12];
    
elseif strcmpi(opts.model,'balland2sticks_3timepoints_noprior')
    
    %% Define useful functions
    displayandlog('Fitting with diffusion model balland2sticks_3timepoints_noprior...');
    fun = @(theta,fields) balland2sticks_3timepoints_noprior_fitmri(theta,fields,Q,opts);
    model_theta0 = @(S0,DT,data,Q,fields,mask) balland2sticks_3timepoints_noprior_theta0(S0,DT,data,Q,mask,opts);
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) balland2sticks_3timepoints_noprior_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
    %% Specify parameters that should be mapped or constrained between 0 and 1
    map_par = [true;true;true;false;true;true;true;false;true;true;true;false;false;false;false;false];
    
    %% Specify parameters that have period 2*pi
    period2pi_par = [false;false;false;false;false;false;false;false;false;false;false;true;true;true;true];
    npar = numel(map_par);
    pos_S0 = [4 8 12];

else
    
    error('No valid model has been selected');
    
end

%% Load fields
if isempty(opts.fieldsfile)
    fields_all = zeros(1,img_dim(1),img_dim(2),img_dim(3));
else
    displayandlog('Loading fields from specified file...');
    fields_all = load_untouch_nii(opts.fieldsfile);
    fields_all = double(permute(fields_all.img,[4 1 2 3]));
end

%% Pre-allocate matrices
theta_all = zeros(npar,img_dim(1),img_dim(2),img_dim(3));
resnorm_all = zeros(1,img_dim(1),img_dim(2),img_dim(3));
logLL_all = zeros(1,img_dim(1),img_dim(2),img_dim(3));

%% Iteration over all slices for initialization and/or least squares fit
displayandlog('Iterating over all slices for initialization and/or voxelwise least squares fit...');

%% Initializations for ml loop over all slices
start_idx = 1;
slices_remaining = true;

while slices_remaining
    
    %% Select slices
    end_idx = start_idx+slices_per_iter-1;
    if end_idx>=img_dim(3)
        end_idx = img_dim(3);
        slices_remaining = false;
    end
    slice_idx = start_idx:end_idx;
    
    %% Load slices
    displayandlog(['Loading slices ' num2str(start_idx) ' to ' num2str(end_idx) ' of ' num2str(size(mask_all,3)) ' slices in total...']);
    data = load_slices(fn_dwi,slice_idx);
    mask = mask_all(:,:,slice_idx);
    noiselvl = noiselvl_all(:,:,:,slice_idx);
    S0 = S0_all(:,:,:,slice_idx);
    DT = DT_all(:,:,:,slice_idx);
    fields = fields_all(:,:,:,slice_idx);
    displayandlog(['Slices ' num2str(start_idx) ' to ' num2str(end_idx) ' contain ' num2str(sum(mask(:))) ' voxels...']);
    
    if sum(mask(:))>0
        
        %% Scale the data for better conditioned optimization
        data = double(data)/globalsignallvl;
        noiselvl = noiselvl/globalsignallvl;
        S0 = S0/globalsignallvl;
        
        %% Compute initial points
        if isempty(opts.initializationfile)
            displayandlog('Computing initial points...');
            theta_cell = model_theta0(S0,DT,data,Q,fields,mask);
        else
            displayandlog('Loading initial points from specified file...');
            theta0 = load_untouch_nii(opts.initializationfile,[],[],[],[],[],slice_idx);
            theta_cell{1} = permute(theta0.img,[4 1 2 3]);
            theta_cell{1}(pos_S0,:) = theta_cell{1}(pos_S0,:)/globalsignallvl;
        end
        
        %% Add priors as additional datapoints
        if opts.number_priors>0
           data = cat(1,data,repmat(opts.prior_offset,[1 size(data,2) size(data,3) size(data,4)]));
        end

        %% Optimize all initial points
        for i=1:numel(theta_cell)
            
            theta_temp = theta_cell{i};
        
            if opts.lsq
                
                %% Put the mapped parameters in the range from [eps to 1-eps] to prevent infs when taking asin/erfinv
                theta_temp(map_par,:) = min(max(theta_temp(map_par,:),eps),1-eps);

                %% Apply mapping
                if strcmpi(opts.lsq_mapping,'erf')
                    %displayandlog('Applying erf mapping for non-linear least squares (Levenberg-Marquardt)...');
                    theta_temp(map_par,:) = erfinv(2*theta_temp(map_par,:)-1);
                    lsq_fun = @(theta,fields) erfmap_fun(theta,fields,fun,map_par);
                elseif strcmpi(opts.lsq_mapping,'sin')
                    %displayandlog('Applying sin mapping for non-linear least squares (Levenberg-Marquardt)...');                    
                    theta_temp(map_par,:) = asin(2*theta_temp(map_par,:)-1);
                    lsq_fun = @(theta,fields) sinmap_fun(theta,fields,fun,map_par);
                elseif strcmpi(opts.lsq_mapping,'none')
                    %displayandlog('Applying no mapping for non-linear least squares (Levenberg-Marquardt)...');
                    lsq_fun = @(theta,fields) fun(theta,fields);
                else
                    error('Invalid mapping selected')
                end
                
                %% Check Jacobian
                if opts.checkjacobian
                    check_Jacobian(lsq_fun,theta_temp,fields,mask);
                    input('Take some time to check the Jacobian :)');
                end

                %% Fit the diffusion data with levenberg_marquardt
                displayandlog('Fitting with voxel-wise non-linear least squares (Levenberg-Marquardt)...');
                theta_temp = levenberg_marquardt_fast(lsq_fun,data,theta_temp,fields,mask);

                %% Undo mapping
                if strcmpi(opts.lsq_mapping,'erf')
                    %displayandlog('Undoing erf mapping for non-linear least squares (Levenberg-Marquardt)...');
                    theta_temp(map_par,:) = erf(theta_temp(map_par,:))/2+0.5;
                elseif strcmpi(opts.lsq_mapping,'sin')
                    %displayandlog('Undoing sin mapping for non-linear least squares (Levenberg-Marquardt)...');
                    theta_temp(map_par,:) = sin(theta_temp(map_par,:))/2+0.5;
                end
            end
            
            theta_cell{i} = theta_temp;
        end
        
        
        %% Select fit with smallest residual
        theta = theta_cell{1};
        for i=2:numel(theta_cell)           
            theta = select_best_fit(fun,data,theta,theta_cell{i},fields,mask,noiselvl,noiselvl,'lsq',['Parameter vector from fit ' num2str(i)]);
        end
        
        %% Scale the data and estimated parameters back
        data(1:(end-opts.number_priors),:,:,:) = data(1:(end-opts.number_priors),:,:,:)*globalsignallvl;
        noiselvl = noiselvl*globalsignallvl;
        theta(pos_S0,:,:,:) = theta(pos_S0,:,:,:)*globalsignallvl;
               
        %% Compute some statistics regarding the fit
        [resnorm,logLL] = compute_fit_statistics(fun,data,theta,fields,mask,noiselvl,opts);
        
        %% Store results in matrices        
        theta_all(:,:,:,slice_idx) = theta;
        mask_all(:,:,slice_idx) = mask;
        resnorm_all(:,:,:,slice_idx) = resnorm;
        logLL_all(:,:,:,slice_idx) = logLL;
        
    else
        displayandlog(['No voxels in mask in slice ' num2str(start_idx) ' to ' num2str(end_idx) '...']);
    end
    
    %% Set new start index
    start_idx = end_idx+1;
    
end

%% Compute (initial) globalnoiselvl for maximum likelihood fit
if strcmpi(opts.globalnoiselvl,'logtensor') % Compute global noiselvl from residual of logtensor fit
    noiselvl_all(:) = globalnoiselvl;
    displayandlog(['Global noiselvl computed from residual of logtensor fit: ' num2str(globalnoiselvl) '...']);
elseif strcmpi(opts.globalnoiselvl,'lsqmodel') % Compute global noiselvl from residual of model lsq-fit
    globalnoiselvl = median(sqrt(resnorm_all(mask_all)/(size(data,1)-npar)));
    noiselvl_all(:) = globalnoiselvl;
    displayandlog(['Global noiselvl computed from residual of lsq ' opts.model '-fit: ' num2str(globalnoiselvl) '...']);
else % Use specified global noiselvl
    globalnoiselvl = str2double(opts.globalnoiselvl);
    noiselvl_all(:) = globalnoiselvl;
    if isnan(globalnoiselvl)
        error('opts.globalnoiselvl should be one of the following options <logtensor>, <lsqmodel> or string that can be converted to a double specifying noiselvl'); 
    end
    displayandlog(['Global noiselvl was specified in opts structure: ' num2str(globalnoiselvl) '...']);
end
displayandlog(['Recomputed SNR (median all voxels): ' num2str(globalsignallvl/globalnoiselvl) '...']);

if opts.ml
    
    %% Pre-allocate matrices
    resnorm_all = zeros(1,img_dim(1),img_dim(2),img_dim(3));
    logLL_all = zeros(1,img_dim(1),img_dim(2),img_dim(3));

    %% Initializations for ml loop over all slices
    start_idx = 1;
    slices_remaining = true;
    
    %% Iteration over all slices for initialization and/or least squares fit
    displayandlog('Iterating over all slices for maximum likelihood fit...');

    while slices_remaining

        %% Select slices
        end_idx = start_idx+slices_per_iter-1;
        if end_idx>=img_dim(3)
            end_idx = img_dim(3);
            slices_remaining = false;
        end
        slice_idx = start_idx:end_idx;

        %% Load slices
        displayandlog(['Loading slices ' num2str(start_idx) ' to ' num2str(end_idx) ' of ' num2str(size(mask_all,3)) ' slices in total...']);
        data = load_slices(fn_dwi,slice_idx);
        mask = mask_all(:,:,slice_idx);
        noiselvl = noiselvl_all(:,:,:,slice_idx);
        fields = fields_all(:,:,:,slice_idx);
        theta = theta_all(:,:,:,slice_idx);
        displayandlog(['Slices ' num2str(start_idx) ' to ' num2str(end_idx) ' contain ' num2str(sum(mask(:))) ' voxels...']);

        if sum(mask(:))>0
            
            %% Scale the data for better conditioned optimization
            data = double(data)/globalsignallvl;
            noiselvl = noiselvl/globalsignallvl;
            theta(pos_S0,:,:,:) = theta(pos_S0,:,:,:)/globalsignallvl;
            
            %% Add priors as additional datapoints
            if opts.number_priors>0
                data = cat(1,data,repmat(opts.prior_offset,[1 size(data,2) size(data,3) size(data,4)]));
            end
            
            %% Put the mapped parameters in the range from [eps to 1-eps] to prevent infs when taking asin/erfinv
            theta(map_par,:) = min(max(theta(map_par,:),eps),1-eps);
            
            %% Apply mapping
            if strcmpi(opts.ml_mapping,'erf')
                %displayandlog('Applying erf mapping for ML estimator...');  
                theta(map_par,:) = erfinv(2*theta(map_par,:)-1);
                ml_fun = @(theta,fields) erfmap_fun(theta,fields,fun,map_par);
                constraints = [];
            elseif strcmpi(opts.ml_mapping,'sin')
                %displayandlog('Applying sin mapping for ML estimator...');
                theta(map_par,:) = asin(2*theta(map_par,:)-1);
                ml_fun = @(theta,fields) sinmap_fun(theta,fields,fun,map_par);
                constraints = [];
            elseif strcmpi(opts.ml_mapping,'none')
                %displayandlog('Applying no mapping for ML estimator...');
                ml_fun = @(theta,fields) fun(theta,fields);
                constraints = [];
            elseif strcmpi(opts.ml_mapping,'constrain')
                %displayandlog('Using constrained optimization for ML estimator...');
                ml_fun = @(theta,fields) fun(theta,fields);
                
                %% Specify the constraints
                lb = -inf(npar,1);
                ub = inf(npar,1);
                lb(map_par) = 0;
                ub(map_par) = 1;
                
                if opts.refinenoiselvl
                    lb = [lb;0];
                    ub = [ub;inf];
                    if opts.regularizationweight~=0
                        error('Constrained optimization and regularization cannot be used simultaneously. Optimizer 1 in fit_MRI supports constrained optimization, but no regularization. Optimizer 3 in fit_MRI supports regularization, but no constrained optimization. Use sin-mapping or erf-mapping instead.');
                    end
                end
                    
                %% Set the constraints
                constraints.lb = lb;
                constraints.ub = ub;
                constraints.Aineq = [];
                constraints.bineq = [];
                constraints.Aeq = [];
                constraints.beq = [];
                constraints.nonlcon = [];
            else
                error('Invalid mapping selected')
            end
            
            %% Check Jacobian
            if opts.checkjacobian
                check_Jacobian(ml_fun,theta,fields,mask);
                input('Take some time to check the Jacobian :)')
            end
            
            if ~opts.refinenoiselvl %% Use globalnoiselvl as fixed noiselvl for ML-estimation
                
                %% Fit the diffusion data with fit_MRI (fit twice for better convergence)
                theta0 = theta;
                displayandlog('Fitting with maximum likelihood estimation (fixed global noiselvl)...');
                theta1 = fit_MRI(ml_fun,data,theta,'fields',fields,'blockSize',opts.fitmriblocksize,'startBlockPos',[0; 0; 0],'mask',mask,'noiseLevel',noiselvl,'optimizer',1,'constraints',constraints,'progressbarUpdateStep',opts.fitmriprogressbar);
                displayandlog('Fitting with maximum likelihood estimation (fixed global noiselvl)...');
                theta2 = fit_MRI(ml_fun,data,theta,'fields',fields,'blockSize',opts.fitmriblocksize,'startBlockPos',floor(-opts.fitmriblocksize/2),'mask',mask,'noiseLevel',noiselvl,'optimizer',1,'constraints',constraints,'progressbarUpdateStep',opts.fitmriprogressbar);
                theta = select_best_fit(ml_fun,data,theta1,theta2,fields,mask,noiselvl,noiselvl,'ml','Parameter vector from fit 2');
                theta = select_best_fit(ml_fun,data,theta,theta0,fields,mask,noiselvl,noiselvl,'ml','Original parameter vector');
                
            else
                
                %% Concatenate parameter vector and noiselvl
                theta = cat(1,theta,noiselvl);                
                lambda = opts.regularizationweight;
                
                if lambda==0 %% Estimate both parameters and noiselvl in ML-estimation (fit twice for better convergence)
                    %% Fit the diffusion data with fit_MRI
                    theta0 = theta;
                    displayandlog('Fitting with maximum likelihood estimation (voxelwise noiselvl estimation)...');
                    theta1 = fit_MRI(ml_fun,data,theta,'fields',fields,'blockSize',opts.fitmriblocksize,'startBlockPos',[0; 0; 0],'mask',mask,'noiseLevel',noiselvl,'optimizer',1,'constraints',constraints,'progressbarUpdateStep',opts.fitmriprogressbar,'numPDFoptpar',1);
                    displayandlog('Fitting with maximum likelihood estimation (voxelwise noiselvl estimation)...');
                    theta2 = fit_MRI(ml_fun,data,theta,'fields',fields,'blockSize',opts.fitmriblocksize,'startBlockPos',floor(-opts.fitmriblocksize/2),'mask',mask,'noiseLevel',noiselvl,'optimizer',1,'constraints',constraints,'progressbarUpdateStep',opts.fitmriprogressbar,'numPDFoptpar',1);
                    theta = select_best_fit(ml_fun,data,theta1,theta2,fields,mask,theta1(npar+1,:,:,:),theta2(npar+1,:,:,:),'ml','Parameter vector from fit 2');
                    theta = select_best_fit(ml_fun,data,theta,theta0,fields,mask,theta(npar+1,:,:,:),theta0(npar+1,:,:,:),'ml','Original parameter vector');
                    
                else %% Estimate both parameters and spatially regularized noiselvl in ML-estimation
                    %% Fit the diffusion data with fit_MRI
                    weights = zeros(npar+1,1);
                    weights(end) = lambda;
                    displayandlog('Fitting with maximum likelihood estimation (spatially regularized noiselvl estimation)...');
                    theta = fit_MRI(ml_fun,data,theta,'fields',fields,'blockSize',opts.fitmriblocksize,'startBlockPos',[0; 0; 0],'blockOverlap',opts.fitmriblockoverlap,'mask',mask,'noiseLevel',noiselvl,'optimizer',3,'progressbarUpdateStep',opts.fitmriprogressbar,'numPDFoptpar',1,'spatialRegularizer','laplacian','voxelSpacing',[1; 1; 1],'spatialRegularizerWeights',weights);
                end 
                
                %% Split noiselvl and parameter vector
                noiselvl = abs(theta(end,:,:,:));
                theta = theta(1:end-1,:,:,:);
                
            end
            
            %% Undo mapping
            if strcmpi(opts.ml_mapping,'erf')
                theta(map_par,:) = erf(theta(map_par,:))/2+0.5;
            elseif strcmpi(opts.ml_mapping,'sin')
                theta(map_par,:) = sin(theta(map_par,:))/2+0.5;
            end

            %% Scale the data and estimated parameters back
            data = data*globalsignallvl;
            noiselvl = noiselvl*globalsignallvl;
            theta(pos_S0,:,:,:) = theta(pos_S0,:,:,:)*globalsignallvl;

            %% Compute some statistics regarding the fit
            [resnorm,logLL] = compute_fit_statistics(fun,data,theta,fields,mask,noiselvl,opts);

            %% Store results in matrices        
            theta_all(:,:,:,slice_idx) = theta;
            mask_all(:,:,slice_idx) = mask;
            noiselvl_all(:,:,:,slice_idx) = noiselvl;
            resnorm_all(:,:,:,slice_idx) = resnorm;
            logLL_all(:,:,:,slice_idx) = logLL;

        else
            displayandlog(['No voxels in mask in slice ' num2str(start_idx) ' to ' num2str(end_idx) '...']);
        end

        %% Set new start index
        start_idx = end_idx+1;

    end

end

%% Write CRLB to disk (requires a lot of memory for large datasets..)
if opts.crlb
    displayandlog('Computing CRLB matrix for all voxels...');
    CRLBmatrix = CramerRaoLowerBound_MRI(theta_all,fun,noiselvl_all,'magnitude',fields_all);
    
    CRLBmatrix(:,~mask_all) = 0;
    CRLBmatrix(CRLBmatrix>1e10) = 1e10;
    CRLBmatrix(CRLBmatrix<-1e10) = -1e10;
    CRLBmatrix(isnan(CRLBmatrix)) = 0;
    
    displayandlog(['Writing CRLB results to disk at location: ' dir_output]);
    crlb2nii(CRLBmatrix,img_header,dir_output,mask);
end

%% Free up some memory
clear data DT_all S0_all
theta_all = theta_all(:,mask_all);
fields_all = fields_all(:,mask_all);
noiselvl_all = noiselvl_all(:,mask_all);
resnorm_all = resnorm_all(:,mask_all);
logLL_all = logLL_all(:,mask_all);

%% Fold back periodical parameters to the interval from 0 to 2*pi
theta_all(period2pi_par,:) = mod(theta_all(period2pi_par,:),2*pi);

%% Write results to disk
displayandlog(['Writing results to disk at location: ' dir_output]);
model_theta2nii(theta_all,fields_all,img_header,dir_output,mask_all,noiselvl_all,resnorm_all,logLL_all);

%% Switch annoying warning on again
if opts.turnofwarning
    warning('on','MATLAB:nearlySingularMatrix');
end

%% Close the log file
fclose('all');

%% Delete uncompressed dwi
if delete_dwi
    delete(fn_dwi);
end

end

function print_opts(opts)
optionfields = fieldnames(opts);
for i=1:numel(optionfields)
    displayandlog([optionfields{i} ': ' num2str(opts.(optionfields{i}))]);
end
end

function boolOut = char2bool(charIn)

if strcmpi(charIn,'true')
    boolOut = true;
elseif strcmpi(charIn,'false')
    boolOut = false;
else
    error(['The string ' charIn ' cannot be converted to a boolean. Use <true> or <false>!'])
end

end

function check_Jacobian(fun,theta,fields,mask)
%% compare the analytical Jacobian with the finite difference Jacobian

theta = theta(:,mask);
theta = theta(:,1:7:min(end,100));
fields = fields(:,mask);
fields = fields(:,1:7:min(end,100));
[S,J] = fun(theta,fields);
Jnumerical = zeros(size(J));

nPar = size(theta,1);
for i=1:nPar
    
    thetadx = theta;
    dx = 1e-3;
    thetadx(i,:) = thetadx(i,:)+dx;
    
    Sdx = fun(thetadx,fields);
    Jnumerical(:,:,i) = (Sdx-S)/dx;
    
end

display('S and Sdx');
[S(:,1), Sdx(:,1)]

display('Theta')
theta(:,1)

for i=1:nPar
    display(['Parameter ' num2str(i) ': ']);
    [J(:,1,i) Jnumerical(:,1,i)]
end

display('S and Sdx');
[S(:,end), Sdx(:,end)]

display('Theta')
theta(:,end)

for i=1:nPar    
    display(['Parameter ' num2str(i) ': ']);
    [J(:,end,i) Jnumerical(:,end,i)]
end

end

function [S,J] = erfmap_fun(theta,fields,fun,map_par)
%% Use erf function to map parameters
theta_map = theta;
theta_map(map_par,:) = erf(theta_map(map_par,:))/2+0.5;

if nargout==1
    S = fun(theta_map,fields);
elseif nargout==2
    [S,J] = fun(theta_map,fields);
    iteration = 1:numel(map_par);
    for i=iteration(map_par)  
        J(:,:,i) = J(:,:,i).*repmat(exp(-theta(i,:).^2)/sqrt(pi),size(S,1),1);
    end
end

end

function [S,J] = sinmap_fun(theta,fields,fun,map_par)
%% Use sin function to map parameters

theta_map = theta;
theta_map(map_par,:) = sin(theta_map(map_par,:))/2+0.5;

if nargout==1
    S = fun(theta_map,fields);
elseif nargout==2
    [S,J] = fun(theta_map,fields);
    iteration = 1:numel(map_par);
    for i=iteration(map_par)        
        J(:,:,i) = J(:,:,i).*repmat(cos(theta(i,:))/2,size(S,1),1);
    end
end

end

function data = load_slices(fn_dwi,slice_idx)
%% Load some slices of the dataset

img = load_untouch_nii(fn_dwi,[],[],[],[],[],slice_idx);
data = img.img;
data = permute(data,[4 1 2 3]);

%% Set zeros and negative values in data to eps to prevent errors in fit_MRI
data(data<=0) = eps;

end

function crlb2nii(CRLBmatrix,dwi_header,dir_output,mask)
%% Write the CRLB matrix to the disk

% Prepare the output location
% Add / or \ to dir_output if necessary
if dir_output(end)~=filesep()
    dir_output = [dir_output filesep()];
end

if ~exist(dir_output,'dir')
    mkdir(dir_output)
end

%% Prepare output file
CRLB_nii.hdr = dwi_header;
CRLB_nii.hdr.dime.dim(1) = 4;
CRLB_nii.hdr.dime.dim(5) = size(CRLBmatrix,1);
CRLB_nii.hdr.dime.datatype = 16;
CRLB_nii.hdr.dime.bitpix = 32;
CRLB_nii.untouch = 1;
CRLBmatrix(:,~mask) = 0;
CRLB_nii.img = permute(CRLBmatrix,[2 3 4 1]);

%% Save the image
fn_niiCRLB = [dir_output 'CRLB_theta.nii.gz'];
save_untouch_nii(CRLB_nii,fn_niiCRLB);

end

function displayandlog(msg)
%% Make a logfile with date and time stamps

global fid
msg = [datestr(now) '   ' msg];
msg = strrep(msg,'\','\\');
fprintf([msg '\n']);
fprintf(fid,[msg '\r\n']);
end

function [resnorm, logLL] = compute_fit_statistics(fun,data,theta,fields,mask,noiselvl,opts)
%% Computes the residual and log likelihood of the fit

dimx = size(theta,2);
dimy = size(theta,3);
dimz = size(theta,4);

data_mask = double(data(:,mask));
data_mask = data_mask(:,:);
theta_mask = theta(:,mask);
theta_mask = theta_mask(:,:);

A_mask = fun(theta_mask,fields(:,mask));
noiselvl_mask = noiselvl(:,:);
noiselvl_mask = noiselvl_mask(:,mask);

if opts.number_priors>0
    data_mask(end-opts.number_priors+1:end,:) = [];
    A_mask(end-opts.number_priors+1:end,:) = [];
end
logLL_mask = logricepdf(data_mask, A_mask, noiselvl_mask, [false, true, false]);
logLL_mask = sum(logLL_mask,1);
residue_mask = data_mask-A_mask;
resnorm_mask = sum(residue_mask.^2,1);

resnorm = zeros(1,dimx*dimy*dimz);
resnorm(mask(:)) = resnorm_mask;
resnorm = reshape(resnorm,dimx,dimy,dimz);

logLL = zeros(1,dimx*dimy,dimz);
logLL(mask(:)) = logLL_mask;
logLL = reshape(logLL,dimx,dimy,dimz);

end

function theta = select_best_fit(fun,data,theta1,theta2,fields,mask,noiselvl1,noiselvl2,type,name)
%% Selects the best parameter vector voxel-wise for a lsq-fit or a ml-fit
size_theta = size(theta1);
theta1 = theta1(:,:);
theta2 = theta2(:,:);
mask = mask(:);
data = double(data(:,:));
data = data(:,mask);
theta1_mask = theta1(:,mask);
theta2_mask = theta2(:,mask);
noiselvl1 = noiselvl1(:,:);
noiselvl1 = noiselvl1(:,mask);
noiselvl2 = noiselvl2(:,:);
noiselvl2 = noiselvl2(:,mask);

if strcmpi(type,'lsq') % select smallest resnorm
    
    A = fun(theta1_mask,fields(:,mask));
    resnorm1 = sum((data-A).^2,1);

    A = fun(theta2_mask,fields(:,mask));
    resnorm2 = sum((data-A).^2,1);

    theta_mask = theta1_mask;
    theta_mask(:,resnorm2<resnorm1) = theta2_mask(:,resnorm2<resnorm1);
    displayandlog([name ' has better fit (lsq) in ' num2str(sum(resnorm2<resnorm1)) ' voxels...']);
    
elseif strcmpi(type,'ml') %select highest log likelihood
    
    A = fun(theta1_mask,fields(:,mask));
    logLL1 = logricepdf(data, A, noiselvl1, [false, true, false]);
    logLL1 = sum(logLL1,1);

    A = fun(theta2_mask,fields(:,mask));
    logLL2 = logricepdf(data, A, noiselvl2, [false, true, false]);
    logLL2 = sum(logLL2,1);

    theta_mask = theta1_mask;
    theta_mask(:,logLL2>logLL1) = theta2_mask(:,logLL2>logLL1);
    displayandlog([name ' has better fit (ml) in ' num2str(sum(logLL2>logLL1)) ' voxels...']);
    
end

theta = theta1;
theta(:,mask) = theta_mask;
theta = reshape(theta,size_theta);

end
