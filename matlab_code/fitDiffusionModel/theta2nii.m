function theta2nii(varargin)
% theta2nii recomputes the relevant parameter files from the parametervector

%% Required options
opts.out = ''; %% Opts.out should contain a string with the directory where the output should be written. If the location does not exist it will be created.
opts.model = ''; % Opts.model should contain a string specifying the model to be fit. Various options are possible.
opts.initializationfile = ''; %% Specify a file to use for the initialization, instead of initializing the model from a log-linear fit.

%% Optional options
opts.d_max = '3e-3'; %% Maximum diffusivity of tensor compartment.
opts.d_iso = '3e-3'; %% Diffusivity of isotropic compartment.
opts.MD = '0.8e-3'; %% Mean diffusivity of the tensor.
opts.fieldsfile = ''; %% Specify a file to use for the fields. This could be required for some diffusion-models.

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
opts.d_max = str2double(opts.d_max);
opts.d_iso = str2double(opts.d_iso);
opts.MD = str2double(opts.MD);

%% Prepare output dir
dir_output = opts.out;
if dir_output(end)~=filesep()
    dir_output = [dir_output filesep()];
end
if ~exist(dir_output,'dir')
    mkdir(dir_output)
end

%% Specify the model
if strcmpi(opts.model,'st_original')
    
    %% Define useful functions
    displayandlog('Selecting diffusion model st_original...');
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) st_original_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
elseif strcmpi(opts.model,'stv2')
    
    %% Define useful functions
    displayandlog('Selecting diffusion model stv2...');
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) stv2_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
elseif strcmpi(opts.model,'bitensor')
    
    %% Define useful functions
    displayandlog('Selecting diffusion model bitensor...');
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) bitensor_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
elseif strcmpi(opts.model,'bitensor_mdc')
    
    %% Define useful functions
    displayandlog('Selecting diffusion model bitensor_mdc...');
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) bitensor_mdc_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
elseif strcmpi(opts.model,'bitensor_adc')
    
    %% Define useful functions
    displayandlog('Selecting diffusion model bitensor_adc...');
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) bitensor_adc_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);

elseif strcmpi(opts.model,'stv1')
    
    %% Define useful functions
    displayandlog('Selecting diffusion model stv1...');
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) stv1_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
elseif strcmpi(opts.model,'dtv1')
    
    %% Define useful functions
    displayandlog('Selecting diffusion model dtv1...');
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) dtv1_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
elseif strcmpi(opts.model,'dtv2')
    
    %% Define useful functions
    displayandlog('Selecting diffusion model dtv2...');
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) dtv2_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
elseif strcmpi(opts.model,'dtv3')
    
    %% Define useful functions
    displayandlog('Selecting diffusion model dtv3...');
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) dtv3_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
elseif strcmpi(opts.model,'dtv4')
    
    %% Define useful functions
    displayandlog('Selecting diffusion model dtv4...');
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) dtv4_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
elseif strcmpi(opts.model,'dtv3_constrainedangles')
    
    %% Define useful functions
    displayandlog('Selecting diffusion model dtv3_constrainedangles...');
    model_theta2nii = @(theta_all,fields_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) dtv3_constrainedangles_theta2nii(theta_all,fields_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
elseif strcmpi(opts.model,'balland1stick')
    
    %% Define useful functions
    displayandlog('Selecting diffusion model balland1stick...');
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) balland1stick_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
elseif strcmpi(opts.model,'balland2sticks')
    
    %% Define useful functions
    displayandlog('Selecting diffusion model balland2sticks...');
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) balland2sticks_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);

elseif strcmpi(opts.model,'balland2sticks_orientationprior')
    
    %% Define useful functions
    displayandlog('Selecting diffusion model balland2sticks...');
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) balland2sticks_orientationprior_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);

elseif strcmpi(opts.model,'balland2sticks_2timepoints')
    
    %% Define useful functions
    displayandlog('Selecting diffusion model balland2sticks_2timepoints...');
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) balland2sticks_2timepoints_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
    
elseif strcmpi(opts.model,'balland2sticks_3timepoints')
    
    %% Define useful functions
    dispAmsterdamlayandlog('Selecting diffusion model balland2sticks_3timepoints...');
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) balland2sticks_3timepoints_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);
  
elseif strcmpi(opts.model,'balland2sticks_3timepoints_noprior')
    
    %% Define useful functions
    displayandlog('Selecting diffusion model balland2sticks_3timepoints_noprior...');
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) balland2sticks_3timepoints_noprior_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);

elseif strcmpi(opts.model,'balland2sticks')
    
    %% Define useful functions
    display('Selecting diffusion model balland2sticks...');
    model_theta2nii = @(theta_all,~,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all) balland2sticks_theta2nii(theta_all,img_header,dir_output,mask_all,sigma_all,resnorm_all,logLL_all,opts);

else
    
    error('No valid model has been selected');
    
end

%% Load the header, and mask
img_header = load_untouch_header_only(opts.initializationfile);
img_dim = img_header.dime.dim(2:5);

%% Load all the data
display('Loading initial points from specified file...');
theta_all = load_untouch_nii(opts.initializationfile);
theta_all = permute(theta_all.img,[4 1 2 3]);

%% Load fields
if isempty(opts.fieldsfile)
    fields_all = zeros(1,img_dim(1),img_dim(2),img_dim(3));
else
    fields_all = load_untouch_nii(opts.fieldsfile);
    fields_all = permute(fields_all.img,[4 1 2 3]);
end

%% Load fields
if isempty(opts.initializationfile)
    error('opts.initializationfile cannot be left empty...!');
else
    
end   

%% Make dummy files
mask_all = squeeze(theta_all(end,:,:,:)>0);
noiselvl_all = zeros(1,img_dim(1),img_dim(2),img_dim(3));
resnorm_all = zeros(1,img_dim(1),img_dim(2),img_dim(3));
logLL_all = zeros(1,img_dim(1),img_dim(2),img_dim(3));

%% Free up some memory
theta_all = theta_all(:,mask_all);
fields_all = fields_all(:,mask_all);
noiselvl_all = noiselvl_all(:,mask_all);
resnorm_all = resnorm_all(:,mask_all);
logLL_all = logLL_all(:,mask_all);

%% Write results to disk
display(['Writing results to disk at location: ' dir_output]);
model_theta2nii(theta_all,fields_all,img_header,dir_output,mask_all,noiselvl_all,resnorm_all,logLL_all);

end