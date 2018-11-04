%% Sample code how to fit various diffusion-tensor models with fitDiffusionModel
% More information about fitDiffusionTensor can be found in the file
% "fitDiffusionModel.m". Code runs fine on Matlab 2013B, compatibility with
% other versions of Matlab is unknown.
% Joor Arkesteijn (joorarkesteijn@gmail.com), Quantitative Imaging Group, 
% Delft University of Technology, 2018.

%% Include fitDiffusionModel into matlab path
addpath(['matlab_code' filesep 'fitDiffusionModel']);
addpath(['matlab_code' filesep 'NIFTI_20130306']); % https://nl.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
addpath(genpath(['matlab_code' filesep 'Tools'])); % http://bigr.nl/people/DirkPoot/

%% Fit a conventional single tensor model
clear opts;
opts.data = ['orig_data' filesep 'data_ecc.nii.gz'];
opts.mask = ['orig_data' filesep 'data_ecc_brain_mask.nii.gz']; % FSL bet can compute a mask
opts.bvals = ['orig_data' filesep 'bval'];
opts.bvecs = ['orig_data' filesep 'bvec'];
opts.model = 'st_original';
opts.out = ['est_data' filesep 'st_original'];
fitDiffusionModel(opts)

%% Fit a ball-and-1-stick model
clear opts;
opts.data = ['orig_data' filesep 'data_ecc.nii.gz'];
opts.mask = ['orig_data' filesep 'data_ecc_brain_mask.nii.gz']; % FSL bet can compute a mask
opts.bvals = ['orig_data' filesep 'bval'];
opts.bvecs = ['orig_data' filesep 'bvec'];
opts.model = 'balland1stick';
opts.out = ['est_data' filesep 'balland1stick'];
fitDiffusionModel(opts)

%% Fit a ball-and-2-sticks model
clear opts;
opts.data = ['orig_data' filesep 'data_ecc.nii.gz'];
opts.mask = ['orig_data' filesep 'data_ecc_brain_mask.nii.gz']; % FSL bet can compute a mask
opts.bvals = ['orig_data' filesep 'bval'];
opts.bvecs = ['orig_data' filesep 'bvec'];
opts.model = 'balland2sticks';
opts.out = ['est_data' filesep 'balland2sticks'];
fitDiffusionModel(opts)

%% Fit a ball-and-2-sticks model with orientation prior
clear opts;
opts.data = ['orig_data' filesep 'data_ecc.nii.gz'];
opts.mask = ['orig_data' filesep 'data_ecc_brain_mask.nii.gz']; % FSL bet can compute a mask
opts.fieldsfile =  ['orig_data' filesep 'D1D2_merged.nii.gz'];
opts.bvals = ['orig_data' filesep 'bval'];
opts.bvecs = ['orig_data' filesep 'bvec'];
opts.model = 'balland2sticks_orientationprior';
opts.orientation_prior = '25';
opts.out = ['est_data' filesep 'balland2sticks_orientationprior'];
fitDiffusionModel(opts)