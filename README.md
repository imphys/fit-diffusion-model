# fit-diffusion-model
Matlab code and atlas data for fitting diffusion models to diffusion-weighted MRI data.

## Introduction
This repository contains Matlab code to fit various diffusion models to diffusion-weighted MRI data. Models include the conventional single tensor model, ball-and-sticks models, and various two-tensor models. Estimation methods include non-linear least squares or maximum-likelihood estimation assuming Rician noise. Furthermore, the repository also contains a model complexity atlas and an orientation atlas computed from 500 elderly subjects [ref], for respectively regularization and model selection purposes. Sample scripts and a sample dataset are provided to demonstrate how to register the atlas data to the subject space, and how to fit various diffusion models to the data.

## Instructions
For registering the atlas data to the subject space:
1. Install [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/) and [DTI-TK](http://dti-tk.sourceforge.net/pmwiki/pmwiki.php)
2. Run the Bash script [Sample_code_coregistration_with_dtitk.sh](Sample_code_coregistration_with_dtitk.sh).

For fitting diffusion models to the diffusion-weighted MRI data:
1. Download ['nifti-tools'](https://nl.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image) from Jimmy Shen and ['estimation tools'](http://bigr.nl/people/DirkPoot/) from Dirk Poot, and put them in your Matlab path. For convenience, these dependencies have also been included in the directory [matlab_code](matlab_code).
2. Run the Matlab script [Sample_code_estimation_with_prior.m](Sample_code_estimation_with_prior.m).

## Questions
For further questions, you may refer to the paper [ref] or contact me at joorarkesteijn(at)gmail.com.
