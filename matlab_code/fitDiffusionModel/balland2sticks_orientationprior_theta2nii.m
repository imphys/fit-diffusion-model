function balland2sticks_orientationprior_theta2nii(theta_mask,dwi_header,dir_output,mask,sigma_mask,resnorm_mask,logLL_mask,opts)
%% balland2sticks_orientationprior_theta2nii
% theta1: C1 on domain [0,1] (fiso=C1)
% theta2: C2 on domain [0,1] (f1=C2-C1*C2, f2=1-C1-C2+C1*C2)
% theta3: C3 on domain [0,1] (d = C3*3e-3)
% theta4: theta1
% theta5: phi1
% theta6: theta2
% theta7: phi2
% theta8: S0

%% Define some useful parameters
nRange = size(theta_mask,2);
theta_mask = real(theta_mask);

%% Prepare the output location
% Add / or \ to dir_output if necessary
if dir_output(end)~=filesep()
    dir_output = [dir_output filesep()];
end

if ~exist(dir_output,'dir')
    mkdir(dir_output)
end

%% Read relevant fields from opts
d_max = opts.d_max;

%% Prepare theta, mask, sigma, RSS and logLL and write to disk
fn_B2Sniitheta = [dir_output 'B2S_theta.nii.gz'];
save_to_disk(dwi_header,fn_B2Sniitheta,theta_mask,mask);

fn_B2Sniimask = [dir_output 'B2S_mask.nii.gz'];
save_to_disk(dwi_header,fn_B2Sniimask,mask(mask)',mask);

fn_B2Sniisigma = [dir_output 'B2S_sigma.nii.gz'];
save_to_disk(dwi_header,fn_B2Sniisigma,sigma_mask,mask);

fn_B2SniiRSS = [dir_output 'B2S_RSS.nii.gz'];
save_to_disk(dwi_header,fn_B2SniiRSS,resnorm_mask,mask);

fn_B2SniilogLL = [dir_output 'B2S_logLL.nii.gz'];
save_to_disk(dwi_header,fn_B2SniilogLL,logLL_mask,mask);

clear sigma resnorm logLL

%% Read parameters from parameter vector
C1 = theta_mask(1,:);
C2 = theta_mask(2,:);
fiso = C1;
f1 = C2-C1.*C2;
f2 = 1-C1-C2+C1.*C2;
C3 = theta_mask(3,:);
d = C3*d_max;
theta1 = theta_mask(4,:);
phi1 = theta_mask(5,:);
theta2 = theta_mask(6,:);
phi2 = theta_mask(7,:);
S0 = max(theta_mask(8,:),1e-10); %% Prevent infs in Rician likelihood

clear C1 C2 C3 theta

%% Assign tensor with largest volume fraction to tensor 1
% switchTensor = (f1<f2);
% f1temp = f1;
% theta1temp = theta1;
% phi1temp = phi1;
% 
% f1(switchTensor) = f2(switchTensor);
% f2(switchTensor) = f1temp(switchTensor);
% theta1(switchTensor) = theta2(switchTensor);
% theta2(switchTensor) = theta1temp(switchTensor);
% phi1(switchTensor) = phi2(switchTensor);
% phi2(switchTensor) = phi1temp(switchTensor);
% 
% clear f1temp lambda_perp1temp theta1temp phi1temp switchTensor

%% Prepare f1, f2, fiso, d and S0 and write to disk
fn_B2Sniif1 = [dir_output 'B2S_f1.nii.gz'];
save_to_disk(dwi_header,fn_B2Sniif1,f1,mask);

fn_B2Sniif2 = [dir_output 'B2S_f2.nii.gz'];
save_to_disk(dwi_header,fn_B2Sniif2,f2,mask);

fn_B2Sniifiso = [dir_output 'B2S_fiso.nii.gz'];
save_to_disk(dwi_header,fn_B2Sniifiso,fiso,mask);

fn_B2Sniid = [dir_output 'B2S_d.nii.gz'];
save_to_disk(dwi_header,fn_B2Sniid,d,mask);

fn_B2SniiS0 = [dir_output 'B2S_S0.nii.gz'];
save_to_disk(dwi_header,fn_B2SniiS0,S0,mask);

clear f1 f2 fiso S0 d

%% Prepare EV1 and write to disk
fn_B2Sniiev1 = [dir_output 'B2S_EV1.nii.gz'];

EV1 = zeros(3,nRange);
EV1(1,:) = sin(theta1).*cos(phi1);
EV1(2,:) = sin(theta1).*sin(phi1);
EV1(3,:) = cos(theta1);

save_to_disk(dwi_header,fn_B2Sniiev1,EV1,mask);

clear theta1 phi1 EV1

%% Prepare EV2 and write to disk
fn_B2Sniiev2 = [dir_output 'B2S_EV2.nii.gz'];

EV2 = zeros(3,nRange);
EV2(1,:) = sin(theta2).*cos(phi2);
EV2(2,:) = sin(theta2).*sin(phi2);
EV2(3,:) = cos(theta2);

save_to_disk(dwi_header,fn_B2Sniiev2,EV2,mask);

clear all

end

function nifti_file = prepare_nii(dwi_header,dim1,dim5,datatype,bitpix)

nifti_file.hdr = dwi_header;
nifti_file.hdr.dime.scl_slope = 1;
nifti_file.hdr.dime.scl_inter = 0;
nifti_file.hdr.dime.dim(1) = dim1;
nifti_file.hdr.dime.dim(5) = dim5;
nifti_file.hdr.dime.datatype = datatype;
nifti_file.hdr.dime.bitpix = bitpix;
nifti_file.untouch = 1;
nifti_file.img = zeros(nifti_file.hdr.dime.dim(2),nifti_file.hdr.dime.dim(3),nifti_file.hdr.dime.dim(4),nifti_file.hdr.dime.dim(5));

end

function save_to_disk(dwi_header,data_fn,data_mask,mask)

% Prepare nifti
dim5 = size(data_mask,1);
dim1 = 3+(dim5>1);
nii = prepare_nii(dwi_header,dim1,dim5,16,32);

% Prepare data
data = zeros([dim5 size(mask)]);
data(:,mask) = data_mask;
nii.img = permute(data,[2 3 4 1]);

% Save data
save_untouch_nii(nii,data_fn);

end

