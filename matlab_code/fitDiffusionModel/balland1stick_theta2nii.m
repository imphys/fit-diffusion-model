function balland1stick_theta2nii(theta_mask,dwi_header,dir_output,mask,sigma_mask,resnorm_mask,logLL_mask,opts)
%% balland1stick_theta2nii
% theta1: C1 on domain [0,1] (fiso=C1)
% theta2: C2 on domain [0,1] (d = C3*3e-3)
% theta3: theta1
% theta4: phi1
% theta5: S0

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
fn_B2Sniitheta = [dir_output 'B1S_theta.nii.gz'];
save_to_disk(dwi_header,fn_B2Sniitheta,theta_mask,mask);

fn_B2Sniimask = [dir_output 'B1S_mask.nii.gz'];
save_to_disk(dwi_header,fn_B2Sniimask,mask(mask)',mask);

fn_B2Sniisigma = [dir_output 'B1S_sigma.nii.gz'];
save_to_disk(dwi_header,fn_B2Sniisigma,sigma_mask,mask);

fn_B2SniiRSS = [dir_output 'B1S_RSS.nii.gz'];
save_to_disk(dwi_header,fn_B2SniiRSS,resnorm_mask,mask);

fn_B2SniilogLL = [dir_output 'B1S_logLL.nii.gz'];
save_to_disk(dwi_header,fn_B2SniilogLL,logLL_mask,mask);

clear sigma resnorm logLL

%% Read parameters from parameter vector
C1 = theta_mask(1,:);
fiso = C1;
f1 = 1-C1;
C3 = theta_mask(2,:);
d = C3*d_max;
theta1 = theta_mask(3,:);
phi1 = theta_mask(4,:);
S0 = max(theta_mask(5,:),1e-10); %% Prevent infs in Rician likelihood

clear C1 C3 theta

%% Prepare f1, fiso, d and S0 and write to disk
fn_B2Sniif1 = [dir_output 'B1S_f1.nii.gz'];
save_to_disk(dwi_header,fn_B2Sniif1,f1,mask);

fn_B2Sniifiso = [dir_output 'B1S_fiso.nii.gz'];
save_to_disk(dwi_header,fn_B2Sniifiso,fiso,mask);

fn_B2Sniid = [dir_output 'B1S_d.nii.gz'];
save_to_disk(dwi_header,fn_B2Sniid,d,mask);

fn_B2SniiS0 = [dir_output 'B1S_S0.nii.gz'];
save_to_disk(dwi_header,fn_B2SniiS0,S0,mask);

clear f1 fiso S0 d

%% Prepare EV1 and write to disk
fn_B2Sniiev1 = [dir_output 'B1S_EV1.nii.gz'];

EV1 = zeros(3,nRange);
EV1(1,:) = sin(theta1).*cos(phi1);
EV1(2,:) = sin(theta1).*sin(phi1);
EV1(3,:) = cos(theta1);

save_to_disk(dwi_header,fn_B2Sniiev1,EV1,mask);

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

