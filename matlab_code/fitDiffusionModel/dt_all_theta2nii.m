function dt_all_theta2nii(theta_mask,dwi_header,dir_output,mask,sigma_mask,resnorm_mask,logLL_mask,opts)
%% dt_all_theta2nii
% theta1: C1 on domain [0,1] (fiso=C1)
% theta2: C2 on domain [0,1] (f1=C2-C1*C2, f2=1-C1-C2+C1*C2)
% theta3: C3 on domain [0,1] (lambda_par = C3*3e-3)
% theta4: C4 on domain [0,1] (lambda_perp1 = C4*lambda_par)
% theta5: C5 on domain [0,1] (lambda_perp2 = C5*lambda_par)
% theta6: theta1
% theta7: phi1
% theta8: theta2
% theta9: phi2
% theta10: S0

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
fn_DTniitheta = [dir_output 'DT_theta.nii.gz'];
save_to_disk(dwi_header,fn_DTniitheta,theta_mask,mask);

fn_DTniimask = [dir_output 'DT_mask.nii.gz'];
save_to_disk(dwi_header,fn_DTniimask,mask(mask)',mask);

fn_DTniisigma = [dir_output 'DT_sigma.nii.gz'];
save_to_disk(dwi_header,fn_DTniisigma,sigma_mask,mask);

fn_DTniiRSS = [dir_output 'DT_RSS.nii.gz'];
save_to_disk(dwi_header,fn_DTniiRSS,resnorm_mask,mask);

fn_DTniilogLL = [dir_output 'DT_logLL.nii.gz'];
save_to_disk(dwi_header,fn_DTniilogLL,logLL_mask,mask);

clear sigma resnorm logLL

%% Read parameters from parameter vector
C1 = theta_mask(1,:);
C2 = theta_mask(2,:);
C3 = theta_mask(3,:);
C4 = theta_mask(4,:);
C5 = theta_mask(5,:);

fiso = C1;
f1 = C2-C1.*C2;
f2 = 1-C1-C2+C1.*C2;
lambda_par = C3*d_max;
lambda_perp1 = C4.*lambda_par;
lambda_perp2 = C5.*lambda_par;

theta1 = theta_mask(6,:);
phi1 = theta_mask(7,:);
theta2 = theta_mask(8,:);
phi2 = theta_mask(9,:);
S0 = max(theta_mask(10,:),1e-10); %% Prevent infs in Rician likelihood

clear C1 C2 C3 C4 C5 theta

%% Assign tensor with largest volume fraction to tensor 1
switchTensor = (f1<f2);
f1temp = f1;
lambda_perp1temp = lambda_perp1;
theta1temp = theta1;
phi1temp = phi1;

f1(switchTensor) = f2(switchTensor);
f2(switchTensor) = f1temp(switchTensor);
lambda_perp1(switchTensor) = lambda_perp2(switchTensor);
lambda_perp2(switchTensor) = lambda_perp1temp(switchTensor);

theta1(switchTensor) = theta2(switchTensor);
theta2(switchTensor) = theta1temp(switchTensor);
phi1(switchTensor) = phi2(switchTensor);
phi2(switchTensor) = phi1temp(switchTensor);

clear f1temp lambda_perp1temp theta1temp phi1temp switchTensor

%% Prepare f1, f2, fiso and S0 and write to disk
fn_DTniif1 = [dir_output 'DT_f1.nii.gz'];
save_to_disk(dwi_header,fn_DTniif1,f1,mask);

fn_DTniif2 = [dir_output 'DT_f2.nii.gz'];
save_to_disk(dwi_header,fn_DTniif2,f2,mask);

fn_DTniifiso = [dir_output 'DT_fiso.nii.gz'];
save_to_disk(dwi_header,fn_DTniifiso,fiso,mask);

fn_DTniiS0 = [dir_output 'DT_S0.nii.gz'];
save_to_disk(dwi_header,fn_DTniiS0,S0,mask);

clear f1 f2 fiso S0

%% Prepare FA1 and MD1 and write to disk
fn_DTniifa1 = [dir_output 'DT_FA1.nii.gz'];
fn_DTniimd1 = [dir_output 'DT_MD1.nii.gz'];

FA1temp1 = sqrt((lambda_par-lambda_perp1).^2+(lambda_perp1-lambda_perp1).^2+(lambda_perp1-lambda_par).^2);
FA1temp2 = sqrt(lambda_par.^2+lambda_perp1.^2+lambda_perp1.^2);
FA1 = sqrt(1/2).*FA1temp1./FA1temp2;
FA1(FA1temp2==0) = 0;
MD1 = (lambda_par+lambda_perp1+lambda_perp1)./3;

save_to_disk(dwi_header,fn_DTniifa1,FA1,mask);
save_to_disk(dwi_header,fn_DTniimd1,MD1,mask);

clear FA1temp1 FA1temp2 FA1 MD1

%% Prepare FA2 and MD2 and write to disk
fn_DTniifa2 = [dir_output 'DT_FA2.nii.gz'];
fn_DTniimd2 = [dir_output 'DT_MD2.nii.gz'];

FA2temp1 = sqrt((lambda_par-lambda_perp2).^2+(lambda_perp2-lambda_perp2).^2+(lambda_perp2-lambda_par).^2);
FA2temp2 = sqrt(lambda_par.^2+lambda_perp2.^2+lambda_perp2.^2);
FA2 = sqrt(1/2).*FA2temp1./FA2temp2;
FA2(FA2temp2==0) = 0;
MD2 = (lambda_par+lambda_perp2+lambda_perp2)./3;

save_to_disk(dwi_header,fn_DTniifa2,FA2,mask);
save_to_disk(dwi_header,fn_DTniimd2,MD2,mask);

clear FA2temp1 FA2temp2 FA2 MD2

%% Prepare EV1 and write to disk
fn_DTniiev1 = [dir_output 'DT_EV1.nii.gz'];

EV1 = zeros(3,nRange);
EV1(1,:) = sin(theta1).*cos(phi1);
EV1(2,:) = sin(theta1).*sin(phi1);
EV1(3,:) = cos(theta1);

save_to_disk(dwi_header,fn_DTniiev1,EV1,mask);

clear theta1 phi1

%% Prepare D1 and write to disk
fn_DTniiDT1 = [dir_output 'DT_DT1.nii.gz'];

D1 = zeros(6,nRange);
D1(1,:) = EV1(1,:).^2.*(lambda_par-lambda_perp1)+lambda_perp1;
D1(2,:) = EV1(2,:).^2.*(lambda_par-lambda_perp1)+lambda_perp1;
D1(3,:) = EV1(3,:).^2.*(lambda_par-lambda_perp1)+lambda_perp1;
D1(4,:) = EV1(1,:).*EV1(2,:).*(lambda_par-lambda_perp1);
D1(5,:) = EV1(1,:).*EV1(3,:).*(lambda_par-lambda_perp1);
D1(6,:) = EV1(2,:).*EV1(3,:).*(lambda_par-lambda_perp1);

save_to_disk(dwi_header,fn_DTniiDT1,D1,mask);

clear EV1 D1

%% Prepare EV2 and write to disk
fn_DTniiev2 = [dir_output 'DT_EV2.nii.gz'];

EV2 = zeros(3,nRange);
EV2(1,:) = sin(theta2).*cos(phi2);
EV2(2,:) = sin(theta2).*sin(phi2);
EV2(3,:) = cos(theta2);

save_to_disk(dwi_header,fn_DTniiev2,EV2,mask);

clear theta2 phi2

%% Prepare D2 and write to disk
fn_DTniiDT2 = [dir_output 'DT_DT2.nii.gz'];

D2 = zeros(6,nRange);
D2(1,:) = EV2(1,:).^2.*(lambda_par-lambda_perp2)+lambda_perp2;
D2(2,:) = EV2(2,:).^2.*(lambda_par-lambda_perp2)+lambda_perp2;
D2(3,:) = EV2(3,:).^2.*(lambda_par-lambda_perp2)+lambda_perp2;
D2(4,:) = EV2(1,:).*EV2(2,:).*(lambda_par-lambda_perp2);
D2(5,:) = EV2(1,:).*EV2(3,:).*(lambda_par-lambda_perp2);
D2(6,:) = EV2(2,:).*EV2(3,:).*(lambda_par-lambda_perp2);

save_to_disk(dwi_header,fn_DTniiDT2,D2,mask);

clear EV2 D2

%% Prepare AD1, AD2, RD1 and RD2 and write to disk
fn_DTniiad1 = [dir_output 'DT_AD1.nii.gz'];
save_to_disk(dwi_header,fn_DTniiad1,lambda_par,mask);

fn_DTniird1 = [dir_output 'DT_RD1.nii.gz'];
save_to_disk(dwi_header,fn_DTniird1,lambda_perp1,mask);

fn_DTniiad2 = [dir_output 'DT_AD2.nii.gz'];
save_to_disk(dwi_header,fn_DTniiad2,lambda_par,mask);

fn_DTniird2 = [dir_output 'DT_RD2.nii.gz'];
save_to_disk(dwi_header,fn_DTniird2,lambda_perp2,mask);

clear all

end

function nifti_file = prepare_nii(dwi_header,dim1,dim5,datatype,bitpix)

% Set nifti header
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

function A = mm3d(B,C,D)

if nargin == 2
    
    A = zeros(size(B));
    A(1,1,:) = B(1,1,:).*C(1,1,:)+B(1,2,:).*C(2,1,:)+B(1,3,:).*C(3,1,:);
    A(1,2,:) = B(1,1,:).*C(1,2,:)+B(1,2,:).*C(2,2,:)+B(1,3,:).*C(3,2,:); 
    A(1,3,:) = B(1,1,:).*C(1,3,:)+B(1,2,:).*C(2,3,:)+B(1,3,:).*C(3,3,:);
    
    A(2,1,:) = B(2,1,:).*C(1,1,:)+B(2,2,:).*C(2,1,:)+B(2,3,:).*C(3,1,:); 
    A(2,2,:) = B(2,1,:).*C(1,2,:)+B(2,2,:).*C(2,2,:)+B(2,3,:).*C(3,2,:); 
    A(2,3,:) = B(2,1,:).*C(1,3,:)+B(2,2,:).*C(2,3,:)+B(2,3,:).*C(3,3,:); 
    
    A(3,1,:) = B(3,1,:).*C(1,1,:)+B(3,2,:).*C(2,1,:)+B(3,3,:).*C(3,1,:); 
    A(3,2,:) = B(3,1,:).*C(1,2,:)+B(3,2,:).*C(2,2,:)+B(3,3,:).*C(3,2,:); 
    A(3,3,:) = B(3,1,:).*C(1,3,:)+B(3,2,:).*C(2,3,:)+B(3,3,:).*C(3,3,:); 
    
elseif nargin == 3
    
    A = mm3d(mm3d(B,C),D);
    
else
    error('Not supported!');
end

end


