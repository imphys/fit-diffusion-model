function dt_all_theta2nii(theta,dwi_header,dir_output,mask,sigma,resnorm,logLL,opts)
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

%% Read relevant fields from opts
d_max = opts.d_max;

%% Find the dimensions in theta, x, y and z
dim_x = size(theta,2);
dim_y = size(theta,3);
dim_z = size(theta,4);

theta = reshape(theta,size(theta,1),dim_x*dim_y*dim_z);

%% Prepare the output location
% Add / or \ to dir_output if necessary
if dir_output(end)~=filesep()
    dir_output = [dir_output filesep()];
end

if ~exist(dir_output,'dir')
    mkdir(dir_output)
end

%% Prepare the output file names
fn_DTniiDT1 = [dir_output 'DT_DT1.nii.gz'];
fn_DTniifa1 = [dir_output 'DT_FA1.nii.gz'];
fn_DTniimd1 = [dir_output 'DT_MD1.nii.gz'];
fn_DTniiv11 = [dir_output 'DT_V11.nii.gz'];
fn_DTniiv21 = [dir_output 'DT_V21.nii.gz'];
fn_DTniiv31 = [dir_output 'DT_V31.nii.gz'];
fn_DTniil11 = [dir_output 'DT_L11.nii.gz'];
fn_DTniil21 = [dir_output 'DT_L21.nii.gz'];
fn_DTniil31 = [dir_output 'DT_L31.nii.gz'];

fn_DTniiDT2 = [dir_output 'DT_DT2.nii.gz'];
fn_DTniifa2 = [dir_output 'DT_FA2.nii.gz'];
fn_DTniimd2 = [dir_output 'DT_MD2.nii.gz'];
fn_DTniiv12 = [dir_output 'DT_V12.nii.gz'];
fn_DTniiv22 = [dir_output 'DT_V22.nii.gz'];
fn_DTniiv32 = [dir_output 'DT_V32.nii.gz'];
fn_DTniil12 = [dir_output 'DT_L12.nii.gz'];
fn_DTniil22 = [dir_output 'DT_L22.nii.gz'];
fn_DTniil32 = [dir_output 'DT_L32.nii.gz'];

fn_DTniif1 = [dir_output 'DT_f1.nii.gz'];
fn_DTniif2 = [dir_output 'DT_f2.nii.gz'];
fn_DTniifiso = [dir_output 'DT_fiso.nii.gz'];
fn_DTniiS0 = [dir_output 'DT_S0.nii.gz'];
fn_DTniitheta = [dir_output 'DT_theta.nii.gz'];

fn_DTniimask = [dir_output 'DT_mask.nii.gz'];
fn_DTniisigma = [dir_output 'DT_sigma.nii.gz'];
fn_DTniiRSS = [dir_output 'DT_RSS.nii.gz'];
fn_DTniilogLL = [dir_output 'DT_logLL.nii.gz'];

%% Prepare the output nifti files
% nifti_file = prepare_nii(dwi_header,dim1,dim5,datatype,bitpix)
DT_DT1 = prepare_nii(dwi_header,4,6,16,32);
DT_fa1 = prepare_nii(dwi_header,3,1,16,32);
DT_md1 = prepare_nii(dwi_header,3,1,16,32);
DT_v11 = prepare_nii(dwi_header,4,3,16,32);
DT_v21 = prepare_nii(dwi_header,4,3,16,32);
DT_v31 = prepare_nii(dwi_header,4,3,16,32);
DT_l11 = prepare_nii(dwi_header,3,1,16,32);
DT_l21 = prepare_nii(dwi_header,3,1,16,32);
DT_l31 = prepare_nii(dwi_header,3,1,16,32);

DT_DT2 = prepare_nii(dwi_header,4,6,16,32);
DT_fa2 = prepare_nii(dwi_header,3,1,16,32);
DT_md2 = prepare_nii(dwi_header,3,1,16,32);
DT_v12 = prepare_nii(dwi_header,4,3,16,32);
DT_v22 = prepare_nii(dwi_header,4,3,16,32);
DT_v32 = prepare_nii(dwi_header,4,3,16,32);
DT_l12 = prepare_nii(dwi_header,3,1,16,32);
DT_l22 = prepare_nii(dwi_header,3,1,16,32);
DT_l32 = prepare_nii(dwi_header,3,1,16,32);

DT_f1 = prepare_nii(dwi_header,3,1,16,32);
DT_f2 = prepare_nii(dwi_header,3,1,16,32);
DT_fiso = prepare_nii(dwi_header,3,1,16,32);
DT_S0 = prepare_nii(dwi_header,3,1,16,32);
DT_theta = prepare_nii(dwi_header,4,size(theta,1),16,32);

DT_mask = prepare_nii(dwi_header,3,1,16,32);
DT_sigma = prepare_nii(dwi_header,3,1,16,32);
DT_RSS = prepare_nii(dwi_header,3,1,16,32);
DT_logLL = prepare_nii(dwi_header,3,1,16,32);

%% Compute all the relevant statistics
nPar = size(theta,1);
nRange = size(theta,2);
theta = real(theta);

%% Define tensors
C1 = theta(1,:);
C2 = theta(2,:);
C3 = theta(3,:);
C4 = theta(4,:);
C5 = theta(5,:);

fiso = C1;
f1 = C2-C1.*C2;
f2 = 1-C1-C2+C1.*C2;
lambda_par = C3*d_max;
lambda_perp1 = C4.*lambda_par;
lambda_perp2 = C5.*lambda_par;

clear C1 C2 C3 C4 C5

theta1 = theta(6,:);
phi1 = theta(7,:);
theta2 = theta(8,:);
phi2 = theta(9,:);
S0 = max(theta(10,:),1e-10); %% Prevent infs in Rician likelihood

clear theta

EW1 = zeros(3,3,nRange);
EW1(1,1,:) = lambda_par;
EW1(2,2,:) = lambda_perp1;
EW1(3,3,:) = lambda_perp1;

EW2 = zeros(3,3,nRange);
EW2(1,1,:) = lambda_par;
EW2(2,2,:) = lambda_perp2;
EW2(3,3,:) = lambda_perp2;

sin_theta1 = sin(theta1);
cos_theta1 = cos(theta1);
sin_phi1 = sin(phi1);
cos_phi1 = cos(phi1);
sin_theta2 = sin(theta2);
cos_theta2 = cos(theta2);
sin_phi2 = sin(phi2);
cos_phi2 = cos(phi2);

R1 = zeros(3,3,nRange);
R1(1,1,:) = sin_theta1.*cos_phi1;
R1(2,1,:) = sin_theta1.*sin_phi1;
R1(3,1,:) = cos_theta1;
R1(1,2,:) = cos_theta1.*cos_phi1;
R1(2,2,:) = cos_theta1.*sin_phi1;
R1(3,2,:) = -sin_theta1;
R1(1,3,:) = -sin_phi1;
R1(2,3,:) = cos_phi1;

R2 = zeros(3,3,nRange);
R2(1,1,:) = sin_theta2.*cos_phi2;
R2(2,1,:) = sin_theta2.*sin_phi2;
R2(3,1,:) = cos_theta2;
R2(1,2,:) = cos_theta2.*cos_phi2;
R2(2,2,:) = cos_theta2.*sin_phi2;
R2(3,2,:) = -sin_theta2;
R2(1,3,:) = -sin_phi2;
R2(2,3,:) = cos_phi2;

D1 = mm3d(R1,EW1,permute(R1,[2 1 3]));
D1 = reshape(D1,9,nRange);
D1 = D1([1 5 9 2 3 6],:);

D2 = mm3d(R2,EW2,permute(R2,[2 1 3]));
D2 = reshape(D2,9,nRange);
D2 = D2([1 5 9 2 3 6],:);

FA1temp1 = sqrt((lambda_par-lambda_perp1).^2+(lambda_perp1-lambda_perp1).^2+(lambda_perp1-lambda_par).^2);
FA1temp2 = sqrt(lambda_par.^2+lambda_perp1.^2+lambda_perp1.^2);
FA1 = sqrt(1/2).*FA1temp1./FA1temp2;
FA1(FA1temp2==0) = 0;
MD1 = (lambda_par+lambda_perp1+lambda_perp1)./3;

FA2temp1 = sqrt((lambda_par-lambda_perp2).^2+(lambda_perp2-lambda_perp2).^2+(lambda_perp2-lambda_par).^2);
FA2temp2 = sqrt(lambda_par.^2+lambda_perp2.^2+lambda_perp2.^2);
FA2 = sqrt(1/2).*FA2temp1./FA2temp2;
FA2(FA2temp2==0) = 0;
MD2 = (lambda_par+lambda_perp2+lambda_perp2)./3;

display('dt_all_theta2nii line183')
whos


%% Assign tensor with largest volume fraction to tensor 1
switchTensor = (f1<f2);
f1temp = f1;
R1temp = R1;
D1temp = D1;
lambda_perp1temp = lambda_perp1;
FA1temp = FA1;
MD1temp = MD1;
f1(switchTensor) = f2(switchTensor);
f2(switchTensor) = f1temp(switchTensor);
R1(:,:,switchTensor) = R2(:,:,switchTensor);
R2(:,:,switchTensor) = R1temp(:,:,switchTensor);
D1(:,switchTensor) = D2(:,switchTensor);
D2(:,switchTensor) = D1temp(:,switchTensor);
lambda_perp1(switchTensor) = lambda_perp2(switchTensor);
lambda_perp2(switchTensor) = lambda_perp1temp(switchTensor);
FA1(switchTensor) = FA2(switchTensor);
FA2(switchTensor) = FA1temp(switchTensor);
MD1(switchTensor) = MD2(switchTensor);
MD2(switchTensor) = MD1temp(switchTensor);

%% Set parameters estimated outside mask to zero
mask = mask(:);
D1(:,~mask) = 0;
D2(:,~mask) = 0;
FA1(~mask) = 0;
FA2(~mask) = 0;
MD1(~mask) = 0;
MD2(~mask) = 0;
R1(:,:,~mask) = 0;
R2(:,:,~mask) = 0;
lambda_par(~mask) = 0;
lambda_perp1(~mask) = 0;
lambda_perp2(~mask) = 0;
f1(~mask) = 0;
f2(~mask) = 0;
fiso(~mask) = 0;
S0(~mask) = 0;
theta(:,~mask) = 0;

%% Reshape, permute, and store
DT_DT1.img = permute(reshape(D1,6,dim_x,dim_y,dim_z),[2 3 4 1]);
DT_fa1.img = permute(reshape(FA1,1,dim_x,dim_y,dim_z),[2 3 4 1]);
DT_md1.img = permute(reshape(MD1,1,dim_x,dim_y,dim_z),[2 3 4 1]);
DT_v11.img = permute(reshape(R1(:,1,:),3,dim_x,dim_y,dim_z),[2 3 4 1]);
DT_v21.img = permute(reshape(R1(:,2,:),3,dim_x,dim_y,dim_z),[2 3 4 1]);
DT_v31.img = permute(reshape(R1(:,3,:),3,dim_x,dim_y,dim_z),[2 3 4 1]);
DT_l11.img = permute(reshape(lambda_par,1,dim_x,dim_y,dim_z),[2 3 4 1]);
DT_l21.img = permute(reshape(lambda_perp1,1,dim_x,dim_y,dim_z),[2 3 4 1]);
DT_l31.img = permute(reshape(lambda_perp1,1,dim_x,dim_y,dim_z),[2 3 4 1]);

DT_DT2.img = permute(reshape(D2,6,dim_x,dim_y,dim_z),[2 3 4 1]);
DT_fa2.img = permute(reshape(FA2,1,dim_x,dim_y,dim_z),[2 3 4 1]);
DT_md2.img = permute(reshape(MD2,1,dim_x,dim_y,dim_z),[2 3 4 1]);
DT_v12.img = permute(reshape(R2(:,1,:),3,dim_x,dim_y,dim_z),[2 3 4 1]);
DT_v22.img = permute(reshape(R2(:,2,:),3,dim_x,dim_y,dim_z),[2 3 4 1]);
DT_v32.img = permute(reshape(R2(:,3,:),3,dim_x,dim_y,dim_z),[2 3 4 1]);
DT_l12.img = permute(reshape(lambda_par,1,dim_x,dim_y,dim_z),[2 3 4 1]);
DT_l22.img = permute(reshape(lambda_perp2,1,dim_x,dim_y,dim_z),[2 3 4 1]);
DT_l32.img = permute(reshape(lambda_perp2,1,dim_x,dim_y,dim_z),[2 3 4 1]);

DT_f1.img = permute(reshape(f1,1,dim_x,dim_y,dim_z),[2 3 4 1]);
DT_f2.img = permute(reshape(f2,1,dim_x,dim_y,dim_z),[2 3 4 1]);
DT_fiso.img = permute(reshape(fiso,1,dim_x,dim_y,dim_z),[2 3 4 1]);
DT_S0.img = permute(reshape(S0,1,dim_x,dim_y,dim_z),[2 3 4 1]);
DT_theta.img = permute(reshape(theta,nPar,dim_x,dim_y,dim_z),[2 3 4 1]);

DT_mask.img = permute(reshape(mask,1,dim_x,dim_y,dim_z),[2 3 4 1]);
DT_sigma.img = sigma;
DT_RSS.img = resnorm;
DT_logLL.img = logLL;

display('dt_all_theta2nii line259')
whos

%% Save dual tensor results
save_untouch_nii(DT_DT1,fn_DTniiDT1);
save_untouch_nii(DT_fa1,fn_DTniifa1);
save_untouch_nii(DT_md1,fn_DTniimd1);
save_untouch_nii(DT_v11,fn_DTniiv11);
save_untouch_nii(DT_v21,fn_DTniiv21);
save_untouch_nii(DT_v31,fn_DTniiv31);
save_untouch_nii(DT_l11,fn_DTniil11);
save_untouch_nii(DT_l21,fn_DTniil21);
save_untouch_nii(DT_l31,fn_DTniil31);

save_untouch_nii(DT_DT2,fn_DTniiDT2);
save_untouch_nii(DT_fa2,fn_DTniifa2);
save_untouch_nii(DT_md2,fn_DTniimd2);
save_untouch_nii(DT_v12,fn_DTniiv12);
save_untouch_nii(DT_v22,fn_DTniiv22);
save_untouch_nii(DT_v32,fn_DTniiv32);
save_untouch_nii(DT_l12,fn_DTniil12);
save_untouch_nii(DT_l22,fn_DTniil22);
save_untouch_nii(DT_l32,fn_DTniil32);

save_untouch_nii(DT_f1,fn_DTniif1);
save_untouch_nii(DT_f2,fn_DTniif2);
save_untouch_nii(DT_fiso,fn_DTniifiso);
save_untouch_nii(DT_S0,fn_DTniiS0);
save_untouch_nii(DT_theta,fn_DTniitheta);

save_untouch_nii(DT_mask,fn_DTniimask);
save_untouch_nii(DT_sigma,fn_DTniisigma);
save_untouch_nii(DT_RSS,fn_DTniiRSS);
save_untouch_nii(DT_logLL,fn_DTniilogLL);

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


