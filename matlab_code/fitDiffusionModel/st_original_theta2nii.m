function st_original_theta2nii(theta_mask,dwi_header,dir_output,mask,sigma_mask,resnorm_mask,logLL_mask,opts)
%% st_original_theta2nii
% theta1: alpha1
% theta2: alpha2
% theta3: alpha3
% theta4: C1 on domain [0,1] (lambda1 = C1*d_max)
% theta5: C2 on domain [0,1] (lambda2 = C2*d_max)
% theta6: C3 on domain [0,1] (lambda3 = C3*d_max)
% theta7: S0

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
fn_STniitheta = [dir_output 'ST_theta.nii.gz'];
save_to_disk(dwi_header,fn_STniitheta,theta_mask,mask);

fn_STniimask = [dir_output 'ST_mask.nii.gz'];
save_to_disk(dwi_header,fn_STniimask,mask(mask)',mask);

fn_STniisigma = [dir_output 'ST_sigma.nii.gz'];
save_to_disk(dwi_header,fn_STniisigma,sigma_mask,mask);

fn_STniiRSS = [dir_output 'ST_RSS.nii.gz'];
save_to_disk(dwi_header,fn_STniiRSS,resnorm_mask,mask);

fn_STniilogLL = [dir_output 'ST_logLL.nii.gz'];
save_to_disk(dwi_header,fn_STniilogLL,logLL_mask,mask);

clear sigma resnorm logLL

%% Read parameters from parameter vector
alpha1 = theta_mask(1,:);
alpha2 = theta_mask(2,:);
alpha3 = theta_mask(3,:);
lambda1 = theta_mask(4,:)*d_max;
lambda2 = theta_mask(5,:)*d_max;
lambda3 = theta_mask(6,:)*d_max;
S0 = max(theta_mask(7,:),1e-10); %% Prevent infs in Rician likelihood

clear theta

%% Prepare S0 and write to disk
fn_STniiS0 = [dir_output 'ST_S0.nii.gz'];
save_to_disk(dwi_header,fn_STniiS0,S0,mask);

clear S0

%% Prepare FA and MD and write to disk
fn_STniifa = [dir_output 'ST_FA.nii.gz'];
fn_STniimd = [dir_output 'ST_MD.nii.gz'];

FAtemp1 = sqrt((lambda1-lambda2).^2+(lambda2-lambda3).^2+(lambda3-lambda1).^2);
FAtemp2 = sqrt(lambda1.^2+lambda2.^2+lambda3.^2);
FA = sqrt(1/2).*FAtemp1./FAtemp2;
FA(FAtemp2==0) = 0;
MD = (lambda1+lambda2+lambda3)./3;

save_to_disk(dwi_header,fn_STniifa,FA,mask);
save_to_disk(dwi_header,fn_STniimd,MD,mask);

clear FAtemp1 FAtemp2 FA MD

%% Prepare EV1, EV2, EV3, L1, L2 and L3 and write to disk 
fn_STniiv1 = [dir_output 'ST_V1.nii.gz'];
fn_STniiv2 = [dir_output 'ST_V2.nii.gz'];
fn_STniiv3 = [dir_output 'ST_V3.nii.gz'];
fn_STniil1 = [dir_output 'ST_L1.nii.gz'];
fn_STniil2 = [dir_output 'ST_L2.nii.gz'];
fn_STniil3 = [dir_output 'ST_L3.nii.gz'];

cos_alpha1 = cos(alpha1);
sin_alpha1 = sin(alpha1);
cos_alpha2 = cos(alpha2);
sin_alpha2 = sin(alpha2);
cos_alpha3 = cos(alpha3);
sin_alpha3 = sin(alpha3);

EV1 = zeros(3,nRange);
EV1(1,:) = cos_alpha2.*cos_alpha3;
EV1(2,:) = cos_alpha1.*sin_alpha3+sin_alpha1.*sin_alpha2.*cos_alpha3;
EV1(3,:) = sin_alpha1.*sin_alpha3-cos_alpha1.*sin_alpha2.*cos_alpha3;

EV2 = zeros(3,nRange);
EV2(1,:) = -cos_alpha2.*sin_alpha3;
EV2(2,:) = cos_alpha1.*cos_alpha3-sin_alpha1.*sin_alpha2.*sin_alpha3;
EV2(3,:) = sin_alpha1.*cos_alpha3+cos_alpha1.*sin_alpha2.*sin_alpha3;

EV3 = zeros(3,nRange);
EV3(1,:) = sin_alpha2;
EV3(2,:) = -sin_alpha1.*cos_alpha2;
EV3(3,:) = cos_alpha1.*cos_alpha2;

clear alpha1 alpha2 alpha3 cos_alpha1 sin_alpha1 cos_alpha2 sin_alpha2 cos_alpha3 sin_alpha3

%% Sort the eigenvalues and eigenvectors
EW_all = cat(1,lambda1,lambda2,lambda3);
[EW_all,Isort] = sort(EW_all,1,'descend');

lambda1 = EW_all(1,:);
lambda2 = EW_all(2,:);
lambda3 = EW_all(3,:);

EVX = cat(1,EV1(1,:),EV2(1,:),EV3(1,:));
EVY = cat(1,EV1(2,:),EV2(2,:),EV3(2,:));
EVZ = cat(1,EV1(3,:),EV2(3,:),EV3(3,:));

for i=1:nRange
    EVX(:,i) = EVX(Isort(:,i),i);
    EVY(:,i) = EVY(Isort(:,i),i);
    EVZ(:,i) = EVZ(Isort(:,i),i);
end

EV1 = cat(1,EVX(1,:),EVY(1,:),EVZ(1,:));
EV2 = cat(1,EVX(2,:),EVY(2,:),EVZ(2,:));
EV3 = cat(1,EVX(3,:),EVY(3,:),EVZ(3,:));

save_to_disk(dwi_header,fn_STniiv1,EV1,mask);
save_to_disk(dwi_header,fn_STniiv2,EV2,mask);
save_to_disk(dwi_header,fn_STniiv3,EV3,mask);
save_to_disk(dwi_header,fn_STniil1,lambda1,mask);
save_to_disk(dwi_header,fn_STniil2,lambda2,mask);
save_to_disk(dwi_header,fn_STniil3,lambda3,mask);

clear EVX EVY EVZ EW_all Isort

%% Prepare DT and write to disk (EV1*EV1'*L1+EV2*EV2'*L2+EV3*EV3'*L3)
fn_STniiDT = [dir_output 'ST_DT.nii.gz'];

D = zeros(6,nRange);
D(1,:) = EV1(1,:).^2.*lambda1 + EV2(1,:).^2.*lambda2 + EV3(1,:).^2.*lambda3;
D(2,:) = EV1(2,:).^2.*lambda1 + EV2(2,:).^2.*lambda2 + EV3(2,:).^2.*lambda3;
D(3,:) = EV1(3,:).^2.*lambda1 + EV2(3,:).^2.*lambda2 + EV3(3,:).^2.*lambda3;
D(4,:) = EV1(1,:).*EV1(2,:).*lambda1 + EV2(1,:).*EV2(2,:).*lambda2 + EV3(1,:).*EV3(2,:).*lambda3;
D(5,:) = EV1(1,:).*EV1(3,:).*lambda1 + EV2(1,:).*EV2(3,:).*lambda2 + EV3(1,:).*EV3(3,:).*lambda3;
D(6,:) = EV1(2,:).*EV1(3,:).*lambda1 + EV2(2,:).*EV2(3,:).*lambda2 + EV3(2,:).*EV3(3,:).*lambda3;

save_to_disk(dwi_header,fn_STniiDT,D,mask);

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


