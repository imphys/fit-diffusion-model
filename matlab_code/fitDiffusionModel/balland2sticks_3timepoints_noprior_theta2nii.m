function balland2sticks_3timepoints_noprior_theta2nii(theta_mask,dwi_header,dir_output,mask,sigma_mask,resnorm_mask,logLL_mask,opts)
% theta1: C1A on domain [0,1] (fiso=C1A)
% theta2: C2A on domain [0,1] (f1=C2A-C1A*C2A, f2=1-C1A-C2A+C1A*C2A)
% theta3: C3A on domain [0,1] (d = C3A*d_max)
% theta4: S0A
% theta5: C1B on domain [0,1] (fiso=C1B)
% theta6: C2B on domain [0,1] (f1=C2B-C1B*C2B, f2=1-C1B-C2B+C1B*C2B)
% theta7: C3B on domain [0,1] (d = C3A*d_max)
% theta8: S0B
% theta9: C1B on domain [0,1] (fiso=C1B)
% theta10: C2B on domain [0,1] (f1=C2B-C1B*C2B, f2=1-C1B-C2B+C1B*C2B)
% theta11: C3B on domain [0,1] (d = C3A*d_max)
% theta12: S0B
% theta13: theta1
% theta14: phi1
% theta15: theta2
% theta16: phi2

%% Prepare the output location
% Add / or \ to dir_output if necessary
if dir_output(end)~=filesep()
    dir_output_A = [dir_output filesep() 'A' filesep()];
    dir_output_B = [dir_output filesep() 'B' filesep()];
    dir_output_C = [dir_output filesep() 'C' filesep()];
else
    dir_output_A = [dir_output 'A' filesep()];
    dir_output_B = [dir_output 'B' filesep()];
    dir_output_C = [dir_output 'C' filesep()];
end

%% Assign tensor with largest volume fraction to tensor 1
C1_tmp = theta_mask(1,:);
C2_tmp = theta_mask(2,:);
f1_sum = C2_tmp-C1_tmp.*C2_tmp;
f2_sum = 1-C1_tmp-C2_tmp+C1_tmp.*C2_tmp;

C1_tmp = theta_mask(5,:);
C2_tmp = theta_mask(6,:);
f1_sum = f1_sum+C2_tmp-C1_tmp.*C2_tmp;
f2_sum = f2_sum+1-C1_tmp-C2_tmp+C1_tmp.*C2_tmp;

C1_tmp = theta_mask(9,:);
C2_tmp = theta_mask(10,:);
f1_sum = f1_sum+C2_tmp-C1_tmp.*C2_tmp;
f2_sum = f2_sum+1-C1_tmp-C2_tmp+C1_tmp.*C2_tmp;

switchTensor = (f1_sum<f2_sum);

% Switch stick fractions
theta_mask(2,switchTensor) = 1-theta_mask(2,switchTensor);
theta_mask(6,switchTensor) = 1-theta_mask(6,switchTensor);
theta_mask(10,switchTensor) = 1-theta_mask(10,switchTensor);

% Switch orientations
theta1temp = theta_mask(13,:);
phi1temp = theta_mask(14,:);
theta_mask(13,switchTensor) = theta_mask(15,switchTensor);
theta_mask(14,switchTensor) = theta_mask(16,switchTensor);
theta_mask(15,switchTensor) = theta1temp(switchTensor);
theta_mask(16,switchTensor) = phi1temp(switchTensor);

clear C1_tmp C2_tmp f1_sum f2_sum theta1temp phi1temp switchTensor

%% Take care that theta2nii function does not apply sorting (otherwise correspondence will be lost)!!!
theta_mask_A = theta_mask([1:3 13:16 4],:);
balland2sticks_orientationprior_theta2nii(theta_mask_A,dwi_header,dir_output_A,mask,sigma_mask,resnorm_mask,logLL_mask,opts);
clear theta_mask_A

%% Take care that theta2nii function does not apply sorting (otherwise correspondence will be lost)!!!
theta_mask_B = theta_mask([5:7 13:16 8],:);
balland2sticks_orientationprior_theta2nii(theta_mask_B,dwi_header,dir_output_B,mask,sigma_mask,resnorm_mask,logLL_mask,opts);
clear theta_mask_B

%% Take care that theta2nii function does not apply sorting (otherwise correspondence will be lost)!!!
theta_mask_C = theta_mask([9:11 13:16 12],:);
balland2sticks_orientationprior_theta2nii(theta_mask_C,dwi_header,dir_output_C,mask,sigma_mask,resnorm_mask,logLL_mask,opts);
clear theta_mask_C

end