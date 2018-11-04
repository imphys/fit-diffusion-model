function balland2sticks_2timepoints_theta2nii(theta_mask,dwi_header,dir_output,mask,sigma_mask,resnorm_mask,logLL_mask,opts)
% theta1: C1A on domain [0,1] (fiso=C1A)
% theta2: C2A on domain [0,1] (f1=C2A-C1A*C2A, f2=1-C1A-C2A+C1A*C2A)
% theta3: C3A on domain [0,1] (d = C3A*d_max)
% theta4: S0A
% theta5: C1B on domain [0,1] (fiso=C1B)
% theta6: C2B on domain [0,1] (f1=C2B-C1B*C2B, f2=1-C1B-C2B+C1B*C2B)
% theta7: C3B on domain [0,1] (d = C3A*d_max)
% theta8: S0B
% theta9: theta1
% theta10: phi1
% theta11: theta2
% theta12: phi2

%% Prepare the output location
% Add / or \ to dir_output if necessary
if dir_output(end)~=filesep()
    dir_output_A = [dir_output filesep() 'A' filesep()];
    dir_output_B = [dir_output filesep() 'B' filesep()];
else
    dir_output_A = [dir_output 'A' filesep()];
    dir_output_B = [dir_output 'B' filesep()];
end

%% Take care that theta2nii function does not apply sorting (otherwise correspondence will be lost)!!!
theta_mask_A = theta_mask([1:3 9:12 4],:);
balland2sticks_orientationprior_theta2nii(theta_mask_A,dwi_header,dir_output_A,mask,sigma_mask,resnorm_mask,logLL_mask,opts);
clear theta_mask_A

%% Take care that theta2nii function does not apply sorting (otherwise correspondence will be lost)!!!
theta_mask_B = theta_mask([5:7 9:12 8],:);
balland2sticks_orientationprior_theta2nii(theta_mask_B,dwi_header,dir_output_B,mask,sigma_mask,resnorm_mask,logLL_mask,opts);
clear theta_mask_B

end