function dtv4_theta2nii(theta_mask,dwi_header,dir_output,mask,sigma_mask,resnorm_mask,logLL_mask,opts)
%% dtv4_theta2nii
% This function writes the theta vector to different nifti files for easy
% visual inspection.

dim_theta = size(theta_mask);
thetaConverted = zeros([10 dim_theta(2:end)]);

thetaConverted(1,:) = 0; %set fiso to zero
thetaConverted(2,:) = 0.5; %set f1 to 0.5
thetaConverted(3,:) = theta_mask(1,:);
thetaConverted(4,:) = theta_mask(2,:);
thetaConverted(5,:) = theta_mask(2,:);
thetaConverted(6,:) = theta_mask(3,:);
thetaConverted(7,:) = theta_mask(4,:);
thetaConverted(8,:) = theta_mask(6,:);
thetaConverted(9,:) = theta_mask(6,:);
thetaConverted(10,:) = theta_mask(7,:);

dt_all_theta2nii(thetaConverted,dwi_header,dir_output,mask,sigma_mask,resnorm_mask,logLL_mask,opts);

end