function dtv2_theta2nii(theta_mask,dwi_header,dir_output,mask,sigma_mask,resnorm_mask,logLL_mask,opts)
%% dtv2_theta2nii
% This function writes the theta vector to different nifti files for easy
% visual inspection.

dim_theta = size(theta_mask);
thetaConverted = zeros([10 dim_theta(2:end)]);

thetaConverted(1,:) = theta_mask(1,:);
thetaConverted(2,:) = theta_mask(2,:);
thetaConverted(3,:) = theta_mask(3,:);
thetaConverted(4,:) = theta_mask(4,:);
thetaConverted(5,:) = theta_mask(4,:);
thetaConverted(6,:) = theta_mask(5,:);
thetaConverted(7,:) = theta_mask(6,:);
thetaConverted(8,:) = theta_mask(7,:);
thetaConverted(9,:) = theta_mask(8,:);
thetaConverted(10,:) = theta_mask(9,:);

dt_all_theta2nii(thetaConverted,dwi_header,dir_output,mask,sigma_mask,resnorm_mask,logLL_mask,opts);

end