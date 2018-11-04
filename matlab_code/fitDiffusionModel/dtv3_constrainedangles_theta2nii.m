function dtv3_constrainedangles_theta2nii(theta_mask,fields_mask,dwi_header,dir_output,mask,sigma_mask,resnorm_mask,logLL_mask,opts)
%% dtv3_constrainedangles_theta2nii
% This function writes the theta vector to different nifti files for easy
% visual inspection.

%% Save some memory
dim_theta = size(theta_mask);
thetaConverted = zeros([10 dim_theta(2:end)]);

thetaConverted(1,:) = 0; %set fiso to zero
thetaConverted(2,:) = theta_mask(1,:);
thetaConverted(3,:) = theta_mask(2,:);
thetaConverted(4,:) = theta_mask(3,:);
thetaConverted(5,:) = theta_mask(3,:);
thetaConverted(6,:) = fields_mask(1,:);
thetaConverted(7,:) = fields_mask(2,:);
thetaConverted(8,:) = fields_mask(3,:);
thetaConverted(9,:) = fields_mask(4,:);
thetaConverted(10,:) = theta_mask(4,:);

dt_all_theta2nii(thetaConverted,dwi_header,dir_output,mask,sigma_mask,resnorm_mask,logLL_mask,opts);

end