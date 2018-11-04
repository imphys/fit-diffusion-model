function stv1_theta2nii(theta_mask,dwi_header,dir_output,mask,sigma_mask,resnorm_mask,logLL_mask,opts)
%% stv1_theta2nii
% This function writes the theta vector to different nifti files for easy
% visual inspection.

dim_theta = size(theta_mask);

thetaConverted = zeros([10 dim_theta(2:end)]);
thetaConverted(1,:) = theta_mask(1,:);
thetaConverted(2,:) = 1; % set f2 to zero
thetaConverted(3,:) = theta_mask(2,:);
thetaConverted(4,:) = theta_mask(3,:);
thetaConverted(5,:) = 0; % lambda_perp2 not defined
thetaConverted(6,:) = theta_mask(4,:);
thetaConverted(7,:) = theta_mask(5,:);
thetaConverted(8,:) = 0; % theta2 not defined
thetaConverted(9,:) = 0; % phi2 not defined
thetaConverted(10,:) = theta_mask(6,:);

dt_all_theta2nii(thetaConverted,dwi_header,dir_output,mask,sigma_mask,resnorm_mask,logLL_mask,opts);

end