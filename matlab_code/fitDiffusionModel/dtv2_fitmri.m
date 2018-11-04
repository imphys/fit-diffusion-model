function [Smodel,Jconverted] = dtv2_fitmri(theta,Q,opts)
% dtv2_fitmri
% theta1: C1 on domain [0,1] (fiso=C1)
% theta2: C2 on domain [0,1] (f1=(1-C1)*C2, f2=(1-C1)(1-C2)
% theta3: C3 on domain [0,1] (lambda_par = C3*3e-3)
% theta4: C4 on domain [0,1] (lambda_perp = C4*lambda_par)
% theta5: theta1 (polar)
% theta6: phi1 (azimuthal)
% theta7: theta2 (polar)
% theta8: phi2 (azimuthal)
% theta9: S0

%% Read relevant fields from opts
theta = real(theta);
nMRI = size(Q,1);
nPar = size(theta,1);
nRange = size(theta,2);

thetaConverted = zeros(10,nRange);
thetaConverted(1,:) = theta(1,:);
thetaConverted(2,:) = theta(2,:);
thetaConverted(3,:) = theta(3,:);
thetaConverted(4,:) = theta(4,:);
thetaConverted(5,:) = theta(4,:);
thetaConverted(6,:) = theta(5,:);
thetaConverted(7,:) = theta(6,:);
thetaConverted(8,:) = theta(7,:);
thetaConverted(9,:) = theta(8,:);
thetaConverted(10,:) = theta(9,:);

if nargout==1
    
    Smodel = dt_all_fitmri(thetaConverted,Q,opts); 
else
    
    [Smodel,J] = dt_all_fitmri(thetaConverted,Q,opts);
    
    Jconverted = zeros(nMRI,nRange,nPar);
    Jconverted(:,:,1) = J(:,:,1);
    Jconverted(:,:,2) = J(:,:,2);
    Jconverted(:,:,3) = J(:,:,3);
    Jconverted(:,:,4) = J(:,:,4)+J(:,:,5);
    Jconverted(:,:,5) = J(:,:,6);
    Jconverted(:,:,6) = J(:,:,7);
    Jconverted(:,:,7) = J(:,:,8);
    Jconverted(:,:,8) = J(:,:,9);
    Jconverted(:,:,9) = J(:,:,10);
    
end

end
