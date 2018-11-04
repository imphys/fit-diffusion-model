function [Smodel,Jconverted] = dtv3_fitmri(theta,Q,opts)
% dtv3_fitmri
% theta1: C2 on domain [0,1] (f1=C2, f2=1-C2)
% theta2: C3 on domain [0,1] (lambda_par = C3*3e-3)
% theta3: C4 on domain [0,1] (lambda_perp = C4*lambda_par)
% theta4: theta1 (polar)
% theta5: phi1 (azimuthal)
% theta6: theta2 (polar)
% theta7: phi2 (azimuthal)
% theta8: S0

%% Read relevant fields from opts
theta = real(theta);
nMRI = size(Q,1);
nPar = size(theta,1);
nRange = size(theta,2);

thetaConverted = zeros(10,nRange);
thetaConverted(1,:) = 0; %set fiso to zero
thetaConverted(2,:) = theta(1,:);
thetaConverted(3,:) = theta(2,:);
thetaConverted(4,:) = theta(3,:);
thetaConverted(5,:) = theta(3,:);
thetaConverted(6,:) = theta(4,:);
thetaConverted(7,:) = theta(5,:);
thetaConverted(8,:) = theta(6,:);
thetaConverted(9,:) = theta(7,:);
thetaConverted(10,:) = theta(8,:);

if nargout==1
    
    Smodel = dt_all_fitmri(thetaConverted,Q,opts); 
else
    
    [Smodel,J] = dt_all_fitmri(thetaConverted,Q,opts);
    
    Jconverted = zeros(nMRI,nRange,nPar);
    Jconverted(:,:,1) = J(:,:,2);
    Jconverted(:,:,2) = J(:,:,3);
    Jconverted(:,:,3) = J(:,:,4)+J(:,:,5);
    Jconverted(:,:,4) = J(:,:,6);
    Jconverted(:,:,5) = J(:,:,7);
    Jconverted(:,:,6) = J(:,:,8);
    Jconverted(:,:,7) = J(:,:,9);
    Jconverted(:,:,8) = J(:,:,10);
    
end

end