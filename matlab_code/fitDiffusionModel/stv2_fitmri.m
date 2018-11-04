function [Smodel,Jconverted] = stv2_fitmri(theta,Q,opts)
% stv2_fitmri
% theta1: C3 on domain [0,1] (lambda_par = C3*3e-3)
% theta2: C4 on domain [0,1] (lambda_perp = C4*lambda_par)
% theta3: theta1 (polar)
% theta4: phi1 (azimuthal)
% theta5: S0

%% Read relevant fields from opts
theta = real(theta);
nMRI = size(Q,1);
nPar = size(theta,1);
nRange = size(theta,2);

thetaConverted = zeros(10,nRange);
thetaConverted(1,:) = 0; % set fiso to zero
thetaConverted(2,:) = 1; % set f1 to one
thetaConverted(3,:) = theta(1,:);
thetaConverted(4,:) = theta(2,:);
thetaConverted(5,:) = 0; % value of lambda_perp1 not used
thetaConverted(6,:) = theta(3,:);
thetaConverted(7,:) = theta(4,:);
thetaConverted(8,:) = 0; % value of theta2 not used
thetaConverted(9,:) = 0; % value of phi2 not used
thetaConverted(10,:) = theta(5,:);

if nargout==1
    
    Smodel = dt_all_fitmri(thetaConverted,Q,opts); 
else
    
    [Smodel,J] = dt_all_fitmri(thetaConverted,Q,opts);
    
    Jconverted = zeros(nMRI,nRange,nPar);
    Jconverted(:,:,1) = J(:,:,3);
    Jconverted(:,:,2) = J(:,:,4)+J(:,:,5);
    Jconverted(:,:,3) = J(:,:,6);
    Jconverted(:,:,4) = J(:,:,7);
    Jconverted(:,:,5) = J(:,:,10);
    
end

end