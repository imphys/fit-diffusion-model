function [Smodel,Jconverted] = stv1_fitmri(theta,Q,opts)
% stv1_fitmri
% theta1: C1 on domain [0,1] (fiso=C1)
% theta2: C3 on domain [0,1] (lambda_par = C3*3e-3)
% theta3: C4 on domain [0,1] (lambda_perp = C4*lambda_par)
% theta4: theta1 (polar)
% theta5: phi1 (azimuthal)
% theta6: S0

%% Read relevant fields from opts
theta = real(theta);
nMRI = size(Q,1);
nPar = size(theta,1);
nRange = size(theta,2);

thetaConverted = zeros(10,nRange);
thetaConverted(1,:) = theta(1,:);
thetaConverted(2,:) = 1; % set f2 to zero
thetaConverted(3,:) = theta(2,:);
thetaConverted(4,:) = theta(3,:);
thetaConverted(5,:) = 0; % lambda_perp2 not defined
thetaConverted(6,:) = theta(4,:);
thetaConverted(7,:) = theta(5,:);
thetaConverted(8,:) = 0; % theta2 not defined
thetaConverted(9,:) = 0; % phi2 not defined
thetaConverted(10,:) = theta(6,:);

if nargout==1
    
    Smodel = dt_all_fitmri(thetaConverted,Q,opts); 
else
    
    [Smodel,J] = dt_all_fitmri(thetaConverted,Q,opts);
    
    Jconverted = zeros(nMRI,nRange,nPar);
    Jconverted(:,:,1) = J(:,:,1);
    Jconverted(:,:,2) = J(:,:,3);
    Jconverted(:,:,3) = J(:,:,4);
    Jconverted(:,:,4) = J(:,:,6);
    Jconverted(:,:,5) = J(:,:,7);
    Jconverted(:,:,6) = J(:,:,10);
    
end

end