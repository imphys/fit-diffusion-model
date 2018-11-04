function [Smodel,Jconverted] = dtv3_constrainedangles_fitmri(theta,fields,Q,opts)
% dtv3_constrainedangles_fitmri
% theta1: C2 on domain [0,1] (f1=C2, f2=1-C2)
% theta2: C3 on domain [0,1] (lambda_par = C3*3e-3)
% theta3: C4 on domain [0,1] (lambda_perp = C4*lambda_par)
% fields1: theta1 (polar)
% fields2: phi1 (azimuthal)
% fields3: theta2 (polar)
% fields4: phi2 (azimuthal)
% theta4: S0

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
thetaConverted(6,:) = fields(1,:);
thetaConverted(7,:) = fields(2,:);
thetaConverted(8,:) = fields(3,:);
thetaConverted(9,:) = fields(4,:);
thetaConverted(10,:) = theta(4,:);

if nargout==1
    
    Smodel = dt_all_fitmri(thetaConverted,Q,opts); 
else
    
    [Smodel,J] = dt_all_fitmri(thetaConverted,Q,opts);
    
    Jconverted = zeros(nMRI,nRange,nPar);
    Jconverted(:,:,1) = J(:,:,2);
    Jconverted(:,:,2) = J(:,:,3);
    Jconverted(:,:,3) = J(:,:,4)+J(:,:,5);
    Jconverted(:,:,4) = J(:,:,10);
    
end

end