function [Smodelprior,Jmodelprior] = balland2sticks_orientationprior_fitmri(theta,fields,Q,opts)
% balland2sticks_orientationprior_fitmri computes the predicted Signal (Smodel) and the Jacobian (J) of the ball-and-2-sticks model.
% theta1: C1 on domain [0,1] (fiso=C1)
% theta2: C2 on domain [0,1] (f1=C2-C1*C2, f2=1-C1-C2+C1*C2)
% theta3: C3 on domain [0,1] (d = C3*d_max)
% theta4: theta1
% theta5: phi1
% theta6: theta2
% theta7: phi2
% theta8: S0

%% Set useful parameters
nRange = size(theta,2);
nMRI = size(Q,1);
nPar = size(theta,1);

%% Compute b and q from Q
b = -Q(:,1)-Q(:,2)-Q(:,3);
q = sqrt(-Q(:,1:3));
q(Q(:,4)>0,2) = -q(Q(:,4)>0,2);
q(Q(:,5)>0,3) = -q(Q(:,5)>0,3);

%% Read relevant fields from opts
d_max = opts.d_max;
prior_offset = opts.prior_offset;
priorAngle = opts.orientation_prior;

%% Load prior on orientations
EV1_prior = fields(1:3,:);
EV2_prior = fields(4:6,:);

%% Compute predicted MRI signal
C1 = theta(1,:);
C2 = theta(2,:);
fiso = C1;
f1 = C2-C1.*C2;
f2 = 1-C1-C2+C1.*C2;
C3 = theta(3,:);
d = C3*d_max;
theta1 = theta(4,:);
phi1 = theta(5,:);
theta2 = theta(6,:);
phi2 = theta(7,:);
S0 = max(theta(8,:),1e-10); %% Prevent infs in Rician likelihood

EV1 = zeros(3,nRange);
EV1(1,:) = sin(theta1).*cos(phi1);
EV1(2,:) = sin(theta1).*sin(phi1);
EV1(3,:) = cos(theta1);

EV2 = zeros(3,nRange);
EV2(1,:) = sin(theta2).*cos(phi2);
EV2(2,:) = sin(theta2).*sin(phi2);
EV2(3,:) = cos(theta2);

S1 = repmat(S0,nMRI,1).*exp(-repmat(d,nMRI,1).*(q*EV1).^2);
S2 = repmat(S0,nMRI,1).*exp(-repmat(d,nMRI,1).*(q*EV2).^2);
Siso = repmat(S0,nMRI,1).*exp(-b*d);
f1S1 = repmat(f1,nMRI,1).*S1;
f2S2 = repmat(f2,nMRI,1).*S2;
fisoSiso = repmat(fiso,nMRI,1).*Siso;
Smodel = (f1S1+f2S2+fisoSiso);

acos_arg1 = max(min(1-1e-10,EV1(1,:).*EV1_prior(1,:)+EV1(2,:).*EV1_prior(2,:)+EV1(3,:).*EV1_prior(3,:)),-1+1e-10);
acos_arg2 = max(min(1-1e-10,EV2(1,:).*EV2_prior(1,:)+EV2(2,:).*EV2_prior(2,:)+EV2(3,:).*EV2_prior(3,:)),-1+1e-10);
scaling_factor = priorAngle/0.04; % sigma of prior divided by sigma of normalized data (0.04, assuming SNR of 25), the larger the scaling factor, the smaller the effect of prior.
prior1 = prior_offset(1) + (acos(acos_arg1)*180/pi)/scaling_factor;
prior2 = prior_offset(2) + (acos(acos_arg2)*180/pi)/scaling_factor;
Smodelprior = cat(1,Smodel,prior1,prior2);

if nargout>1
    
    % theta1: C1 on domain [0,1] (fiso=C1)
    % theta2: C2 on domain [0,1] (f1=C2-C1*C2, f2=1-C1-C2+C1*C2)
    % theta3: C3 on domain [0,1] (d=C3*3e-3)
    % theta4: theta1
    % theta5: phi1
    % theta6: theta2
    % theta7: phi2
    % theta8: S0
    
    J = zeros(nMRI,nRange,nPar);
    
    J(:,:,1) = (Siso-repmat(C2,nMRI,1).*S1+repmat(C2-1,nMRI,1).*S2);
    J(:,:,2) = (repmat(1-C1,nMRI,1).*S1+repmat(C1-1,nMRI,1).*S2);
    
    dS_dd_tensor1 = -(q*EV1).^2.*f1S1;
    dS_dd_tensor2 = -(q*EV2).^2.*f2S2;
    dS_dd_diso = repmat(-b,1,nRange).*fisoSiso;
    
    J(:,:,3) = (dS_dd_tensor1+dS_dd_tensor2+dS_dd_diso)*d_max;
    
    EVtemp = zeros(3,nRange);
    EVtemp(1,:) = cos(theta1).*cos(phi1);
    EVtemp(2,:) = cos(theta1).*sin(phi1);
    EVtemp(3,:) = -sin(theta1);
    
    J(:,:,4) = -repmat(2*d,nMRI,1).*(q*EV1).*(q*EVtemp).*f1S1;
    
    EVtemp = zeros(3,nRange);
    EVtemp(1,:) = sin(theta1).*-sin(phi1);
    EVtemp(2,:) = sin(theta1).*cos(phi1);
    
    J(:,:,5) = -repmat(2*d,nMRI,1).*(q*EV1).*(q*EVtemp).*f1S1;

    EVtemp = zeros(3,nRange);
    EVtemp(1,:) = cos(theta2).*cos(phi2);
    EVtemp(2,:) = cos(theta2).*sin(phi2);
    EVtemp(3,:) = -sin(theta2);
    
    J(:,:,6) = -repmat(2*d,nMRI,1).*(q*EV2).*(q*EVtemp).*f2S2;

    EVtemp = zeros(3,nRange);
    EVtemp(1,:) = sin(theta2).*-sin(phi2);
    EVtemp(2,:) = sin(theta2).*cos(phi2);
    
    J(:,:,7) = -repmat(2*d,nMRI,1).*(q*EV2).*(q*EVtemp).*f2S2;
    
    J(:,:,8) = Smodel./repmat(S0,nMRI,1);
    
    Jprior = zeros(2,nRange,nPar);
    Jprior(1,:,4) = (1/scaling_factor)*(180/pi).*(-1./sqrt(1-acos_arg1.^2)).*(EV1_prior(1,:).*cos(theta1).*cos(phi1)+EV1_prior(2,:).*cos(theta1).*sin(phi1)-EV1_prior(3,:).*sin(theta1));
    Jprior(1,:,5) = (1/scaling_factor)*(180/pi).*(-1./sqrt(1-acos_arg1.^2)).*(-EV1_prior(1,:).*sin(theta1).*sin(phi1)+EV1_prior(2,:).*sin(theta1).*cos(phi1));
    Jprior(2,:,6) = (1/scaling_factor)*(180/pi).*(-1./sqrt(1-acos_arg2.^2)).*(EV2_prior(1,:).*cos(theta2).*cos(phi2)+EV2_prior(2,:).*cos(theta2).*sin(phi2)-EV2_prior(3,:).*sin(theta2));
    Jprior(2,:,7) = (1/scaling_factor)*(180/pi).*(-1./sqrt(1-acos_arg2.^2)).*(-EV2_prior(1,:).*sin(theta2).*sin(phi2)+EV2_prior(2,:).*sin(theta2).*cos(phi2));   
    Jmodelprior = cat(1,J,Jprior);

end

end
