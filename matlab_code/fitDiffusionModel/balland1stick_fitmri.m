function [Smodel,J] = balland1stick_fitmri(theta,Q,opts)
% balland1stick_fitmri computes the predicted Signal (Smodel) and the Jacobian (J) of the ball-and-1-stick model.
% theta1: C1 on domain [0,1] (fiso=C1)
% theta2: C2 on domain [0,1] (d = C2*d_max)
% theta3: theta1
% theta4: phi1
% theta5: S0

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

%% Compute predicted MRI signal
C1 = theta(1,:);
fiso = C1;
f1 = 1-C1;
C2 = theta(2,:);
d = C2*d_max;
theta1 = theta(3,:);
phi1 = theta(4,:);
S0 = max(theta(5,:),1e-10); %% Prevent infs in Rician likelihood

EV1 = zeros(3,nRange);
EV1(1,:) = sin(theta1).*cos(phi1);
EV1(2,:) = sin(theta1).*sin(phi1);
EV1(3,:) = cos(theta1);

S1 = repmat(S0,nMRI,1).*exp(-repmat(d,nMRI,1).*(q*EV1).^2);
Siso = repmat(S0,nMRI,1).*exp(-b*d);
f1S1 = repmat(f1,nMRI,1).*S1;
fisoSiso = repmat(fiso,nMRI,1).*Siso;
Smodel = (f1S1+fisoSiso);

if nargout>1
    
    % theta1: C1 on domain [0,1] (fiso=C1)
    % theta2: C2 on domain [0,1] (d = C2*d_max)
    % theta3: theta1
    % theta4: phi1
    % theta5: S0
    
    J = zeros(nMRI,nRange,nPar);
    
    J(:,:,1) = (Siso-S1);
    
    dS_dd_tensor1 = -(q*EV1).^2.*f1S1;
    dS_dd_diso = repmat(-b,1,nRange).*fisoSiso;
    
    J(:,:,2) = (dS_dd_tensor1+dS_dd_diso)*d_max;
    
    EVtemp = zeros(3,nRange);
    EVtemp(1,:) = cos(theta1).*cos(phi1);
    EVtemp(2,:) = cos(theta1).*sin(phi1);
    EVtemp(3,:) = -sin(theta1);
    
    J(:,:,3) = -repmat(2*d,nMRI,1).*(q*EV1).*(q*EVtemp).*f1S1;
    
    EVtemp = zeros(3,nRange);
    EVtemp(1,:) = sin(theta1).*-sin(phi1);
    EVtemp(2,:) = sin(theta1).*cos(phi1);
    
    J(:,:,4) = -repmat(2*d,nMRI,1).*(q*EV1).*(q*EVtemp).*f1S1;
    
    J(:,:,5) = Smodel./repmat(S0,nMRI,1);

end

end
