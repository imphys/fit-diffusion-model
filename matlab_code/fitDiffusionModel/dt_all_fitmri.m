function [Smodel,J] = dt_all_fitmri(theta,Q,opts)
% dt_all_fitmri computes the predicted Signal (Smodel) and the Jacobian (J) of loosely constrained dual-tensor model.
% theta1: C1 on domain [0,1] (fiso=C1)
% theta2: C2 on domain [0,1] (f1=C2-C1*C2, f2=1-C1-C2+C1*C2)
% theta3: C3 on domain [0,1] (lambda_par = C3*3e-3)
% theta4: C4 on domain [0,1] (lambda_perp1 = C4*lambda_par)
% theta5: C5 on domain [0,1] (lambda_perp2 = C5*lambda_par)
% theta6: theta1
% theta7: phi1
% theta8: theta2
% theta9: phi2
% theta10: S0

%% Read relevant fields from opts
d_max = opts.d_max;
d_iso = opts.d_iso;

%% Define some useful constants
nMRI = size(Q,1);
nPar = size(theta,1);
nRange = size(theta,2);
theta = real(theta);

%% Compute predicted MRI signal
C1 = theta(1,:);
C2 = theta(2,:);
C3 = theta(3,:);
C4 = theta(4,:);
C5 = theta(5,:);

fiso = C1;
f1 = C2-C1.*C2;
f2 = 1-C1-C2+C1.*C2;
lambda_par = C3*d_max;
lambda_perp1 = C4.*lambda_par;
lambda_perp2 = C5.*lambda_par;

theta1 = theta(6,:);
phi1 = theta(7,:);
theta2 = theta(8,:);
phi2 = theta(9,:);
S0 = max(theta(10,:),1e-10); %% Prevent infs in Rician likelihood

EV1 = zeros(3,nRange);
EV1(1,:) = sin(theta1).*cos(phi1);
EV1(2,:) = sin(theta1).*sin(phi1);
EV1(3,:) = cos(theta1);

D1 = zeros(6,nRange);
D1(1,:) = EV1(1,:).^2.*(lambda_par-lambda_perp1)+lambda_perp1;
D1(2,:) = EV1(2,:).^2.*(lambda_par-lambda_perp1)+lambda_perp1;
D1(3,:) = EV1(3,:).^2.*(lambda_par-lambda_perp1)+lambda_perp1;
D1(4,:) = EV1(1,:).*EV1(2,:).*(lambda_par-lambda_perp1);
D1(5,:) = EV1(1,:).*EV1(3,:).*(lambda_par-lambda_perp1);
D1(6,:) = EV1(2,:).*EV1(3,:).*(lambda_par-lambda_perp1);

EV2 = zeros(3,nRange);
EV2(1,:) = sin(theta2).*cos(phi2);
EV2(2,:) = sin(theta2).*sin(phi2);
EV2(3,:) = cos(theta2);

D2 = zeros(6,nRange);
D2(1,:) = EV2(1,:).^2.*(lambda_par-lambda_perp2)+lambda_perp2;
D2(2,:) = EV2(2,:).^2.*(lambda_par-lambda_perp2)+lambda_perp2;
D2(3,:) = EV2(3,:).^2.*(lambda_par-lambda_perp2)+lambda_perp2;
D2(4,:) = EV2(1,:).*EV2(2,:).*(lambda_par-lambda_perp2);
D2(5,:) = EV2(1,:).*EV2(3,:).*(lambda_par-lambda_perp2);
D2(6,:) = EV2(2,:).*EV2(3,:).*(lambda_par-lambda_perp2);

Diso = zeros(6,nRange);
Diso([1 2 3],:) = d_iso;

S1 = repmat(S0,nMRI,1).*exp(Q*D1);
S2 = repmat(S0,nMRI,1).*exp(Q*D2);
f1S1 = repmat(f1,nMRI,1).*S1;
f2S2 = repmat(f2,nMRI,1).*S2;
Siso = repmat(S0,nMRI,1).*exp(Q*Diso);
Smodel = f1S1+f2S2+repmat(fiso,nMRI,1).*Siso;

if nargout>1

    % theta1: C1 on domain [0,1] (fiso=C1)
    % theta2: C2 on domain [0,1] (f1=C2-C1*C2, f2=1-C1-C2+C1*C2)
    % theta3: C3 on domain [0,1] (lambda_par = C3*3e-3)
    % theta4: C4 on domain [0,1] (lambda_perp1 = C4*lambda_par)
    % theta5: C5 on domain [0,1] (lambda_perp2 = C5*lambda_par)
    % theta6: theta1
    % theta7: phi1
    % theta8: theta2
    % theta9: phi2
    % theta10: S0
 
    J = zeros(nMRI,nRange,nPar);
    
    J(:,:,1) = (Siso-repmat(C2,nMRI,1).*S1+repmat(C2-1,nMRI,1).*S2);
    J(:,:,2) = (repmat(1-C1,nMRI,1).*S1+repmat(C1-1,nMRI,1).*S2);
    
    Dtemp = zeros(6,nRange);
    Dtemp(1,:) = EV1(1,:).^2;
    Dtemp(2,:) = EV1(2,:).^2;
    Dtemp(3,:) = EV1(3,:).^2;
    Dtemp(4,:) = EV1(1,:).*EV1(2,:);
    Dtemp(5,:) = EV1(1,:).*EV1(3,:);
    Dtemp(6,:) = EV1(2,:).*EV1(3,:);
    dS_dlambda_par1 = Q*Dtemp.*f1S1;
    
    Dtemp = zeros(6,nRange);
    Dtemp(1,:) = EV2(1,:).^2;
    Dtemp(2,:) = EV2(2,:).^2;
    Dtemp(3,:) = EV2(3,:).^2;
    Dtemp(4,:) = EV2(1,:).*EV2(2,:);
    Dtemp(5,:) = EV2(1,:).*EV2(3,:);
    Dtemp(6,:) = EV2(2,:).*EV2(3,:);
    dS_dlambda_par2 = Q*Dtemp.*f2S2;
    
    Dtemp = zeros(6,nRange);
    Dtemp(1,:) = 1-EV1(1,:).^2;
    Dtemp(2,:) = 1-EV1(2,:).^2;
    Dtemp(3,:) = 1-EV1(3,:).^2;
    Dtemp(4,:) = -EV1(1,:).*EV1(2,:);
    Dtemp(5,:) = -EV1(1,:).*EV1(3,:);
    Dtemp(6,:) = -EV1(2,:).*EV1(3,:);
    dS_dlambda_perp1 = Q*Dtemp.*f1S1;
    
    Dtemp = zeros(6,nRange);
    Dtemp(1,:) = 1-EV2(1,:).^2;
    Dtemp(2,:) = 1-EV2(2,:).^2;
    Dtemp(3,:) = 1-EV2(3,:).^2;
    Dtemp(4,:) = -EV2(1,:).*EV2(2,:);
    Dtemp(5,:) = -EV2(1,:).*EV2(3,:);
    Dtemp(6,:) = -EV2(2,:).*EV2(3,:);
    dS_dlambda_perp2 = Q*Dtemp.*f2S2;
    
    J(:,:,3) = (dS_dlambda_par1+dS_dlambda_par2+dS_dlambda_perp1.*repmat(C4,nMRI,1)+dS_dlambda_perp2.*repmat(C5,nMRI,1))*d_max;
    J(:,:,4) = dS_dlambda_perp1.*repmat(lambda_par,nMRI,1);
    J(:,:,5) = dS_dlambda_perp2.*repmat(lambda_par,nMRI,1);
    
    EVtemp = zeros(3,nRange);
    EVtemp(1,:) = cos(theta1).*cos(phi1);
    EVtemp(2,:) = cos(theta1).*sin(phi1);
    EVtemp(3,:) = -sin(theta1);
    
    Dtemp = zeros(6,nRange);
    Dtemp(1,:) = 2*EVtemp(1,:).*EV1(1,:).*(lambda_par-lambda_perp1);
    Dtemp(2,:) = 2*EVtemp(2,:).*EV1(2,:).*(lambda_par-lambda_perp1);
    Dtemp(3,:) = 2*EVtemp(3,:).*EV1(3,:).*(lambda_par-lambda_perp1);
    Dtemp(4,:) = (EVtemp(1,:).*EV1(2,:)+EVtemp(2,:).*EV1(1,:)).*(lambda_par-lambda_perp1);
    Dtemp(5,:) = (EVtemp(1,:).*EV1(3,:)+EVtemp(3,:).*EV1(1,:)).*(lambda_par-lambda_perp1);
    Dtemp(6,:) = (EVtemp(2,:).*EV1(3,:)+EVtemp(3,:).*EV1(2,:)).*(lambda_par-lambda_perp1);
    
    J(:,:,6) = Q*Dtemp.*f1S1;
    
    EVtemp = zeros(3,nRange);
    EVtemp(1,:) = sin(theta1).*-sin(phi1);
    EVtemp(2,:) = sin(theta1).*cos(phi1);
    
    Dtemp = zeros(6,nRange);
    Dtemp(1,:) = 2*EVtemp(1,:).*EV1(1,:).*(lambda_par-lambda_perp1);
    Dtemp(2,:) = 2*EVtemp(2,:).*EV1(2,:).*(lambda_par-lambda_perp1);
    Dtemp(3,:) = 2*EVtemp(3,:).*EV1(3,:).*(lambda_par-lambda_perp1);
    Dtemp(4,:) = (EVtemp(1,:).*EV1(2,:)+EVtemp(2,:).*EV1(1,:)).*(lambda_par-lambda_perp1);
    Dtemp(5,:) = (EVtemp(1,:).*EV1(3,:)+EVtemp(3,:).*EV1(1,:)).*(lambda_par-lambda_perp1);
    Dtemp(6,:) = (EVtemp(2,:).*EV1(3,:)+EVtemp(3,:).*EV1(2,:)).*(lambda_par-lambda_perp1);
    
    J(:,:,7) = Q*Dtemp.*f1S1;
    
    EVtemp = zeros(3,nRange);
    EVtemp(1,:) = cos(theta2).*cos(phi2);
    EVtemp(2,:) = cos(theta2).*sin(phi2);
    EVtemp(3,:) = -sin(theta2);
    
    Dtemp = zeros(6,nRange);
    Dtemp(1,:) = 2*EVtemp(1,:).*EV2(1,:).*(lambda_par-lambda_perp2);
    Dtemp(2,:) = 2*EVtemp(2,:).*EV2(2,:).*(lambda_par-lambda_perp2);
    Dtemp(3,:) = 2*EVtemp(3,:).*EV2(3,:).*(lambda_par-lambda_perp2);
    Dtemp(4,:) = (EVtemp(1,:).*EV2(2,:)+EVtemp(2,:).*EV2(1,:)).*(lambda_par-lambda_perp2);
    Dtemp(5,:) = (EVtemp(1,:).*EV2(3,:)+EVtemp(3,:).*EV2(1,:)).*(lambda_par-lambda_perp2);
    Dtemp(6,:) = (EVtemp(2,:).*EV2(3,:)+EVtemp(3,:).*EV2(2,:)).*(lambda_par-lambda_perp2);
    
    J(:,:,8) = Q*Dtemp.*f2S2;
    
    EVtemp = zeros(3,nRange);
    EVtemp(1,:) = sin(theta2).*-sin(phi2);
    EVtemp(2,:) = sin(theta2).*cos(phi2);
    
    Dtemp = zeros(6,nRange);
    Dtemp(1,:) = 2*EVtemp(1,:).*EV2(1,:).*(lambda_par-lambda_perp2);
    Dtemp(2,:) = 2*EVtemp(2,:).*EV2(2,:).*(lambda_par-lambda_perp2);
    Dtemp(3,:) = 2*EVtemp(3,:).*EV2(3,:).*(lambda_par-lambda_perp2);
    Dtemp(4,:) = (EVtemp(1,:).*EV2(2,:)+EVtemp(2,:).*EV2(1,:)).*(lambda_par-lambda_perp2);
    Dtemp(5,:) = (EVtemp(1,:).*EV2(3,:)+EVtemp(3,:).*EV2(1,:)).*(lambda_par-lambda_perp2);
    Dtemp(6,:) = (EVtemp(2,:).*EV2(3,:)+EVtemp(3,:).*EV2(2,:)).*(lambda_par-lambda_perp2);

    J(:,:,9) = Q*Dtemp.*f2S2;
    
    J(:,:,10) = Smodel./repmat(S0,nMRI,1);

end

end