function [S,J] = st_original_fitmri(theta,Q,opts)
% st_original_fitmri
% theta1: alpha1
% theta2: alpha2
% theta3: alpha3
% theta4: C1 on domain [0,1] (lambda1 = C1*d_max)
% theta5: C2 on domain [0,1] (lambda2 = C2*d_max)
% theta6: C3 on domain [0,1] (lambda3 = C3*d_max)
% theta7: S0

%% Read relevant fields from opts
d_max = opts.d_max;

%% Set some useful constants
nMRI = size(Q,1);
nPar = size(theta,1);
nRange = size(theta,2);

%% Assign some easy names to the parameter vector
theta = real(theta);
alpha1 = theta(1,:);
alpha2 = theta(2,:);
alpha3 = theta(3,:);
lambda1 = theta(4,:)*d_max;
lambda2 = theta(5,:)*d_max;
lambda3 = theta(6,:)*d_max;
S0 = max(theta(7,:),1e-10); %% Prevent infs in Rician likelihood

%% Compute the eigenvectors (http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.html)
cos_alpha1 = cos(alpha1);
sin_alpha1 = sin(alpha1);
cos_alpha2 = cos(alpha2);
sin_alpha2 = sin(alpha2);
cos_alpha3 = cos(alpha3);
sin_alpha3 = sin(alpha3);

EV1 = zeros(3,nRange);
EV1(1,:) = cos_alpha2.*cos_alpha3;
EV1(2,:) = cos_alpha1.*sin_alpha3+sin_alpha1.*sin_alpha2.*cos_alpha3;
EV1(3,:) = sin_alpha1.*sin_alpha3-cos_alpha1.*sin_alpha2.*cos_alpha3;

EV2 = zeros(3,nRange);
EV2(1,:) = -cos_alpha2.*sin_alpha3;
EV2(2,:) = cos_alpha1.*cos_alpha3-sin_alpha1.*sin_alpha2.*sin_alpha3;
EV2(3,:) = sin_alpha1.*cos_alpha3+cos_alpha1.*sin_alpha2.*sin_alpha3;

EV3 = zeros(3,nRange);
EV3(1,:) = sin_alpha2;
EV3(2,:) = -sin_alpha1.*cos_alpha2;
EV3(3,:) = cos_alpha1.*cos_alpha2;

%% Compute the tensor (EV1*EV1'*L1+EV2*EV2'*L2+EV3*EV3'*L3)
D = zeros(6,nRange);

D(1,:) = EV1(1,:).^2.*lambda1 + EV2(1,:).^2.*lambda2 + EV3(1,:).^2.*lambda3;
D(2,:) = EV1(2,:).^2.*lambda1 + EV2(2,:).^2.*lambda2 + EV3(2,:).^2.*lambda3;
D(3,:) = EV1(3,:).^2.*lambda1 + EV2(3,:).^2.*lambda2 + EV3(3,:).^2.*lambda3;
D(4,:) = EV1(1,:).*EV1(2,:).*lambda1 + EV2(1,:).*EV2(2,:).*lambda2 + EV3(1,:).*EV3(2,:).*lambda3;
D(5,:) = EV1(1,:).*EV1(3,:).*lambda1 + EV2(1,:).*EV2(3,:).*lambda2 + EV3(1,:).*EV3(3,:).*lambda3;
D(6,:) = EV1(2,:).*EV1(3,:).*lambda1 + EV2(2,:).*EV2(3,:).*lambda2 + EV3(2,:).*EV3(3,:).*lambda3;

%% Compute the diffusion signal: S
S = repmat(S0,nMRI,1).*exp(Q*D);

%% Compute the Jacobian: J
if nargout==2
    
    %% Pre-allocate some space to fit Jacobian
    J = zeros(nMRI,nRange,nPar);
    
    %% Derivative to alpha1    
    EV1temp = zeros(3,nRange);
    EV1temp(2,:) = -sin_alpha1.*sin_alpha3+cos_alpha1.*sin_alpha2.*cos_alpha3;
    EV1temp(3,:) = cos_alpha1.*sin_alpha3+sin_alpha1.*sin_alpha2.*cos_alpha3;
    
    EV2temp = zeros(3,nRange);
    EV2temp(2,:) = -sin_alpha1.*cos_alpha3-cos_alpha1.*sin_alpha2.*sin_alpha3;
    EV2temp(3,:) = cos_alpha1.*cos_alpha3-sin_alpha1.*sin_alpha2.*sin_alpha3;
    
    EV3temp = zeros(3,nRange);
    EV3temp(2,:) = -cos_alpha1.*cos_alpha2;
    EV3temp(3,:) = -sin_alpha1.*cos_alpha2;

    Dtemp = zeros(6,nRange);
    Dtemp(1,:) = 2*EV1(1,:).*EV1temp(1,:).*lambda1 + 2*EV2(1,:).*EV2temp(1,:).*lambda2 + 2*EV3(1,:).*EV3temp(1,:).*lambda3;
    Dtemp(2,:) = 2*EV1(2,:).*EV1temp(2,:).*lambda1 + 2*EV2(2,:).*EV2temp(2,:).*lambda2 + 2*EV3(2,:).*EV3temp(2,:).*lambda3;
    Dtemp(3,:) = 2*EV1(3,:).*EV1temp(3,:).*lambda1 + 2*EV2(3,:).*EV2temp(3,:).*lambda2 + 2*EV3(3,:).*EV3temp(3,:).*lambda3;
    Dtemp(4,:) = (EV1temp(1,:).*EV1(2,:)+EV1(1,:).*EV1temp(2,:)).*lambda1 + (EV2temp(1,:).*EV2(2,:)+EV2(1,:).*EV2temp(2,:)).*lambda2 + (EV3temp(1,:).*EV3(2,:)+EV3(1,:).*EV3temp(2,:)).*lambda3;
    Dtemp(5,:) = (EV1temp(1,:).*EV1(3,:)+EV1(1,:).*EV1temp(3,:)).*lambda1 + (EV2temp(1,:).*EV2(3,:)+EV2(1,:).*EV2temp(3,:)).*lambda2 + (EV3temp(1,:).*EV3(3,:)+EV3(1,:).*EV3temp(3,:)).*lambda3;
    Dtemp(6,:) = (EV1temp(2,:).*EV1(3,:)+EV1(2,:).*EV1temp(3,:)).*lambda1 + (EV2temp(2,:).*EV2(3,:)+EV2(2,:).*EV2temp(3,:)).*lambda2 + (EV3temp(2,:).*EV3(3,:)+EV3(2,:).*EV3temp(3,:)).*lambda3;
    
    J(:,:,1) = Q*Dtemp.*S;
    
    %% Derivative to alpha2    
    EV1temp = zeros(3,nRange);
    EV1temp(1,:) = -sin_alpha2.*cos_alpha3;
    EV1temp(2,:) = sin_alpha1.*cos_alpha2.*cos_alpha3;
    EV1temp(3,:) = -cos_alpha1.*cos_alpha2.*cos_alpha3;
    
    EV2temp = zeros(3,nRange);
    EV2temp(1,:) = sin_alpha2.*sin_alpha3;
    EV2temp(2,:) = -sin_alpha1.*cos_alpha2.*sin_alpha3;
    EV2temp(3,:) = cos_alpha1.*cos_alpha2.*sin_alpha3;
    
    EV3temp = zeros(3,nRange);
    EV3temp(1,:) = cos_alpha2;
    EV3temp(2,:) = sin_alpha1.*sin_alpha2;
    EV3temp(3,:) = -cos_alpha1.*sin_alpha2;
    
    Dtemp = zeros(6,nRange);
    Dtemp(1,:) = 2*EV1(1,:).*EV1temp(1,:).*lambda1 + 2*EV2(1,:).*EV2temp(1,:).*lambda2 + 2*EV3(1,:).*EV3temp(1,:).*lambda3;
    Dtemp(2,:) = 2*EV1(2,:).*EV1temp(2,:).*lambda1 + 2*EV2(2,:).*EV2temp(2,:).*lambda2 + 2*EV3(2,:).*EV3temp(2,:).*lambda3;
    Dtemp(3,:) = 2*EV1(3,:).*EV1temp(3,:).*lambda1 + 2*EV2(3,:).*EV2temp(3,:).*lambda2 + 2*EV3(3,:).*EV3temp(3,:).*lambda3;
    Dtemp(4,:) = (EV1temp(1,:).*EV1(2,:)+EV1(1,:).*EV1temp(2,:)).*lambda1 + (EV2temp(1,:).*EV2(2,:)+EV2(1,:).*EV2temp(2,:)).*lambda2 + (EV3temp(1,:).*EV3(2,:)+EV3(1,:).*EV3temp(2,:)).*lambda3;
    Dtemp(5,:) = (EV1temp(1,:).*EV1(3,:)+EV1(1,:).*EV1temp(3,:)).*lambda1 + (EV2temp(1,:).*EV2(3,:)+EV2(1,:).*EV2temp(3,:)).*lambda2 + (EV3temp(1,:).*EV3(3,:)+EV3(1,:).*EV3temp(3,:)).*lambda3;
    Dtemp(6,:) = (EV1temp(2,:).*EV1(3,:)+EV1(2,:).*EV1temp(3,:)).*lambda1 + (EV2temp(2,:).*EV2(3,:)+EV2(2,:).*EV2temp(3,:)).*lambda2 + (EV3temp(2,:).*EV3(3,:)+EV3(2,:).*EV3temp(3,:)).*lambda3;
    
    J(:,:,2) = Q*Dtemp.*S;
    
    %% Derivative to alpha3
    EV1temp = zeros(3,nRange);
    EV1temp(1,:) = -cos_alpha2.*sin_alpha3;
    EV1temp(2,:) = cos_alpha1.*cos_alpha3-sin_alpha1.*sin_alpha2.*sin_alpha3;
    EV1temp(3,:) = sin_alpha1.*cos_alpha3+cos_alpha1.*sin_alpha2.*sin_alpha3;
    
    EV2temp = zeros(3,nRange);
    EV2temp(1,:) = -cos_alpha2.*cos_alpha3;
    EV2temp(2,:) = -cos_alpha1.*sin_alpha3-sin_alpha1.*sin_alpha2.*cos_alpha3;
    EV2temp(3,:) = -sin_alpha1.*sin_alpha3+cos_alpha1.*sin_alpha2.*cos_alpha3;
    
    EV3temp = zeros(3,nRange);
    
    Dtemp = zeros(6,nRange);
    Dtemp(1,:) = 2*EV1(1,:).*EV1temp(1,:).*lambda1 + 2*EV2(1,:).*EV2temp(1,:).*lambda2 + 2*EV3(1,:).*EV3temp(1,:).*lambda3;
    Dtemp(2,:) = 2*EV1(2,:).*EV1temp(2,:).*lambda1 + 2*EV2(2,:).*EV2temp(2,:).*lambda2 + 2*EV3(2,:).*EV3temp(2,:).*lambda3;
    Dtemp(3,:) = 2*EV1(3,:).*EV1temp(3,:).*lambda1 + 2*EV2(3,:).*EV2temp(3,:).*lambda2 + 2*EV3(3,:).*EV3temp(3,:).*lambda3;
    Dtemp(4,:) = (EV1temp(1,:).*EV1(2,:)+EV1(1,:).*EV1temp(2,:)).*lambda1 + (EV2temp(1,:).*EV2(2,:)+EV2(1,:).*EV2temp(2,:)).*lambda2 + (EV3temp(1,:).*EV3(2,:)+EV3(1,:).*EV3temp(2,:)).*lambda3;
    Dtemp(5,:) = (EV1temp(1,:).*EV1(3,:)+EV1(1,:).*EV1temp(3,:)).*lambda1 + (EV2temp(1,:).*EV2(3,:)+EV2(1,:).*EV2temp(3,:)).*lambda2 + (EV3temp(1,:).*EV3(3,:)+EV3(1,:).*EV3temp(3,:)).*lambda3;
    Dtemp(6,:) = (EV1temp(2,:).*EV1(3,:)+EV1(2,:).*EV1temp(3,:)).*lambda1 + (EV2temp(2,:).*EV2(3,:)+EV2(2,:).*EV2temp(3,:)).*lambda2 + (EV3temp(2,:).*EV3(3,:)+EV3(2,:).*EV3temp(3,:)).*lambda3;
    
    J(:,:,3) = Q*Dtemp.*S;
    
    %% Derivate to lambda1    
    Dtemp = zeros(6,nRange);
    Dtemp(1,:) = EV1(1,:).^2;
    Dtemp(2,:) = EV1(2,:).^2;
    Dtemp(3,:) = EV1(3,:).^2;
    Dtemp(4,:) = EV1(1,:).*EV1(2,:);
    Dtemp(5,:) = EV1(1,:).*EV1(3,:);
    Dtemp(6,:) = EV1(2,:).*EV1(3,:);

    J(:,:,4) = Q*Dtemp.*S*d_max;
    
    %% Derivative to lambda2
    Dtemp = zeros(6,nRange);
    Dtemp(1,:) = EV2(1,:).^2;
    Dtemp(2,:) = EV2(2,:).^2;
    Dtemp(3,:) = EV2(3,:).^2;
    Dtemp(4,:) = EV2(1,:).*EV2(2,:);
    Dtemp(5,:) = EV2(1,:).*EV2(3,:);
    Dtemp(6,:) = EV2(2,:).*EV2(3,:);
    
    J(:,:,5) = Q*Dtemp.*S*d_max;
    
    %% Derivative to lambda3
    Dtemp = zeros(6,nRange);
    Dtemp(1,:) = EV3(1,:).^2;
    Dtemp(2,:) = EV3(2,:).^2;
    Dtemp(3,:) = EV3(3,:).^2;
    Dtemp(4,:) = EV3(1,:).*EV3(2,:);
    Dtemp(5,:) = EV3(1,:).*EV3(3,:);
    Dtemp(6,:) = EV3(2,:).*EV3(3,:);
    
    J(:,:,6) = Q*Dtemp.*S*d_max;
    
    %% Derivative to S0
    J(:,:,7) = S./repmat(S0,nMRI,1);
    
end

end