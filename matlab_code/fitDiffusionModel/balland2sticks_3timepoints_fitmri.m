function [Smodelprior,Jmodelprior] = balland2sticks_3timepoints_fitmri(theta,fields,Q,opts)
% Computes the predicted Signal (Smodel) and the Jacobian (J) of longitudinal ball-and-two-sticks model, for three timepoints.
% Smodel is a vector containing [S1; S2]
% Q is a vector containing [Q1; Q2].
% theta1: C1A on domain [0,1] (fiso=C1A)
% theta2: C2A on domain [0,1] (f1=C2A-C1A*C2A, f2=1-C1A-C2A+C1A*C2A)
% theta3: C3A on domain [0,1] (d = C3A*d_max)
% theta4: S0A
% theta5: C1B on domain [0,1] (fiso=C1B)
% theta6: C2B on domain [0,1] (f1=C2B-C1B*C2B, f2=1-C1B-C2B+C1B*C2B)
% theta7: C3B on domain [0,1] (d = C3A*d_max)
% theta8: S0B
% theta9: C1B on domain [0,1] (fiso=C1B)
% theta10: C2B on domain [0,1] (f1=C2B-C1B*C2B, f2=1-C1B-C2B+C1B*C2B)
% theta11: C3B on domain [0,1] (d = C3A*d_max)
% theta12: S0B
% theta13: theta1
% theta14: phi1
% theta15: theta2
% theta16: phi2

%% Set useful parameters
nRange = size(theta,2);
nMRI = size(Q,1);
nPar = size(theta,1);
n1 = (1:(nMRI/3))';
n2 = ((n1(end)+1):(nMRI*2/3))';
n3 = ((n2(end)+1):nMRI)';

%% Compute b and q from Q
% b = -Q(:,1)-Q(:,2)-Q(:,3);
% q = sqrt(-Q(:,1:3));
% q(Q(:,4)>0,2) = -q(Q(:,4)>0,2);
% q(Q(:,5)>0,3) = -q(Q(:,5)>0,3);

%% Read relevant fields from opts
%d_max = opts.d_max;
prior_offset = opts.prior_offset;
priorAngle = opts.orientation_prior;

%% Load prior on orientations
EV1_prior = fields(1:3,:);
EV2_prior = fields(4:6,:);

%% Compute the prior
theta1 = theta(13,:);
phi1 = theta(14,:);
theta2 = theta(15,:);
phi2 = theta(16,:);

EV1 = zeros(3,nRange);
EV1(1,:) = sin(theta1).*cos(phi1);
EV1(2,:) = sin(theta1).*sin(phi1);
EV1(3,:) = cos(theta1);

EV2 = zeros(3,nRange);
EV2(1,:) = sin(theta2).*cos(phi2);
EV2(2,:) = sin(theta2).*sin(phi2);
EV2(3,:) = cos(theta2);

acos_arg1 = max(min(1-1e-10,EV1(1,:).*EV1_prior(1,:)+EV1(2,:).*EV1_prior(2,:)+EV1(3,:).*EV1_prior(3,:)),-1+1e-10);
acos_arg2 = max(min(1-1e-10,EV2(1,:).*EV2_prior(1,:)+EV2(2,:).*EV2_prior(2,:)+EV2(3,:).*EV2_prior(3,:)),-1+1e-10);
scaling_factor = priorAngle/0.04; % sigma of prior (ten degrees) divided by sigma of normalized data (0.04), the larger the scaling factor, the smaller the effect of prior.
prior1 = prior_offset(1) + (acos(acos_arg1)*180/pi)/scaling_factor;
prior2 = prior_offset(2) + (acos(acos_arg2)*180/pi)/scaling_factor;

if nargout==1
    
    %% Compute predicted MRI signal
    theta_A = theta([1:3 13:16 4],:);
    theta_B = theta([5:7 13:16 8],:);
    theta_C = theta([9:11 13:16 12],:);
    S_A = balland2sticks_fitmri(theta_A,Q(n1,:),opts);
    S_B = balland2sticks_fitmri(theta_B,Q(n2,:),opts);
    S_C = balland2sticks_fitmri(theta_C,Q(n3,:),opts);

    %% Construct the signal
    Smodelprior = cat(1,S_A,S_B,S_C,prior1,prior2);

    if ~isreal(Smodelprior)
        disp('hallo');
    end
    
elseif nargout>1
    
    %% Compute predicted MRI signal 
    theta_A = theta([1:3 13:16 4],:);
    theta_B = theta([5:7 13:16 8],:);
    theta_C = theta([9:11 13:16 12],:);
    [S_A,J_A] = balland2sticks_fitmri(theta_A,Q(n1,:),opts);
    [S_B,J_B] = balland2sticks_fitmri(theta_B,Q(n2,:),opts);
    [S_C,J_C] = balland2sticks_fitmri(theta_C,Q(n3,:),opts);
    
    %% Construct the complete signal
    Smodelprior = cat(1,S_A,S_B,S_C,prior1,prior2);
    
    %% Compute the Jacobian of the signal
    J = zeros(nMRI,nRange,nPar);
    J(n1,:,[1:4 13:16]) = J_A(:,:,[1:3 8 4:7]);
    J(n2,:,[5:8 13:16]) = J_B(:,:,[1:3 8 4:7]); 
    J(n3,:,9:16) = J_C(:,:,[1:3 8 4:7]);
    
    %% Compute the Jacobian of the prior
    Jprior = zeros(2,nRange,nPar);
    Jprior(1,:,13) = (1/scaling_factor)*(180/pi).*(-1./sqrt(1-acos_arg1.^2)).*(EV1_prior(1,:).*cos(theta1).*cos(phi1)+EV1_prior(2,:).*cos(theta1).*sin(phi1)-EV1_prior(3,:).*sin(theta1));
    Jprior(1,:,14) = (1/scaling_factor)*(180/pi).*(-1./sqrt(1-acos_arg1.^2)).*(-EV1_prior(1,:).*sin(theta1).*sin(phi1)+EV1_prior(2,:).*sin(theta1).*cos(phi1));
    Jprior(2,:,15) = (1/scaling_factor)*(180/pi).*(-1./sqrt(1-acos_arg2.^2)).*(EV2_prior(1,:).*cos(theta2).*cos(phi2)+EV2_prior(2,:).*cos(theta2).*sin(phi2)-EV2_prior(3,:).*sin(theta2));
    Jprior(2,:,16) = (1/scaling_factor)*(180/pi).*(-1./sqrt(1-acos_arg2.^2)).*(-EV2_prior(1,:).*sin(theta2).*sin(phi2)+EV2_prior(2,:).*sin(theta2).*cos(phi2));   
    
    %% Construct the complete Jacobian
    Jmodelprior = cat(1,J,Jprior);
        
end

end
