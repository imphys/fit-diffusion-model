function [Smodel,Jmodel] = balland2sticks_3timepoints_noprior_fitmri(theta,~,Q,opts)
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

if nargout==1
    
    %% Compute predicted MRI signal
    theta_A = theta([1:3 13:16 4],:);
    theta_B = theta([5:7 13:16 8],:);
    theta_C = theta([9:11 13:16 12],:);
    S_A = balland2sticks_fitmri(theta_A,Q(n1,:),opts);
    S_B = balland2sticks_fitmri(theta_B,Q(n2,:),opts);
    S_C = balland2sticks_fitmri(theta_C,Q(n3,:),opts);

    %% Construct the signal
    Smodel = cat(1,S_A,S_B,S_C);
    
elseif nargout>1
    
    %% Compute predicted MRI signal 
    theta_A = theta([1:3 13:16 4],:);
    theta_B = theta([5:7 13:16 8],:);
    theta_C = theta([9:11 13:16 12],:);
    [S_A,J_A] = balland2sticks_fitmri(theta_A,Q(n1,:),opts);
    [S_B,J_B] = balland2sticks_fitmri(theta_B,Q(n2,:),opts);
    [S_C,J_C] = balland2sticks_fitmri(theta_C,Q(n3,:),opts);
    
    %% Construct the complete signal
    Smodel = cat(1,S_A,S_B,S_C);
    
    %% Compute the Jacobian of the signal
    J = zeros(nMRI,nRange,nPar);
    J(n1,:,[1:4 13:16]) = J_A(:,:,[1:3 8 4:7]);
    J(n2,:,[5:8 13:16]) = J_B(:,:,[1:3 8 4:7]); 
    J(n3,:,9:16) = J_C(:,:,[1:3 8 4:7]);
       
    %% Construct the complete Jacobian
    Jmodel = J;
        
end

end
