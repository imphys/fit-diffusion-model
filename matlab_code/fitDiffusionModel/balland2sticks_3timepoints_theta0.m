function theta0_cell = balland2sticks_3timepoints_theta0(~,~,data,Q,fields,mask,opts)
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

%% Read relevant fields from opts
d_max = opts.d_max;

%% Find indices to seperate timepoints
nMRI = size(Q,1);
n1 = (1:(nMRI/3))';
n2 = ((n1(end)+1):(nMRI*2/3))';
n3 = ((n2(end)+1):nMRI)';

%% Pre-allocate space for parameter matrix
nPar = 16;
theta0a = zeros([nPar size(mask)]);
theta0b = zeros([nPar size(mask)]);

%% Only compute theta0 inside mask to save time
theta0a_mask = theta0a(:,mask(:));
theta0b_mask = theta0b(:,mask(:));
data_mask = data(:,mask(:));
EV1_mask = fields(1:3,mask(:));
EV2_mask = fields(4:6,mask(:));

%% Re-compute S0 and DT seperately for each timepoint
[S0_mask_A,DT_mask_A]=fitSingleTensor(data_mask(n1,:),Q(n1,:));
[S0_mask_B,DT_mask_B]=fitSingleTensor(data_mask(n2,:),Q(n2,:));
[S0_mask_C,DT_mask_C]=fitSingleTensor(data_mask(n3,:),Q(n3,:));

for i=1:size(theta0a_mask,2)
    
    [~,~,~,EW_A] = computeStatisticsFromTensor(DT_mask_A(:,i));
    [~,~,~,EW_B] = computeStatisticsFromTensor(DT_mask_B(:,i));
    [~,~,~,EW_C] = computeStatisticsFromTensor(DT_mask_C(:,i));
    
    EV1 = EV1_mask(1:3,i);
    EV2 = EV2_mask(1:3,i);
    %EV1 = [rand; rand; rand]; EV1 = EV1./norm(EV1);
    %EV2 = [rand; rand; rand]; EV2 = EV2./norm(EV2);
    
    % Put EW in physically sensible range
    EW_A = max(EW_A,0.1e-4);
    EW_A = min(EW_A,d_max);
    EW_B = max(EW_B,0.1e-4);
    EW_B = min(EW_B,d_max);
    EW_C = max(EW_C,0.1e-4);
    EW_C = min(EW_C,d_max);
    
    % Approximate new eigenwaarden
    lambda_par_A = EW_A(1)+EW_A(2)-EW_A(3);
    lambda_par_B = EW_B(1)+EW_B(2)-EW_B(3);
    lambda_par_C = EW_C(1)+EW_C(2)-EW_C(3);
      
    % Approximate C3
    C3A = lambda_par_A/d_max;
    C3A = min(max(C3A,1e-2),0.99);
    lambda_par_A = C3A*d_max;
    
    C3B = lambda_par_B/d_max;
    C3B = min(max(C3B,1e-2),0.99);
    lambda_par_B = C3B*d_max;
    
    C3C = lambda_par_C/d_max;
    C3C = min(max(C3C,1e-2),0.99);
    lambda_par_C = C3C*d_max;
    
    [phi1,theta1] = cart2sph(EV1(1),EV1(2),EV1(3));
    theta1 = pi/2-theta1;
    [phi2,theta2] = cart2sph(EV2(1),EV2(2),EV2(3));
    theta2 = pi/2-theta2;
    
    % Use non-negative least squares to compute the volume fractions
    tensor1 = EV1*EV1'*lambda_par_A;
    tensor1 = tensor1([1 5 9 4 7 8])';
    tensor2 = EV2*EV2'*lambda_par_A;
    tensor2 = tensor2([1 5 9 4 7 8])';
    ball = [lambda_par_A;lambda_par_A;lambda_par_A;0;0;0];
    tensor1_signal = S0_mask_A(i)*exp(Q(n1,:)*tensor1);
    tensor2_signal = S0_mask_A(i)*exp(Q(n1,:)*tensor2);
    ball_signal = S0_mask_A(i)*exp(Q(n1,:)*ball);
    A = [tensor1_signal,tensor2_signal,ball_signal];
    x = lsqnonneg(A,data_mask(n1,i));
    x = max(x,0.01);
    x = x./sum(x);
    C1A = x(3);
    C2A = x(1)/(1-x(3)); 
    
    % Use non-negative least squares to compute the volume fractions
    tensor1 = EV1*EV1'*lambda_par_B;
    tensor1 = tensor1([1 5 9 4 7 8])';
    tensor2 = EV2*EV2'*lambda_par_B;
    tensor2 = tensor2([1 5 9 4 7 8])';
    ball = [lambda_par_B;lambda_par_B;lambda_par_B;0;0;0];
    tensor1_signal = S0_mask_B(i)*exp(Q(n2,:)*tensor1);
    tensor2_signal = S0_mask_B(i)*exp(Q(n2,:)*tensor2);
    ball_signal = S0_mask_B(i)*exp(Q(n2,:)*ball);
    A = [tensor1_signal,tensor2_signal,ball_signal];
    x = lsqnonneg(A,data_mask(n2,i));
    x = max(x,0.01);
    x = x./sum(x);
    C1B = x(3);
    C2B = x(1)/(1-x(3));  
    
    % Use non-negative least squares to compute the volume fractions
    tensor1 = EV1*EV1'*lambda_par_C;
    tensor1 = tensor1([1 5 9 4 7 8])';
    tensor2 = EV2*EV2'*lambda_par_C;
    tensor2 = tensor2([1 5 9 4 7 8])';
    ball = [lambda_par_C;lambda_par_C;lambda_par_C;0;0;0];
    tensor1_signal = S0_mask_C(i)*exp(Q(n3,:)*tensor1);
    tensor2_signal = S0_mask_B(i)*exp(Q(n3,:)*tensor2);
    ball_signal = S0_mask_C(i)*exp(Q(n3,:)*ball);
    A = [tensor1_signal,tensor2_signal,ball_signal];
    x = lsqnonneg(A,data_mask(n3,i));
    x = max(x,0.01);
    x = x./sum(x);
    C1C = x(3);
    C2C = x(1)/(1-x(3));  
    
    %C1A = min(max(C1A,1e-2),0.99);
    %C2A = min(max(C2A,1e-2),0.99);
    %C1B = min(max(C1B,1e-2),0.99);
    %C2B = min(max(C2B,1e-2),0.99);

    theta0a_mask(:,i) = [C1A;C2A;C3A;S0_mask_A(i);C1B;C2B;C3B;S0_mask_B(i);C1C;C2C;C3C;S0_mask_C(i);theta1;phi1;theta2;phi2];
    theta0b_mask(:,i) = [C1A;0.5;0.5;S0_mask_A(i);C1B;0.5;0.5;S0_mask_B(i);C1C;0.5;0.5;S0_mask_B(i);theta1;phi1;theta2;phi2];
    
end

%% Store into pre-allocated matrices
theta0a(:,mask(:)) = theta0a_mask;
theta0b(:,mask(:)) = theta0b_mask;

theta0_cell{1} = theta0a;
theta0_cell{2} = theta0b;

end

function [S0,DT] = fitSingleTensor(data,Q) 
    
    % The logarithm of zero is inf. To prevent these 'infs' all
    % datapoints<eps are set to eps.
    data = max(data,eps);
    
    % A log-linear fit is very sensitive to small values. To prevent
    % misfits due to outliers, each datapoint should at least have a value
    % of 1/100th the mean value in each voxel.
    data = max(repmat(mean(data,1)/100,[size(data,1) 1]),data);
   
    % Prepare the log-linear fit
    y = log(data);
    A = [Q,ones(size(Q,1),1)];
    
    % Heuristically set weights are used to get solutions that should be a
    % bit closer to the true least-squares solution
    weights = exp(-abs(Q(:,1)+Q(:,2)+Q(:,3))*2e-3);
    yw = y.*repmat(weights,[1 size(y,2)]);
    Aw = A.*repmat(weights,[1 7]);
    Aw_inv = inv(Aw'*Aw)*Aw';
    x = Aw_inv*yw;
    
    % Compute tensor and S0
    DT = x(1:6,:);
    S0 = exp(x(7,:));
    
    % If one or more b0 images are present, S0 is computed from the mean
    b0indices = (sum(abs(Q),2)==0);
    if sum(b0indices)>=1
        S0 = mean(data(b0indices,:),1);
    end
    
end