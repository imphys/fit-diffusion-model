function theta0_cell = balland1stick_theta0(S0,DT,data,Q,mask,opts)
%% balland1stick_theta0

%% Read relevant fields from opts
d_max = opts.d_max;

%% Pre-allocate space for parameter matrix
nPar = 5;
theta0a = zeros([nPar size(mask)]);
theta0b = zeros([nPar size(mask)]);

%% Only compute theta0 inside mask to save time
theta0a_mask = theta0a(:,mask(:));
theta0b_mask = theta0b(:,mask(:));
S0_mask = S0(mask(:));
DT_mask = DT(:,mask(:));
data_mask = data(:,mask(:));

for i=1:size(theta0a_mask,2)
    
    [~,~,EV,EW] = computeStatisticsFromTensor(DT_mask(:,i));
    
    % Put EW in physically sensible range
    EW = max(EW,0.1e-4);
    EW = min(EW,d_max);
    
    % Approximate new eigenwaarden
    lambda_par = EW(1)+EW(2)-EW(3);
      
    % Approximate C3
    C3 = lambda_par/d_max;
    C3 = min(max(C3,1e-2),0.99);
    lambda_par = C3*d_max;
    
    % Initialization I
    EV1 = EV(:,1);  
    [phi1,theta1] = cart2sph(EV1(1),EV1(2),EV1(3));
    theta1 = pi/2-theta1;
    
    % Use non-negative least squares to compute the volume fractions
    tensor1 = EV1*EV1'*lambda_par;
    tensor1 = tensor1([1 5 9 4 7 8])';
    ball = [lambda_par;lambda_par;lambda_par;0;0;0];
    tensor1_signal = S0_mask(i)*exp(Q*tensor1);
    ball_signal = S0_mask(i)*exp(Q*ball);
    A = [tensor1_signal,ball_signal];
    x = lsqnonneg(A,data_mask(:,i));
    x = max(x,0.01);
    x = x./sum(x);
    
    theta0a_mask(:,i) = [x(2);C3;theta1;phi1;S0_mask(i)];
    theta0b_mask(:,i) = [0.5;0.5;theta1;phi1;S0_mask(i)];
    
end

%% Store into pre-allocated matrices
theta0a(:,mask(:)) = theta0a_mask;
theta0b(:,mask(:)) = theta0b_mask;

theta0_cell{1} = theta0a;
theta0_cell{2} = theta0b;

end