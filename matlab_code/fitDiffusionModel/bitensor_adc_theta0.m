function theta0_cell = bitensor_adc_theta0(S0,DT,mask,opts)
%% bitensor_adc_theta0
% theta1: alpha1
% theta2: alpha2
% theta3: alpha3
% theta4: C1 on domain [0,1] (lambda2 = C1*lambda1)
% theta5: C2 on domain [0,1] (lambda3 = C2*lambda1)
% theta6: f
% theta7: S0

%% Read relevant fields from opts
d_iso = opts.d_iso;
AD = opts.AD;

%% Pre-allocate space for parameter matrix
nPar = 7;
theta0a = zeros([nPar size(mask)]);
theta0b = zeros([nPar size(mask)]);

%% Only compute theta0 inside mask to save time
theta0a_mask = theta0a(:,mask(:));
theta0b_mask = theta0b(:,mask(:));
S0_mask = S0(mask(:));
DT_mask = DT(:,mask(:));

%% Setting for computing theta0
b = 1000; % Used to compute intial guess for fiso
dmin = 0.01e-3; % Sensible range for diffusivity
dmax = 2.99e-3; % Sensible range for diffusivity

for i=1:size(theta0a_mask,2)
    
    [~,~,EV,EW] = computeStatisticsFromTensor(DT_mask(:,i));
    
    [alpha1,alpha2,alpha3] = R2sph(EV);
    
    % Put EW in physically sensible range
    EW = max(EW,dmin);
    EW = min(EW,dmax);
    
    %% Approximate fiso and f
    fiso = (exp(-b*EW(1))-exp(-b*AD))/(exp(-b*d_iso)-exp(-b*AD));
    f = 1-fiso;
    
    if fiso>0.95 % in CSF
        fiso = 0.95;
        f = 1-fiso;
    elseif fiso<-1  % weird voxels
        fiso = -1;
        f = 1-fiso;
    end
    
    %% Compute lambda23
    lambda2 = log((exp(-b*EW(2))-fiso*exp(-b*d_iso))/(1-fiso))/-b;
    lambda3 = log((exp(-b*EW(3))-fiso*exp(-b*d_iso))/(1-fiso))/-b;
    
    %% Scale lambda23 between dmin and AD
    lambda2 = min(max(lambda2,dmin),AD);
    lambda3 = min(max(lambda3,dmin),AD);
    
    %% Perform asin mapping
    C1 = lambda2/AD;
    C2 = lambda3/AD;
    C1 = min(max(C1,1e-2),0.99); % Do not initialize too close to constraint
    C2 = min(max(C2,1e-2),0.99); % Do not initialize too close to constraint
    
    %% Store in parameter vector
    theta0a_mask(:,i) = [alpha1;alpha2;alpha3;C1;C2;f;S0_mask(i)];
    theta0b_mask(:,i) = [0;0;0;1/3;1/2;0.1;S0_mask(i)];

end

%% Store into pre-allocated matrices
theta0a(:,mask(:)) = theta0a_mask;
theta0b(:,mask(:)) = theta0b_mask;

theta0_cell{1} = theta0a;
theta0_cell{2} = theta0b;

end