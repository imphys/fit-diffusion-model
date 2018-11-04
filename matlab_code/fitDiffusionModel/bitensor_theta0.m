function theta0_cell = bitensor_theta0(S0,DT,mask,opts)
%% bitensor_theta0

%% Read relevant fields from opts
d_max = opts.d_max;
d_iso = opts.d_iso;
mdc = opts.MD;

%% Pre-allocate space for parameter matrix
nPar = 8;
theta0a = zeros([nPar size(mask)]);
theta0b = zeros([nPar size(mask)]);

%% Only compute theta0 inside mask to save time
theta0a_mask = theta0a(:,mask(:));
theta0b_mask = theta0b(:,mask(:));
S0_mask = S0(mask(:));
DT_mask = DT(:,mask(:));

for i=1:size(theta0a_mask,2)
    
    [~,~,EV,EW] = computeStatisticsFromTensor(DT_mask(:,i));
    
    [alpha1,alpha2,alpha3] = R2sph(EV);
    
    % Put EW in physically sensible range
    MD = mean(EW);
    EW = max(EW,0.1e-4);
    EW = min(EW,d_max);
    
    %% Approximate fiso and f
    fiso = (MD-mdc)/(d_iso-mdc);
    fiso = min(max(fiso,0.05),0.95);
    f = 1-fiso;
    
    %% Compute lambda123
    C1 = EW(1)/d_max;
    C2 = EW(2)/d_max;
    C3 = EW(3)/d_max;
    C1 = min(max(C1,1e-2),0.99); % Do not initialize too close to constraint
    C2 = min(max(C2,1e-2),0.99); % Do not initialize too close to constraint
    C3 = min(max(C3,1e-2),0.99); % Do not initialize too close to constraint
    
    %% Store in parameter vector
    theta0a_mask(:,i) = [alpha1;alpha2;alpha3;C1;C2;C3;f;S0_mask(i)];
    theta0b_mask(:,i) = [0;0;0;1/2;1/2;1/2;0.1;S0_mask(i)];

end

%% Store into pre-allocated matrices
theta0a(:,mask(:)) = theta0a_mask;
theta0b(:,mask(:)) = theta0b_mask;

theta0_cell{1} = theta0a;
theta0_cell{2} = theta0b;

end