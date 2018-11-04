function theta0_cell = stv2_theta0(S0,DT,mask,opts)
%% stv2_theta0

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

for i=1:size(theta0a_mask,2)
    
    [~,~,EV,EW] = computeStatisticsFromTensor(DT_mask(:,i));
    
    % Put EW in physically sensible range
    EW = max(EW,0.1e-4);
    EW = min(EW,d_max);
    
    % Approximate new eigenwaarden
    lambda_par = EW(1);
    lambda_perp = 0.5*(EW(2)+EW(3));
    
    C1 = lambda_par/d_max;
    C1 = min(max(C1,1e-2),0.99);
    C2 = lambda_perp/(C1*d_max);
    C2 = min(max(C2,1e-2),0.99);
    
    % Initialization I from EV1
    EV1 = EV(:,1);
    EV1 = EV1./norm(EV1);
    
    [phi1,theta1] = cart2sph(EV1(1),EV1(2),EV1(3));
    theta1 = pi/2-theta1;
    theta0a_mask(:,i) = [C1;C2;theta1;phi1;S0_mask(i)];
    
    % Initialization II from EV2
    EV2 = EV(:,2);
    EV2 = EV2./norm(EV2);
    
    [phi2,theta2] = cart2sph(EV2(1),EV2(2),EV2(3));
    theta2 = pi/2-theta2;
    theta0b_mask(:,i) = [C1;C2;theta2;phi2;S0_mask(i)];
    
end

%% Store into pre-allocated matrices
theta0a(:,mask(:)) = theta0a_mask;
theta0b(:,mask(:)) = theta0b_mask;

theta0_cell{1} = theta0a;
theta0_cell{2} = theta0b;

end