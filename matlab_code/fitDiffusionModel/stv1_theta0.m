function theta0_cell = stv1_theta0(S0,DT,mask,opts)
%% stv1_theta0

%% Read relevant fields from opts
d_max = opts.d_max;
d_iso = opts.d_iso;
mdc = opts.MD;

%% Pre-allocate space for parameter matrix
nPar = 6;
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
    MD = mean(EW);
    EW = max(EW,0.1e-4);
    EW = min(EW,d_max);
    
    % Approximate new eigenwaarden
    lambda_par = EW(1);
    lambda_perp = (EW(2)+EW(3))/2;
    
    C1 = (MD-mdc)/(d_iso-mdc); % fiso    
    C1 = min(max(C1,0.05),0.95);
    C3 = lambda_par/d_max;
    C3 = min(max(C3,1e-2),0.99);
    C4 = lambda_perp/(C3*d_max);
    C4 = min(max(C4,1e-2),0.99);
    
    % Initialization II (perpendicular split)
    EV1 = EV(:,1);
    [phi1,theta1] = cart2sph(EV1(1),EV1(2),EV1(3));
    theta1 = pi/2-theta1;
        
    theta0a_mask(:,i) = [C1;C3;C4;theta1;phi1;S0_mask(i)];
    theta0b_mask(:,i) = [0.1;0.5;0.95;0;0;S0_mask(i)];
    
end
    
%% Store into pre-allocated matrices
theta0a(:,mask(:)) = theta0a_mask;
theta0b(:,mask(:)) = theta0b_mask;
theta0_cell{1} = theta0a;
theta0_cell{2} = theta0b;

end