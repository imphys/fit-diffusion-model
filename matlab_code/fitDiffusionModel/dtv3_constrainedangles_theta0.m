function theta0_cell = dtv3_constrainedangles_theta0(S0,DT,mask,opts)
%% dtv3_constrainedangles_theta0

%% Read relevant fields from opts
d_max = opts.d_max;

%% Pre-allocate space for parameter matrix
nPar = 4;
theta0a = zeros([nPar size(mask)]);
theta0b = zeros([nPar size(mask)]);

%% Only compute theta0 inside mask to save time
theta0a_mask = theta0a(:,mask(:));
theta0b_mask = theta0b(:,mask(:));
S0_mask = S0(mask(:));
DT_mask = DT(:,mask(:));

for i=1:size(theta0a_mask,2)
    
    [~,~,~,EW] = computeStatisticsFromTensor(DT_mask(:,i));
    
    % Put EW in physically sensible range
    EW = max(EW,0.1e-4);
    EW = min(EW,d_max);
    
    % Approximate new eigenwaarden
    lambda_par = EW(1)+EW(2)-EW(3);
    lambda_perp = EW(3);
    
    C1 = lambda_par/d_max;
    C1 = min(max(C1,1e-2),0.99);
    C2 = lambda_perp/(C1*d_max);
    C2 = min(max(C2,1e-2),0.99);
    
    f = 0.3;
    theta0a_mask(:,i) = [f;C1;C2;S0_mask(i)];
    f = 0.7;
    theta0b_mask(:,i) = [f;C1;C2;S0_mask(i)];
    
end

%% Store into pre-allocated matrices
theta0a(:,mask(:)) = theta0a_mask;
theta0b(:,mask(:)) = theta0b_mask;

theta0_cell{1} = theta0a;
theta0_cell{2} = theta0b;

end