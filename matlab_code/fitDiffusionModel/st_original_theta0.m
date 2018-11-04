function theta0_cell = st_original_theta0(S0,DT,mask,opts)
%% st_original_theta0

%% Read relevant fields from opts
d_max = opts.d_max;

%% Pre-allocate space for parameter matrix
nPar = 7;
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
    EW = max(EW,0.01e-3);
    EW = min(EW,2.99e-3);
       
    [alpha1,alpha2,alpha3] = R2sph(EV);
    
    % Compute C1, C2, C3
    C1 = EW(1)/d_max;
    C2 = EW(2)/d_max;
    C3 = EW(3)/d_max;
    C1 = min(max(C1,1e-2),0.99);
    C2 = min(max(C2,1e-2),0.99);
    C3 = min(max(C3,1e-2),0.99);
    
    theta0a_mask(:,i) = [alpha1,alpha2,alpha3,C1,C2,C3,S0_mask(i)];
    theta0b_mask(:,i) = [0,0,0,0.5,0.5,0.5,S0_mask(i)];
    
end

%% Store into pre-allocated matrices
theta0a(:,mask(:)) = theta0a_mask;
theta0b(:,mask(:)) = theta0b_mask;

theta0_cell{1} = theta0a;
theta0_cell{2} = theta0b;

end