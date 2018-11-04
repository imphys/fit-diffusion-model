function theta = levenberg_marquardt_fast(fun,data,theta,fields,mask)
% levenberg_marquardt_fast. This function is a vectorized implementation of
% the levenberg marquardt algorithms for solving non-linear least squares
% optimization problems.


%% Optimization settings
e1 = 1e-8; % stopping condition on norm of theta-theta_update
e2 = 1e-8; % stopping condition on norm of RSS-RSS_update
maxiter = 100; %% The number of iterations
lambda = 100;

%% Define useful constants
npar = size(theta,1);
nmri = size(data,1);
dimx = size(data,2);
dimy = size(data,3);
dimz = size(data,4);
voxels_in_mask = sum(mask(:));

%% Reshape data and theta
data = reshape(data,nmri,dimx*dimy*dimz);
theta = reshape(theta,npar,dimx*dimy*dimz);

%% Pre-allocate space for vectors and matrices in the loop
data = data(:,mask);
theta_mask = theta(:,mask);
fields_mask = fields(:,mask);
finished_mask = false(voxels_in_mask,1);
recompute_mask = false(voxels_in_mask,1);
lambda = ones(1,voxels_in_mask)*lambda;
iteration_vector = 1:voxels_in_mask;

%% Precompute f, J and RSS
[f,J] = fun(theta_mask,fields_mask);
f = f-data;
RSS = sum(f.^2,1);
J = permute(J,[1 3 2]);

for k=1:maxiter
    
    %[k,maxiter,sum(~finished_mask(:))]
       
    %% Compute f and J
    [f_t,J_t] = fun(theta_mask(:,recompute_mask),fields_mask(:,recompute_mask));
    f_t = f_t-data(:,recompute_mask);
    f(:,recompute_mask) = f_t;
    RSS(recompute_mask) = sum(f_t.^2,1);
    J(:,:,recompute_mask) = permute(J_t,[1 3 2]);
    
    %% Compute update (unfortunately quite difficult to paralelize)
    theta_mask_update = theta_mask;
    for i=iteration_vector(~finished_mask)
       
        %% Get jacobian and function
        J_t = J(:,:,i);
        f_t = f(:,i);
        A_t = J_t'*J_t;
        
        % Do not update degenerate parameters in the model
        upd_param = diag(A_t)>1e-10*max(A_t(:));
        if ~min(upd_param)
            J_t = J_t(:,upd_param);
            A_t = J_t'*J_t;
        end
                
        %% Compute update
        g_t = J_t'*f_t;
        h_t = (A_t+lambda(i)*diag(diag(A_t))) \ -g_t; 
        %h_t = (A_t+lambda(i)*eye(sum(upd_param(:)))) \ -g_t;
        
        %% Update theta
        theta_mask_update(upd_param,i) = theta_mask(upd_param,i) + h_t;

    end
    
    %% Decide whether update will be accepted
    f_update = fun(theta_mask_update(:,~finished_mask),fields_mask(:,~finished_mask));
    f_update = f_update-data(:,~finished_mask);
    RSS_update = sum(f_update.^2,1);    
    accept_update = RSS_update<RSS(~finished_mask);
    
    %% Compute step sizes to determine convergence
    norm_RSS_step = abs(RSS_update-RSS(~finished_mask));
    norm_theta_step = sqrt(sum((theta_mask(:,~finished_mask)-theta_mask_update(:,~finished_mask)).^2,1));
    
    %% Update theta if accepted    
    theta_mask(:,~finished_mask) = repmat(~accept_update,npar,1).*theta_mask(:,~finished_mask) + repmat(accept_update,npar,1).*theta_mask_update(:,~finished_mask);
    
    %% Multiply lambda by ten if accepted, divide lambda by ten if not accepted
    lambda(~finished_mask) = (~accept_update.*lambda(~finished_mask)*10) + (accept_update.*lambda(~finished_mask)/10);
    
    %% Update the recompute f and J mask
    recompute_mask(~finished_mask) = accept_update;
    
    %% Update the finished mask
    finished_mask(~finished_mask) = accept_update.*((norm_RSS_step<e1) + (norm_theta_step<e2)) | (lambda(~finished_mask)>1e10);   
    
    %% Do not recompute converged voxels
    recompute_mask(finished_mask) = false;
    
    if sum(~finished_mask(:))==0
        break;
    end        
  
    
end

theta(:,mask) = theta_mask;
theta = reshape(theta,npar,dimx,dimy,dimz);


end
