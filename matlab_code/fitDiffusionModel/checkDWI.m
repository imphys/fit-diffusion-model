function [noiselvl_out,signal_out,noiselvl_perdwi_out,noiselvl_pervoxel_out,S0_out,DT_out] = checkDWI(fn_dwi,fn_mask,Q)

% Load image header
img_header = load_untouch_header_only(fn_dwi);

% Prepare for loop
number_of_slices = img_header.dime.dim(4);
S0_out = zeros([1 img_header.dime.dim(2:4)]);
DT_out = zeros([6 img_header.dime.dim(2:4)]);
noiselvl_pervoxel_out = zeros([1 img_header.dime.dim(2:4)]);
noiselvl_perdwi_out = zeros(size(Q,1),number_of_slices);

% Iterate over over slice and fill S0, DT, noiselvl_pervoxel and noiselvl_perdwi
for slice_idx=1:number_of_slices
    
    % Load single slice from data and mask
    data = load_untouch_nii(fn_dwi,[],[],[],[],[],slice_idx);
    data = permute(data.img,[4 1 2 3]);
    mask = load_untouch_nii(fn_mask,[],[],[],[],[],slice_idx);
    mask = mask.img>0;
    
    % Do some reshaping so different dimensional signals can be treated similar
    dim_data_in = size(data);
    nDWI = dim_data_in(1);
    data = double(reshape(data,nDWI,numel(data)/nDWI));
    
    % The logarithm of zero is inf. To prevent these 'infs' all
    % datapoints<eps are set to eps.
    data_orig = data;
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
    D = x(1:6,:);
    S0 = exp(x(7,:));
    
    % If one or more b0 images are present, S0 is computed from the mean
    b0indices = (sum(abs(Q),2)==0);
    if sum(b0indices)>=1
        S0 = mean(data(b0indices,:),1);
    end
    
    % Predict Sq and noise level
    S_est = repmat(S0,nDWI,1).*exp(Q*D);   
    resnorm_perdwi = sum((S_est(:,mask)-data_orig(:,mask)).^2,2);
    N = max(sum(mask(:)),1);
    noiselvl_perdwi = sqrt(resnorm_perdwi/(N-7*N/nDWI));
    resnorm_pervoxel = sum((S_est-data_orig).^2,1);
    noiselvl_pervoxel = sqrt(resnorm_pervoxel/(nDWI-7));
    
    % Store in big matrix
    noiselvl_perdwi_out(:,slice_idx) = noiselvl_perdwi;
    noiselvl_pervoxel_out(:,:,:,slice_idx) = reshape(noiselvl_pervoxel,[1 dim_data_in(2:end)]);
    DT_out(:,:,:,slice_idx) = reshape(D,[6 dim_data_in(2:end)]);
    S0_out(:,:,:,slice_idx) = reshape(S0,[1 dim_data_in(2:end)]);

end

% Load mask
mask = load_untouch_nii(fn_mask);
mask = mask.img>0;

% Construct output
noiselvl_out = median(noiselvl_pervoxel_out(mask));
noiselvl_perdwi_out = median(noiselvl_perdwi_out,2);
signal_out = median(S0_out(mask));

end