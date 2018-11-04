function apodize(varargin)
% apodize applies an apodization filter to the diffusion-weighted data. It
% assumes that the frequency and phase encoding direction respectively
% correspond to the x and y direction of the matrix loaded by
% load_untouch_nii.
%
% Note:
% 1. It expects niftitools from Jimmy Shen to be in the matlab-path.
% 2. It expects fit_mri tools from Dirk Poot to be in the matlab-path. It
% may be necessary to compile one or more mex-files in fit_mri.
% 3. fitDiffusionModel was implemented in matlab2013b and has not been
% thoroughly tested on other versions.
%
% created by Joor Arkesteijn (g.a.m.arkesteijn@tudelft.nl), Quantitative Imaging Group, Delft.

%% Required options
opts.in = '';
opts.out = '';
opts.w = '';
opts.kx = '';
opts.ky = '';

%% Optional options
opts.sigmaz = '0';

%% If deployed, remove the minus signs in front of the option fields (odd fields)
if isdeployed
    for i=1:2:numel(varargin)
        if strcmpi(varargin{i}(1),'-')
            varargin{i}=varargin{i}(2:end);
        end
    end
end

%% Parse optionvaluepairs
opts = parse_defaults_optionvaluepairs( opts, varargin{:});

%% Output dir
[output_dir,~,~] = fileparts(opts.out);

%% Convert chars to doubles or booleans
opts.kx = ceil(str2double(opts.kx)/2);
opts.ky = ceil(str2double(opts.ky)/2);
opts.sigmaz = str2double(opts.sigmaz);

%% Load image
nii = load_untouch_nii(opts.in);
dim = size(nii.img);
new_nii = double(nii.img);

%% Create a window with dimensions 'zeropad_dimensions' for apodization
window_x = create1Dwindow(dim(1),opts.kx,opts.w);
window_x = repmat(window_x,[1 dim(2)]);
window_y = create1Dwindow(dim(2),opts.ky,opts.w);
window_y = repmat(reshape(window_y,1,[]),[dim(1) 1]);
window = window_x.*window_y;
  
for i=1:size(new_nii,4)
    
    for j=1:size(new_nii,3)
        
        %% Load the i-th diffusion weighted image
        im_temp = double(nii.img(:,:,j,i));
        
        %% Compute the Fourier transform
        F_temp = fftn(im_temp);
        F_temp = fftshift(F_temp);
        
        %% Check if phase and frequency encoding direction are matched
        if j==floor(size(new_nii,3)/2) && i==1
            
            im1 = log(abs(F_temp)+1);
            im2 = log(abs(window)+1);
            im3 = log(abs(F_temp.*window)+1);
            
            imwrite(im1./max(im1(:)),[output_dir filesep() 'apodize_orig.png']);
            imwrite(im2./max(im2(:)),[output_dir filesep() 'apodize_filter.png']);
            imwrite(im3./max(im3(:)),[output_dir filesep() 'apodize_filtered.png']);

        end
        
        %% Apply the apodization window
        F_temp = F_temp.*window;
        F_temp = ifftshift(F_temp)/numel(F_temp)*numel(F_temp);
      
        %% Take the inverse Fourier transform and store
        im_temp = ifftn(F_temp);
        im_temp = sqrt(real(im_temp).^2+imag(im_temp).^2);
        new_nii(:,:,j,i) = im_temp;
    end
    
    %% Apply Gaussian blur in slice selection direction
    if opts.sigmaz>0
       new_nii(:,:,:,i) = gaussblur(new_nii(:,:,:,i),[0 0 opts.sigmaz]); 
    end
    
end

nii.img = new_nii;

%% Adjust header and save to disk
save_untouch_nii(nii,opts.out)
end

function [window,kernel] = create1Dwindow(length,kmax,type)

midpoint = floor(length/2)+1;
k = (1:length)-midpoint;
arg = abs(k)/kmax;

if strcmpi(type,'lukosz')
    window = cos((pi*arg)./(arg+1));
    window(arg>1) = 0;
elseif strcmpi(type,'boxcar')
    window = ones(length,1);
    window(arg>1) = 0;
elseif strcmpi(type,'triangular')
    window = 1-arg;
    window(arg>1) = 0;
elseif strcmpi(type,'gauss')
    sigma_i = kmax; %% alright one exception, gaussian blur has no kmax, unit of sigma is voxels in image space
    sigma_k = length/(2*pi*sigma_i);
    window = exp(-k.^2/(2*sigma_k^2));
elseif strcmpi(type,'staircase')
    window = zeros(length,1);
    for i=1:kmax;
        M = i;
        window(abs(k)<kmax/M) = cos(pi/(M+1));
    end
end

window = window(:);
kernel = fftshift(ifft(ifftshift(window)));

end
