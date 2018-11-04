function [y, dy, d2ydmudmu] = lognormpdf(x,mu,sigma, derivativeselector)
%LOGNORMPDF Normal probability density function (pdf).
%   [Y , dydmu, d2ydmudmu] = LOGNORMPDF(X,MU,SIGMA [,derivativeselector]) returns the logarithm of the pdf of the
%   normal distribution with mean MU and standard deviation SIGMA,
%   evaluated at the values in X. 
%   The size of Y is the common size of the input arguments.  A scalar
%   input functions as a constant matrix of the same size as the other
%   inputs.
%   Optionally the second and third output argument return the first and second derivative of 
%   Y with respect to the 'mu' argument, respectively.
%   derivativeselector : default = [false true false]; should be 3 element logical vector.
%       The first (and second) derivatives are calculated for the selected parameters ['x' 'mu' 'sigma']
%       if derivative wrt to multiple variables is requested, the individual derivatives are concatenated 
%       in dimensions ndims(y)+1.
%
%   Default values for MU and SIGMA are 0 and 1 respectively.
%
%   See also NORMCDF, NORMFIT, NORMINV, NORMLIKE, NORMRND, NORMSTAT.

if nargin<1
    error('stats:normpdf:TooFewInputs','Input argument X is undefined.');
end
if nargin < 2 || isempty(mu)
    mu = 0;
end
if nargin < 3 || isempty(sigma)
    sigma = 1;
end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;

try
%     y = exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);
    resdsigma = bsxfun(@rdivide, (x - mu), sigma);
    y = bsxfun(@minus, (-0.5 * resdsigma.^2) ,  log(sqrt(2*pi) .* sigma));
    if nargout>=2
        if nargin<4 || derivativeselector(2)
            dy = bsxfun(@rdivide, resdsigma,  sigma) ; % (x-mu)/sigma^2
        else
            dy = zeros([size(x) 0]);
        end;
        if nargin>=4 
            if derivativeselector(3)
                dy = cat(ndims(y)+1, dy, bsxfun(@minus, bsxfun(@rdivide, resdsigma.^2, sigma) , 1./sigma) ); % (x-mu)^2/sigma^3 - 1/sigma
            end;
            if derivativeselector(1)
                dy = cat(ndims(y)+1, bsxfun(@rdivide, resdsigma,  -sigma), dy); % -(x - mu) ./ sigma.^2
            end;
        end;
        if nargout>=3
            if nargin<4 || derivativeselector(2)
                d2ydmudmu = bsxfun(@rdivide, -ones(size(y)), sigma.^2) ; % -1/sigma^2
            else
                d2ydmudmu = [];
            end;
            if nargin>=4 
                if derivativeselector(1)
                    d2ydxdx = bsxfun(@rdivide, -ones(size(y)), sigma.^2) ;% -1/sigma^2
                else
                    d2ydxdx = [];
                end;
                if derivativeselector(1) && derivativeselector(2)
                    d2ydxdmu = bsxfun(@rdivide, ones(size(y)), sigma.^2); % 1/sigma^2 ;
                else
                    d2ydxdmu = [];
                end;
                
                if derivativeselector(3)
                    d2ydsigmadsigma = bsxfun(@plus, - 3*bsxfun(@rdivide, resdsigma.^2, sigma.^2) , 1./sigma.^2);  %- 3*(x - mu).^2 ./ sigma.^4 + 1./sigma.^2
                else
                    d2ydsigmadsigma = [];
                end;
                    
                if derivativeselector(1) && derivativeselector(3)
                    d2ydsigmadx =   bsxfun(@rdivide, resdsigma, 0.5*sigma.^2);  % 2*(x - mu) ./ sigma.^3
                else
                    d2ydsigmadx = [];
                end;

                if derivativeselector(2) && derivativeselector(3)
                    d2ydsigmadmu =  bsxfun(@rdivide, resdsigma, -0.5*sigma.^2); % - 2*(x - mu) ./ sigma.^3
                else
                    d2ydsigmadmu = [];
                end;
                d2ydmudmu = cat(ndims(y)+1, d2ydxdx, d2ydxdmu, d2ydmudmu, d2ydsigmadx , d2ydsigmadmu, d2ydsigmadsigma);
            end;
        end;
    end;
catch
    error('stats:normpdf:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end
