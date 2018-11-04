function [y, dy, d2ydmudmu] = logRobustNormpdf(x,mu,sigma, derivativeselector, offset)
%LOGNORMPDF logarithm of a robustified normal probability density function (pdf).
%   [Y , dydmu, d2ydmudmu] = logRobustNormpdf(X,MU,SIGMA [,derivativeselector [,offset=16]]) 
%   returns the logarithm of the pdf of a robustified 
%   normal distribution with mean MU and standard deviation SIGMA,
%   evaluated at the values in X. 
%   For up to abs(x-mu)/sigma < 1 it is almost exactly the normal distribution.
%   for (x-mu)/sigma >>1 the logRobustNormpdf increases as the square root of that value.
%       Y = 16*(sqrt(sqrt(2*(X-MU)/SIGMA.^2+16))-2) + f( SIGMA );
%   f( SIGMA ) approximately normalizes the PDF. It differs from a fully normalized PDF 
%   by a (small) constant g( offset ).
%
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
% Created by Dirk Poot, Erasmus MC, 3-2-2012

if nargin<5 || isempty(offset)
    offset = 16; % offset that determines the range in which the pdf corresponds closely to the normal pdf.
end;
sqoffset = sqrt(sqrt(offset));
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

    resdsigma = bsxfun(@rdivide, (x - mu), sigma);
    resoffssq = (resdsigma.^2+offset);
    resoffssqrr = sqrt(sqrt(resoffssq));
    y = bsxfun(@minus, ((-offset*2/sqoffset) * (resoffssqrr-sqoffset)) ,  log(sqrt(2*pi) .* sigma));
    if nargout>=2
        if nargin<4 || isempty(derivativeselector)
            derivativeselector = [false true false];
        end;
        if derivativeselector(2)
            % dy/dmu = (-offset*2/sqoffset) * (1/4)* 2*resdsigma * -1/sigma /(resdsigma.^2+offset)^(3/4)
            %        = offset/sqoffset *resdsigma/sigma * (resdsigma.^2+offset)^(1/4) / (resdsigma.^2+offset)
            dy = bsxfun(@rdivide, resdsigma.*resoffssqrr./resoffssq,  sigma .* (sqoffset./offset)) ; % (x-mu)/sigma^2
        else
            dy = zeros([size(resdsigma) 0]);
        end;
        if derivativeselector(3)
            % dy/dsigma = 1/sigma + (-offset*2/sqoffset) * 1/4 *1/resoffssq^(3/4) * 2 * resdsigma * -1*(x-mu)/sigma^2
            %           = 1/sigma + (offset/sqoffset) *resoffssqrr/resoffssq * resdsigma^2/sigma 
            dy = cat(ndims(y)+1, dy, bsxfun(@minus, bsxfun(@rdivide, resdsigma.^2.*resoffssqrr./resoffssq, sigma*(sqoffset./offset)) , 1./sigma) ); 
        end;
        if derivativeselector(1)
            %dy/dx = -dy/dmu
            dy = cat(ndims(y)+1, bsxfun(@rdivide, resdsigma.*resoffssqrr./resoffssq,  sigma .* (-sqoffset./offset)) , dy); % -(x - mu) ./ sigma.^2
        end;
        if nargout>=3 
            error('cant yet compute hessian');
            if derivativeselector(2)
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
