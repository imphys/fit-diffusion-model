function [P,dP,hP] = logricepdf_logsigma( x, A, logsigma, selGrHess)
% [P,dP,hP] = logricepdf_logsigma( x, A, logsigma, selGrHess)
% Same as logrice pdf function (which this function calls)
% but with a different parameterisation of sigma.
%
% logricepdf_logsigma( x, a , log(sigma) ) == logricepdf( x , a, sigma)
%
% Using the logarithm of sigma makes regularisation better. 
% (relative changes in sigma should be penalised).
%
% Use in fit_MRI by specifying options:
%       'logPDFfun', @logricepdf_logsigma
%       'imageType','custom'
%
% Created by Dirk Poot, Erasmus MC, 6-12-2012

sigma = exp(logsigma);
if nargout==1
    [P] = logricepdf(x, A, sigma, selGrHess);
elseif nargout>=2
    if nargout>=3
        [P,dP,hP] = logricepdf(x, A, sigma, selGrHess);
    else
        [P,dP] = logricepdf(x, A, sigma, selGrHess);
    end;
    % d/dx f( exp( x) ) = exp(x) *( (d/da f(a)) | a=exp(x))
    if selGrHess(3)
        repf = size(P)./size(sigma);
        if all(repf==1)
            sigma = sigma(:);
        else
            sigma = reshape( repmat(sigma, repf) ,[],1);
        end;
        dP(:,end) = dP(:,end) .* sigma;
        if nargout>=3
            % hessian terms w.r.t. sigma are all grouped last
            hP(:,end-nnz(selGrHess)+1:end) = bsxfun(@times, hP(:,end-nnz(selGrHess)+1:end) , sigma);
            hP(:,end) = hP(:,end) .* sigma + dP(:,end);

        end;   
    end;
end;

function test
%%
x=[1.000000000000000
   0.339595526456833
   0.339595526456833
   0.339595526456833
   0.339595526456833
   0.339595526456833
   0.339595526456833];
A=[ 0.999992345098774
   0.339573016633061
   0.339573120366108
   0.339572796289964
   0.339573040975058
   0.339573256940345
   0.339572880645428];
logsigma = -10.214;
[P,dP,hP] = logricepdf_logsigma( x, A, logsigma, [false true true]);
d = .001;
dPnumA = (logricepdf_logsigma( x, A+d, logsigma, [false true true]) - logricepdf_logsigma( x, A-d, logsigma, [false true true]))/(2*d);
dPnums = (logricepdf_logsigma( x, A, logsigma+d, [false true true]) - logricepdf_logsigma( x, A, logsigma-d, [false true true]))/(2*d);

dP(:,2)-dPnums
