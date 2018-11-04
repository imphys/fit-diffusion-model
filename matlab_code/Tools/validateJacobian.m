function validateJacobian(fun, x, jacobmulfun)
% validateJacobian( fun, x [, jacobmulfun])
% Function that validates the jacobian of the function 'fun' in x.
% Also supports model functions for fit_MRI.
% This can be used to check for errors in the jacobian computation.
% NOTE: this routine is slow as it numerically evaluates the jacobian.
%       Only use for debugging, and try to minimize the problem size.
%
% Evaluates [f, J] = fun(x)
% and checks if the numerical jacobian of fun around x is close enough 
% to the numerical jacobian J.
% If jacobmulfun is provided, it is used to compute the jacobian explicitly.
%    J = jacobmulfun(J, eye ,1); 
%    J'= jacobmulfun(J, eye ,-1); 
%   Additionally J and J' are compared.
% 
% NOTE:
% This function uses jacobianest from the 'DERIVESTsuite' by John D'Errico:
% http://www.mathworks.com/matlabcentral/fileexchange/13490-automatic-numerical-differentiation
%
% Created by Dirk Poot, Erasmus MC, 2-2-2012
[f, J] = fun(x);
if nargin>=3 && ~isempty(jacobmulfun)
    % check jacobmulfun:
    J1 = jacobmulfun(J, eye(numel(x)), 1);
    JT = jacobmulfun(J, eye(numel(f)), -1);
    if max(max(abs(J1-JT')./abs(J1+JT')))>1e-8
        error('J and J'' of Jacobmulfun seem to differ.');
    end;
    r = randn(numel(x),1);
    rH = jacobmulfun(J, r, 0);
    rHs = JT * (J1 * r);
    if max(max(abs(rH-rHs),[],2)./max(abs(rH+rHs),[],2))>1e-6
        error('J''*J of Jacobmulfun seem to differ from multiplying with J'' and J separately.');
    end;
    J = J1;
end;
% Using jacobianest from the 'DERIVESTsuite' by John D'Errico
% http://www.mathworks.com/matlabcentral/fileexchange/13490-automatic-numerical-differentiation
[Jest , errest] = jacobianest(fun, x);
if size(J,3)==size(x,1) && size(J,2)==size(x,2) && size(J,1)*size(J,2)==size(Jest,1)
    % Assume multi time series fit (as for fit_MRI fitting functions)
    Jf = zeros(size(J,1)*size(J,2), size(J,2)*size(J,3));
    strowid = 1;
    stcolid = 1;
    for k=1:size(x,2)
        edrowid = strowid + size(J,1)-1;
        edcolid = stcolid + size(x,1)-1;
        Jf(strowid:edrowid, stcolid:edcolid) = squeeze(J(:,k,:));
        strowid = edrowid+1;
        stcolid = edcolid+1;
    end;
    J= Jf;
end;
if any(any( (Jest-J)./errest >2 & abs(J-Jest)./abs(J+Jest)>1e-4 ))
    error('Numerical and analytic jacobians seem to differ.');
end;
disp('The jacobian appears to be computed correctly.');

