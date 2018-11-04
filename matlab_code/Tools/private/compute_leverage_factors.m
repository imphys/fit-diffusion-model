function [ h ] = compute_leverage_factors( theta, data, opts )
% [ h ] = compute_leverage_factors( theta, data, opts );
% 
% Helper function that computes the leverage factors when 
% fit_MRI is called with 'computeLeverage' set to true.
%
% The leverage factor is the factor by which a small change in the data
% is absorbed by theta and thus also the factor by which the predicted
% value changes due to a change in data.
%  if h = 
%   0 : prediction not affected by that data element
%   1 : any change in data is completely absorbed in the fit. 
%
% Created by Dirk Poot, 10-12-2012, Erasmus MR.


if opts.numPDFoptpar~=0 || ~isempty(opts.project)
    error('(currently) cannot compute leverage factors when estimating noise level or wen projecting');
end;
szTheta = size(theta);
if isempty(opts.fields)
    predfun =opts.function;
end;
if opts.numPDFoptpar==0
    noiseLevel = opts.noiseLevel;
    logpdffun =  @(data, A) opts.logPDFfun(data, A, noiseLevel, [true true false]);  % derivative w.r.t. data and A
else
    error('should not reach.');
end;
logprior = opts.parameterPrior.fun;
h = zeros(size(data));
for voxid = 1 : prod(szTheta(2:end));
    if ~isempty(opts.fields)
        fieldsel = opts.fields( : , voxid );
        predfun = make1arg_anonfun( opts.function, fieldsel );
    end;
    thetasel = theta( : , voxid );
    datasel = data( : , voxid );
    h(:, voxid) = compute_leverage_factor_block( thetasel, datasel, predfun, logpdffun , logprior);
end;
    

function [h] = compute_leverage_factor_block( theta, data, predfun, logpdffun , logprior)
%  h = factor with which a change in data appears in a change in predicted
%      value
% So : 
% h(i) = ( predfun( arg min cost_function( theta , data + delta_i * eps) ) -
%          predfun( theta_hat ) )/eps
%
% approximate by taylor expansion around theta_hat : 
% Assumptions: 
%    - assuming second derivative of function is zero (/negligible) 
%    - gradient at theta_hat is zero (=> it is an optimimum). If not true, 
%          arg min cost_function( theta , data + delta_i * eps) - gradient cost_function( theta_hat ) 
%       is used in the 'predfun' in the first line of the equation above.
%
%  C( theta , data ) =  logp( data , fun( theta ) )
%  = logp + dlogp_ddata * delta_data + dlogp_dA * dfun_dtheta * delta_theta
%         + .5*d2logp_ddataddata * delta_data.^2 
%         + d2logp_ddatadA * delta_data * dfun_dtheta * delta_theta
%         + .5*d2logp_dAA * (dfun_dtheta * delta_theta).^2
% delta_theta for a given delta_data ; assuming in optimum, so dlogp_dA * dfun_dtheta ==0
%  => Solve D[C,delta_theta]==0
% => 0 == d2logp_ddatadA * delta_data * dfun_dtheta + 2 *.5*d2logp_dAA * dfun_dtheta * dfun_dtheta * delta_theta
% => delta_theta = -(d2logp_ddatadA * delta_data * dfun_dtheta) * 
%                    inv( d2logp_dAA * dfun_dtheta * dfun_dtheta )

[f_pred, f_grad] = predfun(theta);
f_grad = reshape( f_grad , [], size(theta,1));
[L, L_G, L_H] = logpdffun( data, f_pred );

H  = f_grad' * bsxfun(@times, L_H(:,3), f_grad );
G = bsxfun(@times, L_H(:,2), f_grad )';
if ~isempty(logprior)
    [p_L, p_G, p_H] = logprior( theta );
    H = H + p_H;
end;
delta_theta = -(H\G);
h = sum(f_grad.*delta_theta' ,2);