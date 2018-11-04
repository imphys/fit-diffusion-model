function [f, g, h] = lognormalVec(par, mu, neginvsigma, hessIJlin)
% [f, g, h] = lognormalVec(par, mu, neginvsigma)
% Computes (unnormalized) log likelihood of vector normal distribution
% and derivatives w.r.t. par
% Formaly evaluates
%  f = 0.5* trace( (par-mu)' * neginvsigma * (par-mu) )
%
% INPUTS:
%  par : column vector, or ND-array with column vectors
%  mu  : 'mean' value of the normal distribution
%        expanded with bsxfun to size par
%  neginvsigma : -inv(Sigma), where Sigma is the covariance matrix of the normal distribution.
%  hessIJlin   : column vector with the linear index of the elements of 
%                the per-vector hessian that are returned:
%                hi = hess_f( par(:,i) )
%                h(:,i) = hi( hessIJlin );
%
% OUTPUTS:
%   f : scalar function value
%   g : gradient of f w.r.t. par
%   h : compactified (with hessIJlin) hessian of f w.r.t. par.
%
% Created by Dirk Poot, Erasmus MC.
% Extracted from fit_MRI: 27-1-2012

r = bsxfun(@minus, par, mu);
g = neginvsigma*r;
f = sum(r.*g,1);
f = .5*sum(f(:));
if nargout>=3
    szpar = size(par);
    h = neginvsigma(hessIJlin)*ones(1,prod(szpar(2:end)));
end;
