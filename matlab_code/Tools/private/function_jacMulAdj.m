function [Jx] = function_jacMulAdj( J, x)
% [Jx] = function_jacMulAdj( J, x)
% Function that multiplies adjoint jacobian of function with x
%
% J = nimg x  ntraces  x  nparameters
% x = nimg x ntraces
%
% Jx = nparameters x ntraces.
% Created by Dirk Poot, Erasmus MC, 2-2-2012

Jx = permute( sum( bsxfun(@times, J  , x  ) , 1) , [3 2 1]);
