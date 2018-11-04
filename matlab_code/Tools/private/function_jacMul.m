function [Jx] = function_jacMul( J, x)
% function [Jx] = function_jacMul( J, x)
% Function that multiplies jacobian of function with x
% J = nimg  x  ntraces  x  nparameters
% x = nparameters x ntraces
% 
% Jx = nimg x ntraces
%
% Created by Dirk Poot, Erasmus MC, 2-2-2012

Jx = sum( bsxfun(@times, J  , permute( x , [3 2 1]) ) , 3);
