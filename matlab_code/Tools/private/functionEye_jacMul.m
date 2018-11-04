function x = functionEye_jacMul(dum, x)
% x = functionEye_jacMul(dum, x)
% just pass x from input to output, ignore dum.
% Dummy function that helps to most efficiently evaluate the special 'eye' 
% function argument of fit_MRI. 
% (i.e. when the string 'eye' is passed as first argument to fit_MRI)
%
% Created by Dirk Poot, Erasmus MC, 2-2-2012
