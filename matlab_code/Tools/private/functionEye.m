function [x,g,h] = functionEye(x)
% [x,g,h] = functionEye(x)
% Jx = J * x = x
% Just pass x from input to output, set g and h to empty.
% Dummy function that helps to most efficiently evaluate the special 'eye' 
% function argument of fit_MRI. 
% (i.e. when the string 'eye' is passed as first argument to fit_MRI)
%
% Created by Dirk Poot, Erasmus MC, 2-2-2012
g =[];h=[];

