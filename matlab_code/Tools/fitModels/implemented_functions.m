function [listout] = implemented_functions()
% [list] = implemented_functions()
% List is 2 column cell array that specifies the mfile names
%
% Created by Dirk Poot, 25-11-2013, TUDelft

persistent list

if isempty(list)
% function_handle ,  functionstring ,  coefficients , independent_variables, 
% example: @sin   ,   'a*x^2+b*x+c' , {'a','b','c'} ,     {'x'}
list = {
    @T2pred , 'A*exp(-TE*R2)' , {'A','R2'} , {'TE'} ;
    @predict_FAVT1T2 ,'abs( (A * sin( alpha ) * ( 1 - exp(- TR * (R1/1000)) / (1 - cos( alpha ) *exp(- TR * (R1/1000))))*exp(-TE * R2/1000)) )', {'R1','A','R2'},{'alpha','TR','TE'}
};
end;

listout = list;