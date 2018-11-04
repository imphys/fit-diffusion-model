function [ Mx ] = mul_default_preconditioner( R, x)
%[ Mx ] = mul_default_preconditioner( R, x)
% Multiplies with the preconditioner created by make_default_preconditioner
% Assumes a cholesky (like) decomposition of H.
% that is : R'*R =(approx)= H
% and R easily invertible.
% 
% Created by Dirk Poot, Erasmus MC
% 6-2-2013

Mx = R \(R' \ x);