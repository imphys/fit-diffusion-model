/*
 * This TotalVariation_Vector_fc is a MEX-file for MATLAB.
 * Created by Dirk Poot, Erasmus MC
 * Build with: (add -g for debugging (add -o for optimizations while debugging); -v for verbose mode )
 * mex TotalVariation_Vector_fc.cpp -largeArrayDims

 [ TV, grad, hessinfo ] = TotalVariation_Vector_fc( x, spacing, weights, offset )
	% see totalVariationVecRegularizer.m for explaination

	TV = sum_i sqrt(  offset^2 + num_n  (x_i-x_n)'* weights * (x_i-x_n) ) - offset.
	% where n = all forward neigbors of element i.

 */
#define ONLYFORWARDGRADIENT
#include "TotalVariation_Vector_c.cpp"