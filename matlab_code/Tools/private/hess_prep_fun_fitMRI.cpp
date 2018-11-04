/*
 * This is a MEX-file for MATLAB.
 * compile: 
 *   mex hess_prep_fun_fitMRI.cpp
 * 
 * [hesscoeff] = hess_prep_fun_fitMRI( hess_pdf, d_fun , funHessian_I, funHessian_J);
 *
 * Evaluates the same as:
 *  size(hess_pdf ) = [numim, numtr ] 
 *  hessLng = hess_pdf(:,:,ones(1,nI)) .* (d_fun(:,:,funHessian_I).*d_fun(:,:,funHessian_J)); 
 *  hesscoeff = reshape( -sum( hessLng ,1) , numtr, nI)'; 
 * 
 *
 * Created by Dirk Poot, Erasmus MC
 */
#include <iterator>
#include "mex.h"

using std::iterator_traits;

template<typename T> void hess_prep_mul( T hesscoeff_d , T hess_pdf, T d_fun_p, T funHessian_I_p, T funHessian_J_p, int numim, int numtr, int numvar, int nhessel )
{

    typedef typename iterator_traits< T >::value_type value_type;
	for (int hid = 0 ; hid< nhessel; ++hid) {
        int I = funHessian_I_p[ hid ]-1;
        int J = funHessian_J_p[ hid ]-1;
        if ((I<0) | (J<0) | (I>=numvar) | (J>=numvar)) 
            mexErrMsgTxt("Invalid value in funHessian_I or funHessian_J.");
        
        T hess_pdf_l = hess_pdf;
        T d_fun_p_I = d_fun_p + I * numim * numtr ;
        T d_fun_p_J = d_fun_p + J * numim * numtr ;
        
        for (int tr = 0 ; tr<numtr; ++tr) {
            value_type cursum = 0;
            for (int imnr = 0 ; imnr < numim ; ++imnr) {
                cursum += (*hess_pdf_l) * (*d_fun_p_I) * (*d_fun_p_J);
                ++hess_pdf_l;
                ++d_fun_p_I;
                ++d_fun_p_J;
            }
            *hesscoeff_d = -cursum;
            ++hesscoeff_d;
        }
    }
};


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  /* Check for proper number of arguments. */
  if(nrhs!=4) {
    mexErrMsgTxt("Four inputs required. [hesscoeff] = hess_prep_fun_fitMRI( hess_pdf, d_fun , funHessian_I, funHessian_J); ");
  } else if(nlhs!=1) {
    mexErrMsgTxt("One output required.");
  }
  
  const mxArray *hess_pdf     = prhs[0];
  const mxArray *d_fun        = prhs[1];
  const mxArray *funHessian_I = prhs[2];
  const mxArray *funHessian_J = prhs[3];
          
  /* The input must be a noncomplex double.*/
  if( !mxIsDouble(hess_pdf) || mxIsComplex(hess_pdf) )
      mexErrMsgTxt("hess_pdf should be non complex double matrix");
  double * hess_pdf_p  = mxGetPr(hess_pdf);
  int numim  = mxGetM(hess_pdf);
  int numtr  = mxGetN(hess_pdf);
   
  if( !mxIsDouble(d_fun) || mxIsComplex(d_fun) )
      mexErrMsgTxt("d_fun should be non complex double matrix");
  const int * d_fun_sz = mxGetDimensions(d_fun);
  if ( (numim!=d_fun_sz[0]) || (numtr!=d_fun_sz[1]) ) {
      mexErrMsgTxt("size of d_fun should match size of hess_pdf in fist 2 dimensions.");
  }
  int numvar = 1;
  if (mxGetNumberOfDimensions(d_fun)>=3) {
      numvar = d_fun_sz[2];
  };
  double * d_fun_p = mxGetPr(d_fun);
  
  int nhessel  =  mxGetNumberOfElements( funHessian_I);
  if( !mxIsDouble(funHessian_I) || mxIsComplex(funHessian_I) || 
      !mxIsDouble(funHessian_J) || mxIsComplex(funHessian_J) ||
      nhessel != mxGetNumberOfElements( funHessian_J)    )
      mexErrMsgTxt("funHessian_I and funHessian_J should be non complex double vectors of equal length.");
  double * funHessian_I_p =  mxGetPr(funHessian_I);
  double * funHessian_J_p =  mxGetPr(funHessian_J);
          
  // Create output:
  plhs[0] = mxCreateDoubleMatrix( numtr, nhessel , mxREAL);
  double * hesscoeff_p = mxGetPr(plhs[0]);
  
  /* Call the subroutine. */
  hess_prep_mul<double *>( hesscoeff_p , hess_pdf_p, d_fun_p, funHessian_I_p, funHessian_J_p, numim, numtr, numvar, nhessel );
}
