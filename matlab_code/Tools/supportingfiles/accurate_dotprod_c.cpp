/* mex accurate_dotprod_c.cpp
 * optional arguments : -v (verbose) , -g (debug symbols)
 * NOTE: this function is only efficient on machines with (approximately) 11 or more XMM registers, so 64-bit x86; NOT on 32-bit x86.
 *
 * [C, Cl] = accurate_dotprod_c(A, B)
 * This routine computes C = A'*B
 * But with higher accuracy than standard (MATLAB) matrix multiplication.
 * It is slower though. 
 *
 * a demonstration of difference:
 *	A = ones(1000,4); A = ones(1000,4) + 1e-7*randn(size(A));
 *	B = ones(size(A)).*(1-2*(rand(size(A))>.5)) + 1e-7*randn(size(A));
 *	B(end,:) = -sum(B(1:end-1,:));
 *	Ci=A'*B,C=accurate_dotprod_c(A,B),
 *	(Ci-C)./C % relative error in Ci
 *
 * The relative error in Ci is quite large (1e-7). C is most likely even correctly rounded.
 */
#include "mex.h"
//#include  <stdlib.h>
#include <cmath>
#include <algorithm>
#include "emm_wrap.cpp"


void dotprod_double(double pa[], double pb[], double dsth[], double dstl[], int hAB, int stepcolAB, int wA, int wB)
{
	/* Compute the dotproduct of A and B accurately.
	 */
    const int vlen = 2;
    const double upscaleforround = 134217728.; // 2^27
    typedef vec<double, vlen> v;
	v upscaleforroundv = v(upscaleforround);
    double tempstore[2*vlen+(vlen-1)*3];
    int kb =0; // current column index in B
	for ( ; kb <wB; kb++) {
        int ka = 0; // current column index in A
		for (; ka <wA; ka++) {
            v accum_h = vec<double, vlen>(0.0);
            v accum_l = vec<double, vlen>(0.0);
			double * palp = pa + ka*stepcolAB;
			double * pblp = pb + kb*stepcolAB;
			int i=0;
			for ( ; i<hAB-1 ; i+=vlen) { // loop over 1 column of A and B

				// Note that to get high accuracy, we need to do the multiplication and adding a bit difficult.
				// We split A_i and B_i in 2 parts, a 'high' and 'low' part.
				// The high part is created such that it can be multiplied exactly. (i.e. without roundoff error)
				// Also the high-low parts can be multiplied exactly.
				// by adding and then subtracting again, the roundoff error due to the adding can be found.
				// 

				// Note that an fused multiply add instruction will be very usefull, when the intermediate result has infinite precision.
				// In this respect, the _mm_dp_pd intrinsic might have been usefull, but unfortunately also does round the intermediate result.

				// How many registers are needed?
				// upscaleforroundv, accum_(h/l) , cur(A/B)(_h), tmp (cur(A/B)_sc), curprodh, (curprodhl, curprodl, curprodh_backup can reuse curA/B(_h)), condition for the swap, accum_h_backup
				//       1         +   2         +      4      +    1             +   1     +                                                            +           1           +         1 
				//  = 11
                v curA = v(palp);
                v curB = v(pblp);
                v curA_sc = curA * upscaleforroundv;
                v curA_h = (curA + curA_sc) - curA_sc ; // round curA
                curA -= curA_h; // now only low part.

                v curB_sc = curB * upscaleforroundv;
                v curB_h = (curB + curB_sc) - curB_sc ; // round curA
                curB -= curB_h; // now only low part.
                 
                v curprodh = curA_h * curB_h;  // exact
                v curprodhl = curA * curB_h + curA_h * curB;// exact
                v curprodl = curA * curB;// exact
                v curprodh_backup = curprodh ;
                curprodh += curprodhl;
                // a + b + c = (a + b)  + (c + ((a-(a+b)) +b))
                curprodl += (( curprodh_backup - curprodh) + curprodhl); // add remainder of curprodhl to curprodl
                
                conditional_swap(accum_h, curprodh, abs(accum_h)<abs(curprodh)); 
                // abs(accum_h) >= abs(curprodh)
                v accum_h_backup = accum_h;
                accum_h = accum_h+curprodh; // inexact, but roundoff recovered in next line
                v remaccum = (accum_h_backup-accum_h) + curprodh; // exact; recover roundoff error in new accum_h computation
                accum_l += remaccum + curprodl; // inexact.

                palp +=vlen;
				pblp +=vlen;
			};
            accum_h.store(tempstore);
            accum_l.store(tempstore+vlen);
			int ntempstore = 2*vlen;
			for ( ; i<hAB ; i++) {
				double curA = palp[0];
                double curB = pblp[0];
                double curA_sc = curA * upscaleforround;
                double curA_h = (curA + curA_sc) - curA_sc ; // round curA
                curA -= curA_h; // now only low part.

                double curB_sc = curB * upscaleforround;
                double curB_h = (curB + curB_sc) - curB_sc ; // round curA
                curB -= curB_h; // now only low part.
                 
                tempstore[ntempstore] = curA_h * curB_h;  // exact
                tempstore[ntempstore+1] = curA * curB_h + curA_h * curB;// exact
                tempstore[ntempstore+2] = curA * curB;// exact

				ntempstore +=3;
                palp ++;
				pblp ++;
			}

            double accum_h_final = tempstore[0];
            double accum_l_final = 0;
            for (int k = 1; k<=ntempstore; k++) {
				double ldval;
				if (k==ntempstore) { // in last iteration, add low value to provide 'perfect' rounding.
					ldval = accum_l_final;
					accum_l_final =0;
				} else {
					ldval = tempstore[k];
				}
				using std::abs;
				if (abs(accum_h_final)<abs(ldval)) {
					double tmp = ldval;
					ldval = accum_h_final;
					accum_h_final = tmp;
				}
				// accum_h_final has largest magnitude.
                double oldaccum = accum_h_final;
                accum_h_final += ldval;
                accum_l_final += (oldaccum - accum_h_final) + ldval;
            };
			
            dsth[ ka + wA * kb ] = accum_h_final;
            if (dstl) {
				// if requested, also store low part. Note that this typically contains substantial roundoff errors.
                dstl[ ka + wA * kb ] = accum_l_final;
            }
		};

	};
    

};

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
/* out = point_dist_c(A, B)
*/
	double *A,*B;
	int lenA, lenB;

	/* Check for proper number of arguments. */
	if(nrhs!=2) {
		mexErrMsgTxt("Two inputs required: out = accurate_dotprod_c(A, B). ");
	} else if((nlhs!=1)&& (nlhs!=2)) {
		mexErrMsgTxt("One or two outputs required.");
	}

	/* parse inputs */
	/* // A: */ 
	if ((mxGetNumberOfDimensions(prhs[1])!=2) | (mxGetClassID(prhs[0])!=mxDOUBLE_CLASS) ) { 
		mexErrMsgTxt("incorrect Input A");
	}
	int wA = mxGetN(prhs[0]);
	lenA   = mxGetM(prhs[0]);
 	/*if ((lenA & 1) !=0) {
 		mexErrMsgTxt("lenA should be a multiple of 2.");
 	}*/
	A = mxGetPr(prhs[0]);

	lenB   = mxGetM(prhs[1]);
    int wB = mxGetN(prhs[1]);
	if ((mxGetNumberOfDimensions(prhs[1])!=2) | (mxGetClassID(prhs[1])!=mxDOUBLE_CLASS) | (lenB!=lenA) ) { 
		mexErrMsgTxt("incorrect Input B");
	}
//  	if ((lenB & 1) !=0) {
//  		mexErrMsgTxt("lenB should be a multiple of 2.");
//  	}
	B = mxGetPr(prhs[1]);

	plhs[0] = mxCreateDoubleMatrix(wA, wB, mxREAL);
	double * y = mxGetPr(plhs[0]);
	double * yl = NULL;
	if (nlhs>=2) {
		plhs[1] = mxCreateDoubleMatrix(wA, wB, mxREAL);
		yl = mxGetPr(plhs[1]);
	}
    dotprod_double(A, B, y, yl, lenA, lenA, wA, wB);
}

