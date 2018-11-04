#ifndef COMPACTEDHESSIANMULTIPLY_C
#define COMPACTEDHESSIANMULTIPLY_C
/*
 * This compactedHessianMultiply_c.cpp is a MEX-file for MATLAB.
 * Created by Dirk Poot, Erasmus Medical Center 23-9-2013, 
 * - Build with: (add -g for debugging (add -O for optimizations while debugging); -v for verbose mode )
 *   mex compactedHessianMultiply_c.cpp -I".\Tools\supportingfiles\" 
 * Or on Linux:
 *   mex compactedHessianMultiply_c.cpp -g CXXFLAGS="\$CXXFLAGS -mssse3 -mtune=core2"
 *
 * This function computes:
 * HY = compactedHessianMultiply_c( hessinfo, Y )
 *
 * rdindx = 1;
 * for k1 = 1:size(Y,1)
 *   for k2 = 1:k1-1
 *       HY(k1,:) = HY(k1,:) + hessinfo(rdindx,:).*Y(k2,:);
 *       HY(k2,:) = HY(k2,:) + hessinfo(rdindx,:).*Y(k1,:);
 *       rdindx = rdindx+1;
 *   end;
 *   HY(k1,:) = HY(k1,:) + hessinfo(rdindx,:).*Y(k1,:);
 *   rdindx = rdindx+1;
 * end;
 */

#include "mex.h"
#include "step_iterators.cpp"
#include "emm_vec.cpp"

template < typename HYtype, typename Htype, typename Ytype >
 inline void evaluateCompactedHessianMul1( HYtype HY, Htype H, const Ytype Y, int sizey )  {
    // core routine doing the computations
    typedef typename std::iterator_traits< HYtype >::value_type T;
    for (int k1 = 0 ; k1< sizey; ++k1) {
        T HYk1( 0.0 );
        T Yk1  = Y[k1];
        for (int k2 = 0; k2<k1; ++k2) {
            HYk1   += ( (*H) * (Y[k2]) );
            HY[k2] += ( (*H) * (Yk1)  );
            ++H;
        }
        HY[k1] = HYk1 + *H * Yk1  ;
        ++H;
    }
}

template < typename HYtype, typename Htype, typename Ytype >
        void evaluateCompactedHessianMulN( HYtype HY, Htype H, Ytype Y, int sizey , int Nits)  {
    // iterates over core routine to perform multiple hessian multiplications
    for(; Nits>0; --Nits) {
        evaluateCompactedHessianMul1< typename HYtype::value_type, typename Htype::value_type, typename Ytype::value_type>( (*HY), (*H), (*Y), sizey );
        ++HY;++H;++Y;
    }
}
    
template < typename datatype > 
    void evaluateCompactedHessianMul( datatype * HY, const datatype * H , const datatype * Y, int sizey, int sizeh, int Nits) {
    // routine that from pointers constructs the iterator types for the internal routines
    const int vlen = 4;
    if (Nits>=vlen) {
        // vectorized case
        typedef vecptr_step< const datatype * , vlen > readType_v;
        typedef vecptr_step<       datatype * , vlen > writeType_v;
        typedef step_iterator2D< readType_v > readType_v2;
        typedef step_iterator2D< writeType_v > writeType_v2;
        int Nitsv = Nits / vlen;
        writeType_v2 HYt( writeType_v( HY , sizey ), sizey * vlen );
        readType_v2  Yt(  readType_v(  Y  , sizey ), sizey * vlen  );
        readType_v2  Ht(  readType_v(  H  , sizeh ), sizeh * vlen  );
        evaluateCompactedHessianMulN( HYt, Ht, Yt, sizey, Nitsv );
        Nits -= Nitsv * vlen;
        HY += sizey * Nitsv * vlen;
        H  += sizeh * Nitsv * vlen;
        Y  += sizey * Nitsv * vlen;
    };
    if (Nits>0) {
        // end case.
        typedef step_iterator2D< const datatype * > readType_v2;
        typedef step_iterator2D< datatype *> writeType_v2;
        
        writeType_v2 HYt( HY , sizey );
        readType_v2 Ht( H , sizeh );
        readType_v2 Yt( Y , sizey );
        
        evaluateCompactedHessianMulN( HYt, Ht, Yt, sizey, Nits );
    }
}
        
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
   
	/* Check for proper number of arguments. */

	if ((nrhs!=2)) {
		mexErrMsgTxt("Two inputs required: [out ] = compactedHessianMultiply_c( H, Y ).");
	} else if (nlhs!=1) {
		mexErrMsgTxt("One output required.");
	}
	const mxArray * H	= prhs[0]; // N-D array, single or double.
	const mxArray * Y	= prhs[1]; // N-D array, single or double.
	
	// parse inputs:
	size_t sizeh = mxGetM( H );
	size_t sizey = mxGetM( Y );
    size_t Nits  = mxGetN( H );
    if (  (mxGetClassID(H)!=mxDOUBLE_CLASS) || (mxGetClassID(Y)!=mxDOUBLE_CLASS) 
        || mxIsComplex(H) || mxIsComplex(Y) 
        || (Nits!=mxGetN( Y) ) ) {
        mexErrMsgTxt("Inputs should be non complex double with same number of columns.");
    }
    if (sizey*(sizey+1)/2 != sizeh) {
        mexErrMsgTxt("Invalid number of rows in H or Y.");
    }
	const double * Hp = mxGetPr(H);
	const double * Yp = mxGetPr(Y);
    
    // Create output
    int ndimsIn = mxGetNumberOfDimensions( Y );
    const mwSize * sizeY = mxGetDimensions( Y );

    plhs[0] = mxCreateNumericArray((int) ndimsIn, sizeY , mxDOUBLE_CLASS, mxREAL);
    double* HYp = mxGetPr(plhs[0]);

    // compute:
    evaluateCompactedHessianMul( HYp, Hp , Yp, sizey, sizeh, Nits);
}
#endif