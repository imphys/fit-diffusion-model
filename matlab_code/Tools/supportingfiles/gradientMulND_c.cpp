/*
 * This gradientMulND_c is a MEX-file for MATLAB.
 * Created by Dirk Poot, Erasmus MC
 * Build with: (add -g for debugging (add -o for optimizations while debugging); -v for verbose mode )
 * mex gradientMulND_c.cpp

 */
#include "mex.h"
//#include "emm_wrap.cpp"
//#include <complex>
//#include <omp.h>

#ifdef ANISOTROPICVOXELS
    #define VOXELSPACINGARGS , const imgtype * voxelspacing
    #define VOXELSPACING(dim) voxelspacing[ dim ]
#else
    #define VOXELSPACINGARGS 
    #define VOXELSPACING(dim) 1
#endif
template <typename imgtype, int ndims >
inline void gradientfilter_core( const imgtype * __restrict source, imgtype * __restrict dest, 
						   const mwSize * __restrict size, const mwSize * __restrict stepSrc, const mwSize * __restrict stepDst,
                           imgtype lambda  VOXELSPACINGARGS )
/* Evaluate gradient over valid area. 
	Requires:
	size[i]>=2
	(otherwise, you can remove that dim; it will not contribute to gradient)
 */
{
    int pos[ndims];
    // initialize pos:
	imgtype lambdasc = lambda * ((imgtype) (-.25)); // scale all contributions by -4 to undo this scaling of lambda.
    imgtype centersc = 0;	// .5 *ndim * (-4)
    for (int dim=0;dim<ndims;dim++) {
        pos[dim]=0;
		centersc += -2* VOXELSPACING( dim );
    }
    // compute gradient:
    for (;pos[ndims-1]<size[ndims-1];) {
        // Compute pointers to start of line:
        const imgtype * __restrict sourcep = source;
        imgtype * __restrict destp = dest;
		bool edgecase=false;
        for (int dim=0;dim<ndims;dim++) {
            sourcep += stepSrc[dim]*pos[dim];
            destp += stepDst[dim]*pos[dim];
			edgecase |= ((dim>0) && ((pos[dim]<2) || (pos[dim]>size[dim]-3)));
        }
        // evaluate gradient of line:
		if (!edgecase) {
			// no edge cases (Except for pos0)
			if (size[0]==2) {
				// special case if size[0]==2, gradient differs.
				// G = [ 2 -2
				//	    -2  2 ]
				// pos0 == 0:
				{imgtype reslt = (centersc - 6 * VOXELSPACING( 0 )) * sourcep[0]; // -6 = 2*-4 - (-2) 
				reslt += 8* VOXELSPACING( 0 ) * sourcep[stepSrc[0]]; // 8=-2 *-4
				for (int dim=1;dim<ndims;dim++) {
					reslt += (sourcep[-2*stepSrc[dim]] + sourcep[2*stepSrc[dim]]) * VOXELSPACING( dim ) ;
				}
				destp[0] = reslt * lambdasc;
				sourcep += stepSrc[0]; destp += stepDst[0];}

				// pos0 == 1:
				{imgtype reslt = (centersc - 6 * VOXELSPACING( 0 )) * sourcep[0]; // -6 = 2*-4 - (-2) 
				reslt += 8* VOXELSPACING( 0 ) * sourcep[-stepSrc[0]]; // 8=-2 *-4
				for (int dim=1;dim<ndims;dim++) {
					reslt += (sourcep[-2*stepSrc[dim]] + sourcep[2*stepSrc[dim]]) * VOXELSPACING( dim ) ;
				}
				destp[0] = reslt * lambdasc;
				sourcep += stepSrc[0]; destp += stepDst[0];}
			} else {
				// size[0]>=3
				// first and last element is normal:
				
				// pos0 == 0:
				// G = [1.25 -1 -.25]
				{imgtype reslt = (centersc - 3 * VOXELSPACING( 0 )) * sourcep[0]; // -3 = 1.25*-4 - (-2)
				reslt += 4 *VOXELSPACING( 0 ) * sourcep[stepSrc[0]];	  // 4=-1 *-4
				reslt +=    VOXELSPACING( 0 ) * sourcep[2*stepSrc[0]];	  // 1=-.25 *-4
				for (int dim=1;dim<ndims;dim++) {
					reslt += (sourcep[-2*stepSrc[dim]] + sourcep[2*stepSrc[dim]]) * VOXELSPACING( dim ) ;
				}
				destp[0] = reslt * lambdasc;
				sourcep += stepSrc[0]; destp += stepDst[0];}

				// pos0 == 1:
				// G = [-1 2 1] or [-1 1.25 0 -.25]
				{imgtype reslt;
				if (size[0]==3) {
					// second point special if size[0] ==3
					reslt = (centersc - 6 * VOXELSPACING( 0 )) * sourcep[0]; // -6 = 2*-4 - (-2)
					reslt += 4 * sourcep[-stepSrc[0]];	  // 4=-1 *-4
					reslt += 4 * sourcep[ stepSrc[0]];	  // 4=-1 *-4
				} else {
					// normal second point:
					reslt = (centersc - 3 * VOXELSPACING( 0 )) * sourcep[0]; // -3 = 1.25*-4 - (-2)
					reslt += 4 * VOXELSPACING( 0 )* sourcep[  -stepSrc[0] ];	  // 4=-1 *-4
					reslt +=     VOXELSPACING( 0 )* sourcep[ 2*stepSrc[0] ];	  // 1=-.25 *-4
				}
				for (int dim=1;dim<ndims;dim++) {
					reslt += (sourcep[-2*stepSrc[dim]] + sourcep[2*stepSrc[dim]]) * VOXELSPACING( dim ) ;
				}
				destp[0] = reslt * lambdasc;
				sourcep += stepSrc[0]; destp += stepDst[0];}
			}

			// main part, pos0= 2 to size[0]-3:
			//  convolve with [-.25 0 .5 0 -.25]
			for (int pos0=2;pos0<size[0]-2;pos0++) {
				imgtype reslt = centersc * sourcep[0];
				for (int dim=0;dim<ndims;dim++) {
					reslt += (sourcep[-2*stepSrc[dim]] + sourcep[2*stepSrc[dim]]) * VOXELSPACING( dim ) ;
				}
				destp[0] = reslt * lambdasc;
				sourcep += stepSrc[0]; destp += stepDst[0];
			}

			// do point end-1, if size[0]>3
			// G = [-.25 0 1.25 -1]
			if (size[0]>3) {
				imgtype reslt= (centersc - 3 * VOXELSPACING( 0 )) * sourcep[0]; // -5 = 1.25*-4
				reslt +=     VOXELSPACING( 0 ) * sourcep[-2*stepSrc[0]];	  // 1=-.25 *-4
				reslt += 4 * VOXELSPACING( 0 ) * sourcep[   stepSrc[0]];	  // 4=-1 *-4
				for (int dim=1;dim<ndims;dim++) {
					reslt += (sourcep[-2*stepSrc[dim]] + sourcep[2*stepSrc[dim]]) * VOXELSPACING( dim ) ;
				}
				destp[0] = reslt * lambdasc;
				sourcep += stepSrc[0]; destp += stepDst[0];
			}
			// last point, if size[0]>2
			// G = [-.25 -1 1.25]
			if (size[0]>2) {
				imgtype reslt= (centersc - 3 * VOXELSPACING( 0 )) * sourcep[0]; // -5 = 1.25*-4
				reslt += 4 * VOXELSPACING( 0 ) * sourcep[-1*stepSrc[0]];	  // 4=-1 *-4
				reslt +=     VOXELSPACING( 0 ) * sourcep[-2*stepSrc[0]];	  // 1=-.25 *-4
				for (int dim=1;dim<ndims;dim++) {
					reslt += (sourcep[-2*stepSrc[dim]] + sourcep[2*stepSrc[dim]]) * VOXELSPACING( dim ) ;
				}
				destp[0] = reslt * lambdasc;
			}
		} else {
			// one or more edge cases.
			for (int pos0=0;pos0<size[0];pos0++) {
				pos[0] = pos0;
				imgtype reslt = 0;
				for (int dim=0;dim<ndims;dim++) {
					if (pos[dim]==0) {
						if (size[dim]==2) {
							reslt += (- 8* sourcep[0]+ 8*sourcep[1*stepSrc[dim]]) * VOXELSPACING( dim ) ;
						} else {
							reslt += (sourcep[2*stepSrc[dim]] - 5* sourcep[0]+ 4*sourcep[1*stepSrc[dim]]) * VOXELSPACING( dim ) ;
						}
					} else if (pos[dim]==1) {
						if (size[dim]==2) {
							reslt += (- 8* sourcep[0]+ 8*sourcep[-1*stepSrc[dim]]) * VOXELSPACING( dim ) ;
						} else if (size[dim]==3) {
							reslt += (4*sourcep[-1*stepSrc[dim]] - 8* sourcep[0]+ 4*sourcep[1*stepSrc[dim]]) * VOXELSPACING( dim ) ;
						} else {
							reslt += (4*sourcep[-1*stepSrc[dim]] - 5* sourcep[0]+ 1*sourcep[2*stepSrc[dim]]) * VOXELSPACING( dim ) ;
						}
					} else if (pos[dim]==size[dim]-2) {
						reslt += (4*sourcep[+1*stepSrc[dim]] - 5* sourcep[0]+ 1*sourcep[-2*stepSrc[dim]]) * VOXELSPACING( dim ) ;
					} else if (pos[dim]==size[dim]-1) {
						reslt += (sourcep[-2*stepSrc[dim]] - 5* sourcep[0]+ 4*sourcep[-1*stepSrc[dim]]) * VOXELSPACING( dim ) ;
					} else {
						reslt += (sourcep[-2*stepSrc[dim]] - 2* sourcep[0]+ sourcep[2*stepSrc[dim]]) * VOXELSPACING( dim ) ;
					}
				}
				destp[0] = reslt * lambdasc;
				sourcep += stepSrc[0]; destp += stepDst[0];
			}
		}

        // move to next column:
		if (ndims>1) {
			int dim =0;
			do {
				pos[dim]=0;
				dim++;
				pos[dim]++;
			} while ((pos[dim]>=size[dim]) && (dim<ndims-1));
		} else { // do not loop if not more than 1 dimension:
			break;
		}
    }
}

#undef VOXELSPACINGARGS 
#undef VOXELSPACING

#ifndef ANISOTROPICVOXELS    
#define ANISOTROPICVOXELS
#include __FILE__
#undef ANISOTROPICVOXELS

template <typename imgtype, int ndims>
void gradientfilter_d1( const imgtype * __restrict source,   imgtype * __restrict dest, 
						   const mwSize * __restrict size, const mwSize * __restrict stepSrc, const mwSize * __restrict stepDst,
                           imgtype lambda , const imgtype * voxelspacing)
{
    if (voxelspacing==NULL) {
        gradientfilter_core<imgtype, ndims>( source,  dest, size, stepSrc, stepDst, lambda);
    } else {
        gradientfilter_core<imgtype, ndims>( source,  dest, size, stepSrc, stepDst, lambda, voxelspacing);
    }
}
template <typename imgtype>
void gradientfilter( const imgtype * __restrict source,  imgtype * __restrict dest, 
					const mwSize * __restrict size , const mwSize * __restrict stepSrc, const mwSize * __restrict stepDst,
                           imgtype lambda, const imgtype * voxelspacing, int ndims)
{
    if (ndims==4) {
        gradientfilter_d1<imgtype, 4>( source,  dest, size, stepSrc, stepDst, lambda, voxelspacing);
    } else if (ndims==3) {
        gradientfilter_d1<imgtype, 3>( source,  dest, size, stepSrc, stepDst, lambda, voxelspacing);
    } else if (ndims==2) {
        gradientfilter_d1<imgtype, 2>( source,  dest, size, stepSrc, stepDst, lambda, voxelspacing);
    } else if (ndims==1) {
        gradientfilter_d1<imgtype, 1>( source,  dest, size, stepSrc, stepDst, lambda, voxelspacing);
    } else 
        mexErrMsgTxt("unsuported #dims.");
}
/*
struct tempspace {
	char * base;
	mwSize
}
*/
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
/*  out = gradientMulND_c(src, lambda)
*/
 
	// * Check for proper number of arguments. *
	if ((nrhs!=2)&&(nrhs!=3)) {
		mexErrMsgTxt("Two or inputs required: out = gradientMulND_c(src, lambda [, reciprocalvoxelspacing] ).");
	} else if(nlhs!=1) {
		mexErrMsgTxt("One outputs required.");
	}

	// * parse inputs * /
	//in:
	int ndims = mxGetNumberOfDimensions(prhs[0]);

	if ((ndims>=5) || (mxGetNumberOfElements(prhs[1])!=1)) { 
		mexErrMsgTxt("Lambda not scalar or too many dimensions in input image, (increase maximum in the code of this function).");
	}
	mxClassID imgclass = mxGetClassID(prhs[0]);
    
	if ((imgclass!=mxDOUBLE_CLASS) & (imgclass!= mxSINGLE_CLASS) ) {
		mexErrMsgTxt("Image input should be single or double. (Other types can easily be added)");
	}
    if (mxIsComplex(prhs[0])) {
        mexErrMsgTxt("Image input should non complex");
    }
    
	const mwSize * sizein = mxGetDimensions(prhs[0]);
	if ((ndims==2) && (sizein[1]==1)) {
		ndims = 1;
	}
	double * inR = mxGetPr(prhs[0]);
    //mxComplexity = mxGetComplexity
    double lambda = mxGetScalar(prhs[1]);
    double * voxelspacing = NULL;
    if (nrhs>=3) {
        if ((imgclass!= mxGetClassID(prhs[2])) || (mxGetNumberOfElements(prhs[2])!=ndims) ) {
            mexErrMsgTxt("Reciprocal Voxel Spacing should be of same class as image and have ndims elements.");
        }
        voxelspacing = mxGetPr(prhs[2]);
    }
    
    mwSize cumsz = 1;
    for (int dim = 0; dim<ndims;dim++){
        cumsz *= sizein[dim];
		if (sizein[dim]<2) {
			mexErrMsgTxt("scalar dimensions currently not allowed, (use squeeze, or transpose row vectors).");
		}
    }
    plhs[0] = mxCreateNumericArray(ndims, sizein, imgclass, mxREAL);
    double * outR = mxGetPr(plhs[0]);
    
	mwSize n_tempspace_bytes = sizeof(mwSize) * ( (ndims) );
    char * tempspace = (char *) mxMalloc( n_tempspace_bytes );
    mwSize * stepDim = (mwSize *) tempspace;
    char * nexttempspace = tempspace + sizeof(mwSize) * ndims;
	mwSize next_tempspace_bytes = n_tempspace_bytes - sizeof(mwSize) * ndims;
    stepDim[0] = 1;
    for (int dim = 1; dim<ndims; dim++) {
        stepDim[dim] = stepDim[dim-1] * sizein[dim-1];
    }
    if (imgclass==mxDOUBLE_CLASS) {
		/* DEBUG:
		for (int k =0; k<cumsz;k++) {
			outR[k]=2;
		} // */
        gradientfilter<double>(         inR,           outR, sizein, stepDim, stepDim,         lambda,           voxelspacing, ndims);
    } else if (imgclass==mxSINGLE_CLASS) {
        gradientfilter<float>((float *) inR, (float *) outR, sizein, stepDim, stepDim, (float) lambda, (float *) voxelspacing, ndims);
    } else 
        mexErrMsgTxt("Unsupported image type.");

    mxFree(tempspace);
}
#endif

