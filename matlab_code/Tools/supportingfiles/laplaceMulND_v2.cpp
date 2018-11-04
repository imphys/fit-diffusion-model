/*
 * This laplaceMulND_v2 is a MEX-file for MATLAB.
 * Created by Dirk Poot, Erasmus MC
 * Build with: (add -g for debugging (add -o for optimizations while debugging); -v for verbose mode )
 * mex laplaceMulND_v2.cpp -I".\Tools\supportingfiles\"


 */

#define ADDGUARDS  // if ADDGUARDS is defined, tempSpaceManaging adds guard's around the temp variables, which can later be checked. 
				      // this file does check all variables, but only if ADDGUARDS is defined (to avoid potential overhead). 

#include "mex.h"
#include "tempSpaceManaging.cpp"
//#include "emm_wrap.cpp"
//#include <complex>
//#include <omp.h>
//#include <iterator>
//#define CHECKPOINTERBOUNDS
#include "guardpointer.cpp"


#ifdef ANISOTROPICVOXELS
    #define VOXELSPACINGARGS , const imgType * voxelspacing
    #define VOXELSPACING(dim) voxelspacing[ dim ]
#else
    #define VOXELSPACINGARGS 
	#define VOXELSPACING(dim) (dim==0 ? -nsteps : 1)  // condition should be optimized away.
#endif



template <int nsteps , typename sourceType, typename destType, typename imgType>
inline void laplacefilter_line( sourceType source, destType dest, 
								mwSize length, const mwSize * __restrict stepSrc, const mwSize stepDst,
								int boundaryMethod, 
						        imgType lambda  VOXELSPACINGARGS )
/* Evaluate laplacian of a line
 * boundaryMethod  :  0 = no boundary handling (in line dimension)
 *                    1 = gradient edges
 *                    2 = cyclic. 
 * for boundary cases:
 *  stepSrc[1] = forward step along line
 *  stepSrc[2] = backward step along line. 
 */
{
//	using std::iterator_traits;
//typedef typename iterator_traits< destType >::value_type        imgType;

    // initialize pos:
	imgType centersc = VOXELSPACING( 0 );
 
	int pos0 = 0;
	imgType firstlineval;
	if ( boundaryMethod>0 ) {
		--length;
		imgType reslt;
		if (boundaryMethod==2) {
			firstlineval = source[ 0 ];
			// cyclic domain:
			reslt = centersc * firstlineval + source[ stepSrc[1] ] * VOXELSPACING( 1 ) + source[ stepSrc[0] * length ] * VOXELSPACING( 2 ) ;
		} else {
			// gradient edges:
			reslt = (centersc+VOXELSPACING( 2 )) * source[0] + source[ stepSrc[1] ] * VOXELSPACING( 1 ) ;
			//reslt = centersc_e * source[0];
		}
		for (int step = 3; step<nsteps ; ++step) {
			reslt += source[ stepSrc[step] ] * VOXELSPACING( step ) ;
		}
		*dest = reslt*lambda;
		source += stepSrc[0]; 
		dest += stepDst;
		++pos0;
	} 

    for (; pos0 < length; ++pos0 ) {
        imgType reslt = centersc * source[0];
        for (int step = 1; step<nsteps ; ++step) {
            reslt += source[ stepSrc[step] ] * VOXELSPACING( step ) ;
        }
        *dest = reslt*lambda;
		source += stepSrc[0]; 
		dest += stepDst;
    }

	if ( boundaryMethod>0 ) {
		imgType reslt;
		if (boundaryMethod==2) {
			// cyclic domain:
			reslt = centersc * source[0] + firstlineval * VOXELSPACING( 1 ) ;
		} else {
			// gradient edges:
			reslt = (centersc+VOXELSPACING( 1 )) * source[0] ;
		}
		for (int step = 2; step<nsteps ; ++step) {
			reslt += source[ stepSrc[step] ] * VOXELSPACING( step ) ;
		}
		*dest = reslt*lambda;
//		source += stepSrc[0]; 
//		destp += stepDst;
//		++pos0;
	} 

}


#undef VOXELSPACINGARGS 
#undef VOXELSPACING

#ifndef ANISOTROPICVOXELS    
#define ANISOTROPICVOXELS
#include __FILE__
#undef ANISOTROPICVOXELS

template <typename sourceType, typename destType, typename imgType >
void laplacefilter_full( sourceType source,  destType  dest, 
						 const mwSize * __restrict size, const mwSize * __restrict stepSrc, const mwSize * __restrict stepDst,
						 const int ndims, const bool iscyclic,
                         imgType lambda , const imgType * voxelspacing , tempspace tmp)
{
	//using std::iterator_traits;
	// get type of what we get for a voxel:
	//typedef typename iterator_traits< Mtype >::value_type        imgType;

	mwSize * pos = tmp.get<mwSize>( ndims+1 );
	for (int dim = 0; dim < ndims; dim++) {
		pos[dim] = 0;
	};
	mwSize numloops = 1;
	for (int dim = 1; dim < ndims; dim++) {
		numloops *= size[dim];
	};
    mwSize * stepSrc_lp = tmp.get<mwSize>(1+2*ndims);
    imgType * voxelspacing_lp = tmp.get<imgType>(1+2*ndims);
	int boundaryMethod = (iscyclic ? 2 : 1 );
	int boundaryMethod_lp = boundaryMethod;
	if ( (voxelspacing!=NULL) && (voxelspacing[0]==0) ) { 
		boundaryMethod_lp = 0; // if voxelspacing is zero in first dimension, we don't step in that dimension, so also don't do any edge case handling. 
	}
	bool changedEdges = true;
	sourceType source_lp;
	destType   dest_lp;
	int nsteps;

	for ( ; numloops>0 ; --numloops ) {

		if (changedEdges ) {
			nsteps = 0;
			stepSrc_lp[ nsteps ] = stepSrc[ 0 ]; 
			source_lp = source;
			dest_lp =   dest;
			++nsteps;
			for ( int dim =0 ; dim < ndims ; ++dim ) { 
				source_lp += stepSrc[ dim ] * pos[dim];
				dest_lp   += stepDst[ dim ] * pos[dim];

#define Set_Step( stepsize ) \
				stepSrc_lp[ nsteps ] = stepsize;\
				voxelspacing_lp[ nsteps ] = (voxelspacing==NULL ? 1 : voxelspacing[dim]);\
				++nsteps

				if ((voxelspacing==NULL) || (voxelspacing[dim]!=0)) {
					if (pos[dim]<size[ dim ]-1) { // test if forward step can be made
						Set_Step( stepSrc[ dim ] );
					} else if (iscyclic) {
						//boundary wraps around: 
						Set_Step( -(size[ dim ]-1)* stepSrc[ dim ] );
					}
					if ((dim==0) || (pos[dim]>0)) { // test if backward step can be made
						Set_Step( -stepSrc[ dim ] );
					} else if (iscyclic) {
						//boundary
						Set_Step( (size[ dim ]-1) * stepSrc[ dim ] );
					}
				}
			};
			voxelspacing_lp[ 0 ] =0;
			for (int step =1 ; step < nsteps; ++step ) {
				voxelspacing_lp[ 0 ] -= voxelspacing_lp[ step ];
			};

		};

		if (nsteps<=6) {
		if (nsteps==2) {
			laplacefilter_line<2>( source_lp, dest_lp, size[0], stepSrc_lp, stepDst[0], boundaryMethod_lp, lambda,   voxelspacing_lp );
		} else if (nsteps==3) {
			laplacefilter_line<3>( source_lp, dest_lp, size[0], stepSrc_lp, stepDst[0], boundaryMethod_lp, lambda,   voxelspacing_lp );
		} else if (nsteps==4) {
			laplacefilter_line<4>( source_lp, dest_lp, size[0], stepSrc_lp, stepDst[0], boundaryMethod_lp, lambda,   voxelspacing_lp );
		} else if (nsteps==5) {
			laplacefilter_line<5>( source_lp, dest_lp, size[0], stepSrc_lp, stepDst[0], boundaryMethod_lp, lambda,   voxelspacing_lp );
		} else if (nsteps==6) {
			laplacefilter_line<6>( source_lp, dest_lp, size[0], stepSrc_lp, stepDst[0], boundaryMethod_lp, lambda,   voxelspacing_lp );
		} else 
			mexErrMsgTxt("Unsupported image dimension");
		} else if (nsteps==7) {
			laplacefilter_line<7>( source_lp, dest_lp, size[0], stepSrc_lp, stepDst[0], boundaryMethod_lp, lambda,   voxelspacing_lp );
		} else if (nsteps==8) {
			laplacefilter_line<8>( source_lp, dest_lp, size[0], stepSrc_lp, stepDst[0], boundaryMethod_lp, lambda,   voxelspacing_lp );
		} else if (nsteps==9) {
			laplacefilter_line<9>( source_lp, dest_lp, size[0], stepSrc_lp, stepDst[0], boundaryMethod_lp, lambda,   voxelspacing_lp );
		} else {
			mexErrMsgTxt("Unsupported image dimension");
		}

		// update position:
		{
		int dim =1;
		if (pos[dim]<1 || pos[dim]>=size[dim]-2) {
			changedEdges =true;
			dim=0;
			do {
				pos[dim]=0;
				dim++;
				pos[dim]++;
			} while ((pos[dim]>=size[dim]) && (dim<ndims-1));
		} else {
			// in dim 1, pos = 1... size[dim]-2 have the same edges
			// so if current and next pos fall in this range, use shortcut:
			changedEdges =false;
			++pos[dim];
			source_lp += stepSrc[ dim ];
			dest_lp   += stepDst[ dim ];
		}
		};
	}

};


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
/*  out = laplaceMulND_c(src, lambda)
*/
 
	// * Check for proper number of arguments. *
	if ((nrhs!=2)&&(nrhs!=3)&&(nrhs!=4)) {
		mexErrMsgTxt("Two or inputs required: out = laplaceMulND_c(src, lambda [, reciprocalvoxelspacing, iscyclic ] ).");
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
    bool iscyclic = false;
    if (nrhs>=4) {
        if (mxGetNumberOfElements(prhs[3])!=1) {
            mexErrMsgTxt("iscyclic should be a scalar (boolean).");
        }
        iscyclic = (mxGetScalar(prhs[3])!=0);
    }
    plhs[0] = mxCreateNumericArray(ndims, sizein, imgclass, mxREAL);
    double * outR = mxGetPr(plhs[0]);
    
	// Compute required temporary space, and allocate it:
    mwSize ntempdata = ndims;			 // tempiom in laplacefilter_add
	mwSize ntempsize = (2 * (1+2*ndims) + 3*ndims+1);   // stepTmp in laplacefilter_add
												 // nextsize, nextsize in laplacefilter_full
												 // stepDim in main
	mwSize n_tempspace_bytes = mxGetElementSize(prhs[0]) * ntempdata + sizeof(mwSize) * ntempsize + sizeof(double) * ndims;
#ifdef ADDGUARDS // check if variables are not corrupted:
	n_tempspace_bytes += 100; // add some bytes to have enough space for the guards.
#endif
    char * tempspaceptr = (char *) mxMalloc( n_tempspace_bytes );
	tempspace tmp = tempspace(tempspaceptr, n_tempspace_bytes);

	mwSize * size_sq = tmp.get<mwSize>( ndims ); 
	int ndims_nonscalar =0;
	for (int dim = 0; dim < ndims; dim++){
		if ( sizein[dim]!=1 ) {
			size_sq[ ndims_nonscalar ] = sizein[dim];
			++ndims_nonscalar;
		}
    }
    mwSize * stepDim = tmp.get<mwSize>( ndims_nonscalar ); 

	stepDim[0] = 1;
    for (int dim = 1; dim<ndims_nonscalar; dim++) {
        stepDim[dim] = stepDim[dim-1] * size_sq[dim-1];
    }
    if (imgclass==mxDOUBLE_CLASS) {
		typedef double baseT;
		typedef const baseT *  __restrict sourceType ;
		typedef       baseT *  __restrict destType ;

		mxGetPtr( sourceType, source_ptr,   prhs[0] );
		mxGetPtr( destType,   dest_ptr,     plhs[0] );
		
		destType voxelspacing_sq = (destType) voxelspacing;
		if (voxelspacing!=NULL) {
			voxelspacing_sq = tmp.get< baseT>( ndims_nonscalar );
			int d = 0;
			for (int dim = 0 ; dim < ndims; ++dim){
				if ( sizein[dim]!=1 ) {
					voxelspacing_sq[d] = ((baseT *) voxelspacing)[dim];
					++d;
				}
			}
		}
		
        laplacefilter_full<source_ptr_type, dest_ptr_type, baseT>( source_ptr, dest_ptr, size_sq, stepDim, stepDim, ndims_nonscalar, iscyclic , lambda, voxelspacing_sq, tmp);
    } else if (imgclass==mxSINGLE_CLASS) {
		typedef float baseT;
		typedef const baseT *  __restrict sourceType ;
		typedef       baseT *  __restrict destType ;

		mxGetPtr( sourceType, source_ptr,   prhs[0] );
		mxGetPtr( destType,   dest_ptr,     plhs[0] );

		destType voxelspacing_sq = (destType) voxelspacing;
		if (voxelspacing!=NULL) {
			voxelspacing_sq = tmp.get< baseT>( ndims_nonscalar );
			int d = 0;
			for (int dim = 0 ; dim < ndims; ++dim){
				if ( sizein[dim]!=1 ) {
					voxelspacing_sq[d] = ((baseT *) voxelspacing)[dim];
					++d;
				}
			}
		}

		laplacefilter_full<source_ptr_type, dest_ptr_type, baseT>( source_ptr, dest_ptr, size_sq, stepDim, stepDim, ndims_nonscalar, iscyclic , lambda, voxelspacing_sq, tmp);
    } else 
        mexErrMsgTxt("Unsupported image type.");

#ifdef ADDGUARDS // check if variables are not corrupted:
	tmp.checkGuard( stepDim  , ndims_nonscalar ); 
#endif    
	mxFree(tempspaceptr);
}
#endif

