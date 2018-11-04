/*
 * This laplaceMulND_c is a MEX-file for MATLAB.
 * Created by Dirk Poot, Erasmus MC
 * Build with: (add -g for debugging (add -o for optimizations while debugging); -v for verbose mode )
 * mex laplaceMulND_c.cpp

 */

#define ADDGUARDS  // if ADDGUARDS is defined, tempSpaceManaging adds guard's around the temp variables, which can later be checked. 
				      // this file does check all variables, but only if ADDGUARDS is defined (to avoid potential overhead). 

#include "mex.h"
#include "tempSpaceManaging.cpp"
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
inline void laplacefilter_core_p1( const imgtype * __restrict source, imgtype * __restrict dest, 
						   const mwSize * __restrict size, const mwSize * __restrict stepSrc, const mwSize * __restrict stepDst,
                           imgtype lambda  VOXELSPACINGARGS )
/* Evaluate laplacian over valid area. 
 * Source is read for (almost) all sum_(dim==0..ndims-1) stepSrc[ pos_in[dim] ]
 *  with 0 <= pos_in[  dim ] < size[ dim ]
 *
 * Dest is written to for all sum_(dim==0..ndims-1) stepDest[ pos_out[dim] ]
 *  with 1 <= pos_out[ dim ] < size[ dim ] - 1
 */
{
    int pos[ndims]; // pos_out
    // initialize pos:
    imgtype centersc = 0;
    for (int dim=0;dim<ndims;dim++) {
        pos[dim]=1;
        centersc += -2 * VOXELSPACING( dim );
		if (pos[dim]>=size[ndims-1]-1) {
			return;
		}
    }
    // compute laplacian first round:
    for (;pos[ndims-1]<size[ndims-1]-1;) {
        // Compute pointers to start of line:
        const imgtype * __restrict sourcep = source;
        imgtype * __restrict destp = dest;
        for (int dim=0;dim<ndims;dim++) {
            sourcep += stepSrc[dim]*pos[dim];
            destp += stepDst[dim]*pos[dim];
        }
        // evaluate laplacian of line:
        for (int pos0=1;pos0<size[0]-1;pos0++) {
            imgtype reslt = centersc * sourcep[0];
            for (int dim=0;dim<ndims;dim++) {
                reslt += (sourcep[-stepSrc[dim]] + sourcep[stepSrc[dim]]) * VOXELSPACING( dim ) ;
            }
            destp[0] = reslt*lambda;
			sourcep += stepSrc[0]; destp += stepDst[0];
        }
        // update position:
		if (ndims>1) {
			int dim =0;
			do {
				pos[dim]=1;
				dim++;
				pos[dim]++;
			} while ((pos[dim]>=size[dim]-1) && (dim<ndims-1));
		} else { // do not loop if not more than 1 dimension:
			break;
		}
    }
}

#ifndef ANISOTROPICVOXELS
template<typename imgtype> inline void write(  imgtype * storeto, const imgtype value) {
	storeto[0] = value;
};
template<typename imgtype> inline void accum(  imgtype * storeto, const imgtype value) {
// /*
	storeto[0] += value;
/*/
	imgtype total = storeto[0] + value;
	if ( !( (total>-1e20) && (total<1e20) ) ) {
		total = 123456;
	}
	storeto[0] = total;
//*/
};
#endif

template <typename imgtype, int ndims,  void store(imgtype*, imgtype) >
inline void laplacefilter_core_p2( const imgtype * source,  imgtype * dest, 
								   const mwSize * __restrict size, const mwSize * __restrict stepSrc, const mwSize * __restrict stepDst  VOXELSPACINGARGS )
/* Evaluate laplacian (adjoint) over valid area. 
 * Source is read for all sum_(dim==0..ndims-1) stepSrc[ pos_in[dim] ]
 *  with 1 <= pos_in[  dim ] < size[ dim ] -1
 *
 * Dest is written to for (almost) all sum_(dim==0..ndims-1) stepDest[ pos_out[dim] ]
 *  with 0 <= pos_out[ dim ] < size[ dim ] 
 */
{
    int pos[ndims]; // pos_out

    // initialize pos:
    imgtype centersc = 0;
    for (int dim=0;dim<ndims;dim++) {
        pos[dim]=0;
        centersc += -2 * VOXELSPACING( dim );
    }
    
    // compute laplacian second round:
    for (;pos[ndims-1]<size[ndims-1];) {
        // Compute pointers to start of line:
        const imgtype * sourcep = source;
        imgtype * destp = dest;
        int boundarycnt = 0;
        int boundarydim = 0;
        for (int dim=1;dim<ndims;dim++) { // start at dim==1 since pos[0] ==0 & that should not activate isboundary
            sourcep += stepSrc[dim]*pos[dim];
            destp += stepDst[dim]*pos[dim];
            if ((pos[dim]==0) || (pos[dim]>=size[dim]-1)) {
                boundarycnt += ndims;
                boundarydim = dim;
            } else if ((pos[dim]==1) || (pos[dim]==size[dim]-2)) {
				boundarycnt ++;
			}
        }
        // evaluate laplacian of line: (pos[0] = 0 .. size[0]-1 )
        if (boundarycnt>0) {
            if (boundarycnt>=2*ndims) {
                // empty when in 2 or more boundaries:
                for (int pos0 = 0;  pos0 < size[0] ; ++pos0 , destp += stepDst[0] ) {
                    store( destp, 0);
                }
            } else if (boundarycnt>=ndims) {
                // first and last element of line are in 2 boundary zones: set to zero
                store(destp, 0);
                sourcep += stepSrc[0];
                destp += stepDst[0];
				// step to the only line of src that can be accessed from the current line of dest:
                if (pos[boundarydim]==0) {
                    sourcep += stepSrc[boundarydim];
                } else { // end of line.
                    sourcep -= stepSrc[boundarydim];
                }
                for (int pos0 = 1 ; pos0 < size[0]-1; ++pos0 ) {
                    store(destp, sourcep[0] * VOXELSPACING( boundarydim ));
					sourcep += stepSrc[0];
					destp += stepDst[0];
                }
                store(destp,0);
			} else { // not in outer boundary, but some of the extending points might fall outside:
                
				store(destp, sourcep[stepSrc[0]] * VOXELSPACING( 0) );
				sourcep += stepSrc[0];
				destp   += stepDst[0];
				for (int pos0 = 1 ; pos0 < size[0]-1 ; ++pos0) {
                    imgtype reslt =  centersc * sourcep[0];
					pos[0] = pos0;
					for (int dim=0;dim<ndims;dim++) {
                        if (pos[dim]>1)
                            reslt += sourcep[-stepSrc[dim]] * VOXELSPACING( dim );
						if (pos[dim]<(size[dim]-2))
                            reslt += sourcep[ stepSrc[dim]] * VOXELSPACING( dim );
					}
					store(destp, reslt);
					sourcep += stepSrc[0]; destp += stepDst[0];
				}
                store(destp, sourcep[-stepSrc[0]] * VOXELSPACING( 0 ) );
			}
        } else { // not boundary:
			// first element in line:
             store(destp, sourcep[stepSrc[0]] * VOXELSPACING( 0 ) );
             sourcep += stepSrc[0];
             destp += stepDst[0];
			// second element in line:
			 {	 imgtype reslt =  centersc * sourcep[0];
				 reslt += sourcep[stepSrc[0]] * VOXELSPACING( 0 );
                 for (int dim=1 ; dim<ndims ; ++dim ) {
                     reslt += ( sourcep[-stepSrc[dim]] + sourcep[stepSrc[dim]] ) * VOXELSPACING( dim );
                 }
                 store(destp,reslt);
				 sourcep += stepSrc[0]; destp += stepDst[0];
			 }
			 // main part:
             for (int pos0 = 2; pos0 < size[0]-2 ; ++pos0) {
                 imgtype reslt =  centersc * sourcep[0];
                 for (int dim = 0 ; dim < ndims ; ++dim ) {
                     reslt += ( sourcep[-stepSrc[dim]] + sourcep[stepSrc[dim]] ) * VOXELSPACING( dim );
                 }
                 store(destp, reslt);
				 sourcep += stepSrc[0]; destp += stepDst[0];
             }
			 // one-but-last element in line
			 if (size[0]>3) {	 
				 imgtype reslt =  centersc * sourcep[0];
				 reslt += sourcep[-stepSrc[0]] * VOXELSPACING( 0 );
                 for (int dim = 1 ; dim < ndims ; ++dim ) {
                     reslt += (sourcep[-stepSrc[dim]] + sourcep[stepSrc[dim]]) * VOXELSPACING( dim );
                 }
                 store(destp,reslt);
				 sourcep += stepSrc[0]; destp += stepDst[0];
			 };
			 // last element in line:
             store(destp, sourcep[-stepSrc[0]] * VOXELSPACING( 0 ));
        }
        // update position:
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

template <typename imgtype, int ndims, void store(imgtype *, imgtype)>
void laplacefilter_core( const imgtype * __restrict source,  imgtype * temp, imgtype *  dest, 
						   const mwSize * __restrict size, const mwSize * __restrict stepSrc, const mwSize * __restrict stepTmp, const mwSize * __restrict stepDst,
                           imgtype lambda , const imgtype * voxelspacing)
{
	// offset temp origin just in case it is mapped to dest.
	for (int dim = 0; dim<ndims; dim++) {
		temp += stepTmp[dim];
	}
    if (voxelspacing==NULL) {
        laplacefilter_core_p1<imgtype, ndims>( source,  temp,       size, stepSrc, stepTmp,         lambda);
        laplacefilter_core_p2<imgtype, ndims, store>(   temp, dest, size,          stepTmp, stepDst );
    } else {
        laplacefilter_core_p1<imgtype, ndims>( source,  temp,       size, stepSrc, stepTmp,         lambda, voxelspacing);
        laplacefilter_core_p2<imgtype, ndims, store>(   temp, dest, size,          stepTmp, stepDst,        voxelspacing);
    }
}
template <typename imgtype>
void laplacefilter( const imgtype * __restrict source,  imgtype * __restrict dest, 
					const mwSize * __restrict size , const mwSize * __restrict stepSrc, const mwSize * __restrict stepDst,
                           imgtype lambda, const imgtype * voxelspacing, int ndims)
{
    if (ndims==4) {
        laplacefilter_core<imgtype, 4, write<imgtype> >( source,  dest, dest, size, stepSrc, stepDst, stepDst, lambda, voxelspacing);
    } else if (ndims==3) {
        laplacefilter_core<imgtype, 3, write<imgtype> >( source,  dest, dest, size, stepSrc, stepDst, stepDst, lambda, voxelspacing);
    } else if (ndims==2) {
        laplacefilter_core<imgtype, 2, write<imgtype> >( source,  dest, dest, size, stepSrc, stepDst, stepDst, lambda, voxelspacing);
    } else if (ndims==1) {
        laplacefilter_core<imgtype, 1, write<imgtype> >( source,  dest, dest, size, stepSrc, stepDst, stepDst, lambda, voxelspacing);
    } else 
        mexErrMsgTxt("unsuported #dims.");
}

template <typename imgtype>
void laplacefilter_add( const imgtype * __restrict source,  imgtype * __restrict dest, 
						const  mwSize * __restrict size , const mwSize * __restrict stepSrc, const mwSize * __restrict stepDst,
                           imgtype lambda, const imgtype * voxelspacing, int ndims,
                           tempspace tmp)
{
    mwSize cumsz =1;
    for (int dim = 0;dim<ndims; dim++){
        cumsz *= size[dim];
    }
    mwSize * stepTmp = tmp.get<mwSize>( ndims );
    stepTmp[0]=1;
    for (int dim = 1 ; dim < ndims ; ++dim ) {
        stepTmp[dim]=stepTmp[dim-1] * size[dim-1];
    }
	imgtype * tempim  = tmp.get<imgtype>( cumsz );
    if (ndims==3) {
        laplacefilter_core<imgtype, 3, accum<imgtype> >( source, tempim, dest, size, stepSrc, stepTmp, stepDst, lambda, voxelspacing);
    } else if (ndims==2) {
        laplacefilter_core<imgtype, 2, accum<imgtype> >( source, tempim, dest, size, stepSrc, stepTmp, stepDst, lambda, voxelspacing);
    } else if (ndims==1) {
        laplacefilter_core<imgtype, 1, accum<imgtype> >( source, tempim, dest, size, stepSrc, stepTmp, stepDst, lambda, voxelspacing);
    } else 
        mexErrMsgTxt("unsuported #dims.");
#ifdef ADDGUARDS // check if variables are not corrupted:
	tmp.checkGuard( stepTmp , ndims ); 
	tmp.checkGuard( tempim , cumsz ); 
#endif
}

template <typename imgtype>
void laplacefilter_full(const imgtype * __restrict source,  imgtype * __restrict dest, 
                        const mwSize * __restrict size , const mwSize * __restrict stepDim,
                           imgtype lambda, const imgtype * voxelspacing , int ndims, bool doadd, 
                           tempspace tmp ) 
{
	bool dolapl = true;
	for (int dim = 0 ; dim < ndims ; ++dim) {
		dolapl = dolapl && (size[dim]>2);
	}
	if (dolapl) {
		if (doadd) {
			laplacefilter_add( source, dest, size , stepDim, stepDim, lambda, voxelspacing, ndims, tmp);
		} else { 
			laplacefilter(     source, dest, size , stepDim, stepDim, lambda, voxelspacing, ndims);
		}
	}
	if (ndims>1) {
        // handle edge planes. (they are not fully incorporated in laplacefilter)
        mwSize * nextsize = tmp.get<mwSize>( ndims-1 );
        mwSize * nextstep = tmp.get<mwSize>( ndims-1 );
        
        imgtype * nextvoxelspacing = NULL;
        if (voxelspacing!=NULL) {
            nextvoxelspacing = tmp.get<imgtype>(ndims-1);
            for (int dim =0; dim < ndims-1 ; ++dim) {
                nextvoxelspacing[dim] = voxelspacing[dim+1];
            }
        }
		// First fix pos[0] to 0 or size[0]-1; then dim 1, and so on.
        for (int dim =0; dim<ndims-1 ; ++dim ) {
            nextsize[dim] = size[dim+1];
            nextstep[dim] = stepDim[dim+1];
        }
        for (int dim =0; dim < ndims ; ++dim) {
            laplacefilter_full( source                           , dest                           , nextsize , nextstep, lambda, nextvoxelspacing, ndims-1, true, tmp);
            laplacefilter_full( source+(size[dim]-1)*stepDim[dim], dest+(size[dim]-1)*stepDim[dim], nextsize , nextstep, lambda, nextvoxelspacing, ndims-1, true, tmp);
			if ( dim < ndims-1 ) { // don't write last iteration as this would write out of bounds.
				nextsize[dim] = size[dim];
				nextstep[dim] = stepDim[dim];
				if (voxelspacing!=NULL) {
					nextvoxelspacing[dim] = voxelspacing[dim];
				}
			}
        }
#ifdef ADDGUARDS // check if variables are not corrupted:
		tmp.checkGuard( nextsize , ndims-1 ); 
		tmp.checkGuard( nextstep , ndims-1 ); 
		if (voxelspacing!=NULL) {
            tmp.checkGuard( nextvoxelspacing , ndims-1);
		}
#endif
    }
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
/*  out = laplaceMulND_c(src, lambda)
*/
 
	// * Check for proper number of arguments. *
	if ((nrhs!=2)&&(nrhs!=3)) {
		mexErrMsgTxt("Two or inputs required: out = laplaceMulND_c(src, lambda [, reciprocalvoxelspacing] ).");
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
    
    mwSize cumsz = 1;
    mwSize minsz = sizein[0];
    for (int dim = 0; dim<ndims;dim++){
        cumsz *= sizein[dim];
        if (minsz>sizein[dim])
            minsz=sizein[dim];
    }
    plhs[0] = mxCreateNumericArray(ndims, sizein, imgclass, mxREAL);
    double * outR = mxGetPr(plhs[0]);
    
    mwSize ntempdata = cumsz / minsz;			 // tempiom in laplacefilter_add
	mwSize ntempsize = (3 * (ndims-1)+ ndims);   // stepTmp in laplacefilter_add
												 // nextsize, nextsize in laplacefilter_full
												 // stepDim in main
	mwSize n_tempspace_bytes = mxGetElementSize(prhs[0]) * ntempdata + sizeof(mwSize) * ntempsize + sizeof(double) * ndims;
#ifdef ADDGUARDS // check if variables are not corrupted:
	n_tempspace_bytes += 100; // add some bytes to have enough space for the guards.
#endif
    char * tempspaceptr = (char *) mxMalloc( n_tempspace_bytes );
	tempspace tmp = tempspace(tempspaceptr, n_tempspace_bytes);

    mwSize * stepDim = tmp.get<mwSize>( ndims ); 
    stepDim[0] = 1;
    for (int dim = 1; dim<ndims; dim++) {
        stepDim[dim] = stepDim[dim-1] * sizein[dim-1];
    }
    if (imgclass==mxDOUBLE_CLASS) {
		/* DEBUG:
		for (int k =0; k<cumsz;k++) {
			outR[k]=2;
		} // */
        laplacefilter_full<double>(         inR,           outR, sizein, stepDim,         lambda, voxelspacing, ndims, false, tmp);
    } else if (imgclass==mxSINGLE_CLASS) {
        laplacefilter_full<float>((float *) inR, (float *) outR, sizein, stepDim, (float) lambda, (float *) voxelspacing, ndims, false, tmp);
    } else 
        mexErrMsgTxt("Unsupported image type.");

#ifdef ADDGUARDS // check if variables are not corrupted:
	tmp.checkGuard( stepDim  , ndims ); 
#endif    
	mxFree(tempspaceptr);
}
#endif

