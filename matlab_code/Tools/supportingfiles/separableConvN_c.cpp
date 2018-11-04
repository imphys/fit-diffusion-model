/*
 * This separableConvN_c is a MEX-file for MATLAB.
 * Created by Dirk Poot, Erasmus MC
 * Build with: (add -g for debugging (add -o for optimizations while debugging); -v for verbose mode )
 * mex separableConvN_c.cpp
 * mex separableConvN_c.cpp -v CXXFLAGS="\$CXXFLAGS -mssse3 -mtune=core2" 

 * WITH MULTITHREADING SUPPORT: (compile with MATLAB 2008a, as our MATLAB 2011b does not support multi-threading compilation)
 * mex separableConvN_c.cpp -v CXXFLAGS="\$CXXFLAGS -fopenmp -mssse3 -mtune=core2" LDFLAGS="\$LDFLAGS -fopenmp"
 */

#include "mex.h"
#include "emm_wrap.cpp"
#include <complex>

#ifdef _OPENMP
#define USEOMP
#endif

#ifdef USEOMP
#include <omp.h>
#endif

template<typename imgtype, typename int_t> struct imReaderType1  { imgtype* im;                int_t step_in_Tdim;};
template<typename imgtype, typename int_t> struct imReaderType1c { imgtype* imR; imgtype* imI; int_t step_in_Tdim;};
template<typename imgtype, typename int_t> struct imReaderType2  { imgtype* im;                int_t step_in_Tdim; int_t step_in_sdim;};
template<typename imgtype, typename int_t> struct imReaderType2c { imgtype* imR; imgtype* imI; int_t step_in_Tdim; int_t step_in_sdim;};

template<typename imgtype, typename int_t> imgtype imReaderFun0(int_t i, const imgtype * im) {
	return im[i];
};
template<typename imgtype, typename int_t> imgtype imReaderFun1(int_t i, const imReaderType1<imgtype,int_t> v) {
	return v.im[i*v.step_in_Tdim];
};
template<typename imgtype, typename filttype, typename int_t> std::complex<filttype> imReaderFun1c(int_t i, const imReaderType1c<imgtype,int_t> v) {
	return std::complex<filttype>(v.imR[i*v.step_in_Tdim], v.imI[i*v.step_in_Tdim]);
};
template<typename imgtype, typename int_t, int vlen> inline vec<imgtype,vlen> imReaderFunvec1(int_t i, const imReaderType1<imgtype,int_t> v) {
	return vec<imgtype,vlen>(v.im+i*v.step_in_Tdim);
};
template<typename imgtype, typename int_t, int vlen> inline vec<imgtype,vlen> imReaderFunvec2(int_t i, const imReaderType2<imgtype,int_t> v) {
	return vec<imgtype,vlen>(v.im+i*v.step_in_Tdim,v.step_in_sdim);
};
			
template<typename imgtype, typename int_t> void imWriterFun1(int_t i, imgtype newval, const imReaderType1<imgtype,int_t> v) {
	v.im[i*v.step_in_Tdim] = newval;
};
template<typename imgtype, typename filttype, typename int_t> void imWriterFun1c(int_t i, std::complex<filttype> newval, const imReaderType1c<imgtype,int_t> v) {
	v.imR[i*v.step_in_Tdim] = (imgtype) newval.real();
	v.imI[i*v.step_in_Tdim] = (imgtype) newval.imag();
};

template<typename imgtype, typename int_t, int vlen> void imWriterFunvec1(int_t i, vec<imgtype,vlen> newval, const imReaderType1<imgtype,int_t> v) {
    newval.store(v.im + i*v.step_in_Tdim);
};

template<typename imgtype, typename int_t, int vlen> void imWriterFunvec2(int_t i, vec<imgtype,vlen> newval, const imReaderType2<imgtype,int_t> v) {
    newval.store(v.im + i*v.step_in_Tdim, v.step_in_sdim);
};

template<typename imgtype, typename int_t> void imAccumFun0(int_t i, imgtype newval, imgtype * im) {
	im[i] += newval;
};
template<typename imgtype, typename int_t> void imAccumFun1(int_t i, imgtype newval, const imReaderType1<imgtype,int_t> v) {
	v.im[i*v.step_in_Tdim] += newval;
};
template<typename imgtype, typename filttype, typename int_t> void imAccumFun1c(int_t i, std::complex<filttype> newval, const imReaderType1c<imgtype,int_t> v) {
	v.imR[i*v.step_in_Tdim] += newval.real();
	v.imI[i*v.step_in_Tdim] += newval.imag();
};

template <typename outtype, typename filttype, typename int_t, typename readerargsT, outtype reader(int_t, const readerargsT), typename writerargsT, void writer(int_t, outtype, const writerargsT)>
inline void filterLineWithFusetemp( const readerargsT readerargs, const writerargsT writerargs, 
						 int_t ninTdim ,
						 const filttype * __restrict F, int_t  mF, 
						 outtype * __restrict tempspace) 

/* filters the elements provided by reader with the filter F. 
NOTE: if doAdjoint==true: writer should accumulate the results (so it should be like '+=' )
*/
{
    int m_offset = -(mF-1)/2;
    for( int_t m = 0 ; m< ninTdim; m++ ) {
        tempspace[m] = reader( m, readerargs );
    }
	for( int_t m = 0 ; m< ninTdim; m++ ) {
		// i_start : start position in filter; 
		int_t i_start = 0; 
		if (i_start + m + m_offset < 0) {
            // solve i_start + m + m_offset ==0:
			i_start = -( m + m_offset); // adjust start position due to begin of image 
		}
		int_t i_end = mF;
		if (m + m_offset + i_end > ninTdim) { // does end of filter extend after allowed range in source image
            // solve m + m_offset + i_end == ninTdim
			i_end = ninTdim - (m + m_offset);  // adjust so that end samples max_i, last sample: i_end-1
        } 
		outtype	tempOutR = (outtype) 0.0; 
        for (int_t i = i_start; i< i_end; i++) {   
            tempOutR += tempspace[ m + m_offset + i ] * F[mF-1-i];
        }
        writer(m , tempOutR , writerargs);
	} // end loop m
}

template <typename imgtype, typename int_t >
inline void filterLineWithFusetemp_delegate( imgtype * img, int_t step_Tdim, int_t step_sdim, 
						 int_t ninTdim , int_t ninsdim, 
						 const imgtype * __restrict F, int_t  mF, 
						 imgtype * __restrict tempspace) 
{
    int s =0;
    const int vlen = 4;
    for( ; s< ninsdim-3; s+=vlen ) {
        typedef vec<imgtype,vlen> outtypev;
        if (step_sdim==1) {
            typedef imReaderType1<imgtype, int_t> rwType1;
            rwType1 rw;
            rw.im = img + s * step_sdim;
            rw.step_in_Tdim = step_Tdim;
            filterLineWithFusetemp< outtypev, imgtype, int_t, rwType1, imReaderFunvec1<imgtype, int_t, vlen> , rwType1, imWriterFunvec1<imgtype, int_t, vlen> >
                      ( rw, rw, ninTdim , F, mF, (outtypev * ) tempspace);
        } else {
            typedef imReaderType2<imgtype, int_t> rwType2;
            rwType2 rw;
            rw.im = img + s * step_sdim;
            rw.step_in_Tdim = step_Tdim;
            rw.step_in_sdim = step_sdim;
             filterLineWithFusetemp< outtypev, imgtype, int_t, rwType2, imReaderFunvec2<imgtype, int_t, vlen> , rwType2, imWriterFunvec2<imgtype, int_t, vlen> >
                       ( rw, rw, ninTdim , F, mF, (outtypev * ) tempspace);
        }
        
    }
    for( ; s< ninsdim; s++ ) {
        typedef imReaderType1<imgtype, int_t> rwType;
        rwType rw;
        rw.im = img + s * step_sdim;
        rw.step_in_Tdim = step_Tdim;
        typedef imgtype outtype;
        filterLineWithFusetemp<outtype, imgtype, int_t, rwType, imReaderFun1 , rwType, imWriterFun1>
                  ( rw, rw, ninTdim , F, mF, (outtype*) tempspace);
        
    }
                         
}
#define MaxNrOfDims 6
#ifdef _MSC_VER
typedef signed __int64 size_int;
#else
typedef int64_t size_int;
#endif
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
/*  out = separableConvN_c(in, filt)
//OR: out = separableConvN_c(in, filt, startend)
// Algorithm:
out = in;
for dim=1:ndims(in)
	out = conv(out, filt{dim}, dim); % if conv would support dimension selecting dim argument.
end;
*/
/*	double *outR, *outI, *F;
	mxArray * curFilt;
	const mwSize * sizeCompArea;
	int sizeout[MaxNrOfDims];
	size_int curArPos[MaxNrOfDims];
	size_int cumsizeout[MaxNrOfDims];
	size_int cumsizein[MaxNrOfDims];
	size_int nCompArea, n_curArea, curindxIn,curindxOut;
	size_int  ndims, ndimCompArea, mF, nFsteps;
	mxClassID imgclass;
	mxClassID filtclass;
	bool iscomplex, emptyplane;
	mxComplexity iscomplexArg;
	bool doerror;
	int k, dimFilt, dimSteps, arealp;
	double curpos;*/
 
	// * Check for proper number of arguments. *
	if ((nrhs!=2)&&(nrhs!=2)) {
		mexErrMsgTxt("Two or three inputs required: out = separableConvN_c(in, filt).");
	} else if(nlhs!=1) {
		mexErrMsgTxt("One outputs required.");
	}

	// * parse inputs * /
	//in:
	int ndims = mxGetNumberOfDimensions(prhs[0]);

	int numfilt = mxGetNumberOfElements(prhs[1]);
	if (ndims<numfilt) { 
		mexErrMsgTxt("Too many filters specified.");
	}
	mxClassID imgclass = mxGetClassID(prhs[0]);
	if ((imgclass!=mxDOUBLE_CLASS) & (imgclass!= mxSINGLE_CLASS)) {
		mexErrMsgTxt("Image input should be single or double. (Other types can easily be added)");
	}
	if (mxGetClassID(prhs[1])!=mxCELL_CLASS) {
		mexErrMsgTxt("Filt should be a cell array.");
	}
	plhs[0] = mxDuplicateArray(prhs[0]);
    //prhs[0] = NULL;
//    mxDestroyArray(prhs[0]);
	bool iscomplex= mxIsComplex(plhs[0]);
	const mwSize * sizein = mxGetDimensions(plhs[0]);
	double * inR = mxGetPr(plhs[0]);
	double * inI = NULL;
	if (iscomplex) {
        //mexErrMsgTxt("Complex currently not supported.");
		inI = mxGetPi(plhs[0]);
    }

	int maxsz = 0;
	for(int filtdim =0;filtdim<numfilt;filtdim++) {
		if (maxsz<sizein[filtdim]) 
			maxsz=sizein[filtdim];
	}
#ifdef USEOMP
    int olddynamic =  omp_get_dynamic();
    omp_set_dynamic(1);
    omp_set_num_threads(omp_get_num_procs());	
	int max_use_threads = omp_get_max_threads();
#else
	int max_use_threads = 1;
#endif
	// DEBUG:
	//mexPrintf("number of threads: %d" , max_use_threads);


	int n_tempspace_bbytes = maxsz  * sizeof(vec<double,4>)* max_use_threads;
    char * tempspace_b = (char *) mxMalloc( n_tempspace_bbytes );

    for (int realandim_idx = 0 ; realandim_idx<2; ++realandim_idx ){
        if ( realandim_idx > 0 ) 
            if (iscomplex) {
                inR = inI;
            } else {
                continue;
            }
	for(int filtdim =0;filtdim<numfilt;filtdim++) {
		mxArray * curFilt = mxGetCell(prhs[1],filtdim);
        if (mxIsEmpty(curFilt)) {
            continue;
        }
            
		if (mxGetClassID(curFilt) != imgclass && !mxIsComplex(curFilt)) {
			mexErrMsgTxt("All filters should be non complex and of the same class as the image.");
		}
		mwSize n_sdim = 1;
		for (int tdim = 0 ; tdim< filtdim; tdim++ ) {
			n_sdim *= sizein[tdim];
		}
		mwSize cumsz = 1;
		for (int tdim = filtdim+1; tdim< ndims; tdim++ ) {
			cumsz *= sizein[tdim];
		}
		double * F = mxGetPr(curFilt);
		int mF = mxGetNumberOfElements(curFilt);

        int step_col = sizein[filtdim] * n_sdim;
        int step_sdim = 1;
        int step_Tdim = n_sdim;
		int step_col_in_sdim = 0;

		if (cumsz == 1) {
			// just one column, not efficient for splitting over columns.
			// Therefore, split step_dim in  max_use_threads blocks.
			cumsz  = max_use_threads;
			step_col_in_sdim = ((n_sdim + (max_use_threads-1))/max_use_threads) ;
			step_col = step_col_in_sdim * step_sdim;
		}

        int ninTdim = sizein[filtdim];
		//template < typename imgtype , typename filttype, typename postype, typename int_t >
        int col;
#ifdef USEOMP
#pragma omp parallel for default(none) shared(cumsz, imgclass,inR,step_col,step_Tdim, step_sdim, ninTdim , n_sdim, F, mF,tempspace_b,n_tempspace_bbytes, step_col_in_sdim, max_use_threads) //private()        
#endif
		for (col = 0; col<cumsz; col++ ) {
			int do_n_sdim = n_sdim;
#ifdef USEOMP
			int thread_num = omp_get_thread_num();
			if (step_col_in_sdim!=0) {
				// for OMP we splitted in sdim.
				do_n_sdim = step_col_in_sdim; // reduce number in sdim.
				if ((col+1)*step_col_in_sdim > n_sdim) {
					// fix up last element.
					do_n_sdim = n_sdim - col*step_col_in_sdim;
				}
			}
#else
			int thread_num = 0;
#endif
            int n_tempspace_bytes = n_tempspace_bbytes /max_use_threads;
            char * tempspace = tempspace_b + n_tempspace_bytes * thread_num;
            if ( imgclass==mxDOUBLE_CLASS) {
                filterLineWithFusetemp_delegate( inR + col*step_col, step_Tdim, step_sdim, 
						 ninTdim , do_n_sdim, 
						 F, mF, (double *) tempspace);
            } else { // single class:
                float * inR_ = (float *) inR;
                float * F_   = (float *) F;
                filterLineWithFusetemp_delegate( inR_ + col*step_col, step_Tdim, step_sdim, 
						 ninTdim , do_n_sdim, 
						 F_, mF, (float *) tempspace);
            }
/*			affineTransformCorePlane_delegate(inR + sizein[filtdim] * n_sdim, (inI==NULL? NULL: inI + col * sizein[filtdim] * n_sdim), sizein[filtdim], n_sdim, n_sdim, 1,//imgtype * __restrict inR, imgtype * __restrict inI, int_t n_in_Tdim, int_t n_sdim, int_t step_in_Tdim, int_t step_in_sdim,
				(double) -(mF-1)/2, (double) 1, (double) 0, //postype pos, postype dpos_Tdim, postype dpos_sdim, 
				inR + sizein[filtdim] * n_sdim, (inI==NULL? NULL: inI + col*sizein[filtdim] * n_sdim), sizein[filtdim], n_sdim, 1,// imgtype * __restrict outR, imgtype * __restrict outI, int_t n_out_Tdim, int_t step_out_Tdim, int_t step_out_sdim, 
				area, 2, (double) 0, //postype * compOutArea, int_t n_compOutArea, postype compOutAreaOffset,
				F, mF, 1, //filttype * F, int_t mF, int_t nFsteps,
				area,0, 0, tempspace, n_tempspace, //filttype * __restrict AR, int_t nAR, int_t ARrunoutLen, tempspace, int_t n_tempspace,
				false);*/
		}
	}}
#ifdef USEOMP
	omp_set_dynamic(olddynamic);    
#endif
	mxFree(tempspace_b);
}
