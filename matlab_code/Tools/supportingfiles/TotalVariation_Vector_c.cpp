/*
 * This TotalVariation_Vector_c is a MEX-file for MATLAB.
 * Created by Dirk Poot, Erasmus MC
 * Build with: (add -g for debugging (add -o for optimizations while debugging); -v for verbose mode )
 * mex TotalVariation_Vector_c.cpp -largeArrayDims

 [ TV, grad, hessinfo ] = TotalVariation_Vector_c( x, spacing, weights, offset )
	% see totalVariationVecRegularizer.m for explaination

	TV = sum_i sqrt(  offset^2 + num_n  (x_i-x_n)'* weights * (x_i-x_n) ) - offset.
	% where n = all neigbors of element i.
 *
 *
 * More explicitly defined:
 *  

 */

// #define EXPLICITSPARSEHESSIAN  // if define EXPLICITSPARSEHESSIAN : return explicit hessian as sparse matrix.	
							      //  if not defined                 : return hessian multiplication preparation matrix. 
#define ADDGUARDS 

#include "mex.h"
#include "tempSpaceManaging.cpp"
//#include "emm_wrap.cpp"
//#include <complex>
//#include <omp.h>
#include <vector>
#include <iterator>
#include <math.h>
#include <algorithm>



#ifndef VOXGRADNORMCLASSES
#define VOXGRADNORMCLASSES
using std::iterator_traits;
using std::vector;

template <typename normType, typename vecType> class dist_w_vecw {
private:
	typedef typename iterator_traits< vecType >::value_type vecElType;
	const vecType scale;
	const int veclen;
public:
	typedef mwSize * ThessianIndexVec; 
	static const bool constantHessian = true;
	dist_w_vecw( const vecType scale_, int veclen_) : scale(scale_), veclen(veclen_) { 
//		scale=scale_ ;
//		veclen=veclen_;
	};

	inline normType computeNorm( const vecType g ) {
		normType norm =0;
		for (int k=0;k<veclen;k++) {
			norm += g[k]*g[k]*scale[k];
		};
		return norm;
	};
	inline void derivativeNorm( const vecType g, vecType dg) {
		for (int k=0;k<veclen;k++) {
			dg[k] = 2*g[k]*scale[k];
		};
	};
	int hessianNorm(const vecType g, vecType & Hg, ThessianIndexVec & offsetH, tempspace & tmp) {
		Hg = tmp.get<vecElType>(veclen);
		offsetH = tmp.get<mwSize>(veclen);
		for (int k = 0 ; k < veclen; k++ ) {
			Hg[k] = scale[k]*2;
			offsetH[k] = k*(veclen+1);

		}
		return veclen;
	}
	void hessianMul(const vecType g, const vecType x, vecType Hx) {
		for (int k=0;k<veclen;k++) {
			Hx[k] = 2*x[k]*scale[k];
		};
	}
};



template <typename normType, typename vecType > class dist_w_matw {
private:
	typedef typename iterator_traits< vecType >::value_type vecElType;
	vector< vecElType > mat;
	int veclen;
public:
	typedef mwSize * ThessianIndexVec; 
	static const bool constantHessian = true;
	dist_w_matw( const vecType mat_, int veclen_) {
		veclen = veclen_;
		mat = vector< vecElType >( veclen * veclen ); 
		
		for (int k1 = 0; k1<veclen; k1++) {
			for (int k2=k1+1; k2<veclen; k2++ ) {
				mat[k1 + veclen * k2] = mat_[k1 + veclen * k2] + mat_[k2 + veclen * k1];
			}
			mat[k1 + veclen * k1] = mat_[k1 + veclen * k1];
		}
	
	};
	inline normType computeNorm(const vecType g ) {
		normType norm =0;
		vecElType * matb = &mat[0];
		for (int k1=0;k1<veclen;k1++) {
			for (int k2=0; k2<=k1; k2++ ) {
				norm += g[k1]*g[k2]*matb[k2];
			}
			matb += veclen;
		}
		return norm;
	}
	inline void derivativeNorm(const vecType  g, vecType  dg) {
		vecElType * matb = &mat[0];
		for (int k1=0;k1<veclen;k1++) {
			dg[k1] = 0;
			for (int k2=0; k2<=k1; k2++ ) {
				dg[k1] += g[k2]*matb[k2];
				dg[k2] += g[k1]*matb[k2];
			}
			matb += veclen;
		}
	}
	int hessianNorm(const vecType g, vecType & Hg, ThessianIndexVec & offsetH, tempspace & tmp) {
		Hg = tmp.get<vecElType>(veclen*veclen);
		offsetH = tmp.get<mwSize>(veclen*veclen);
		vecElType * matb = &mat[0];
		for (int k=0;k<veclen*veclen;k++) {
			Hg[k]=0;
			offsetH[k] = k;
		}
		for (int k1=0;k1<veclen;k1++) {
			for (int k2=0; k2<=k1; k2++ ) {
				Hg[k1*veclen+k2] += matb[k2];
				Hg[k2*veclen+k1] += matb[k2];
			}
			matb += veclen;
		}
		return veclen*veclen;
	}
	void hessianMul(const vecType g, const vecType x, vecType Hx) {
		vecElType * matb = &mat[0];
		for (int k1=0;k1<veclen;k1++) {
			Hx[k1] = 0;
			for (int k2=0; k2<=k1; k2++ ) {
				Hx[k1] += matb[k2] * x[k2];
				Hx[k2] += matb[k2] * x[k1];
			}
			matb += veclen;
		}
	}
};


template< typename base > class assignmenthelper {
public:
	typedef assignmenthelper< base > self;
	typedef typename base::int_T int_t ;
	typedef typename base::value_type value_type;
private:
	base * obj;
	int_t idx;

public:
	assignmenthelper( base & obj_, int_t idx_) : obj(&obj_), idx(idx_) {};
	inline void operator=(value_type newval) {
		obj->write(idx,newval);
	}
	inline void operator=(self newval) {
		obj->write(idx,newval.obj->read(newval.idx));
	}
	/*inline void operator+=(value_type newval) {
		obj->accumulate(idx,newval);
	}*/
	/*inline value_type operator*() {
		return obj.read(idx);
	}*/
	inline operator value_type() const { 
		return obj->read(idx); 
	};
	inline void swap(const self & it2) {
		//if (obj==it2.obj) {
			obj->swap(idx, it2.idx);
		/*} else {
			obj->swap(idx, it2.obj, it2.idx);
		}*/
	}
};
namespace std
{
    template<typename base>
    void swap(assignmenthelper<base>& lhs, assignmenthelper<base>& rhs)
    {
       lhs.swap(rhs);
    }
}

template< typename base > class iteratorhelper {
public:
	typedef typename std::random_access_iterator_tag iterator_category;
	typedef typename base::value_type value_type;
    typedef ptrdiff_t difference_type;
    typedef ptrdiff_t distance_type;
    typedef typename base::pointer pointer;
    typedef typename base::reference reference;

	typedef iteratorhelper< base > self;
	typedef typename base::int_T int_t ;
private:
	base * obj;
public:
	int_t idx;
public:
	iteratorhelper() {};
	iteratorhelper( base & obj_, int_t idx_) : obj(&obj_), idx(idx_) {};
/*		obj = obj_;
		idx = idx_;
	} // */
	inline self& operator++( ) { // prefix form (increment first, then fetch).
		idx++;
		return *this;
	};
	inline self operator++( int ) { // postfix form (fetch, then increment).
		self ans = *this;
		++*this;
		return ans;
	};
	
	inline self& operator--( ) { //prefix
		idx--;
		return *this;
	};
	inline self operator--(int) { //postfix
		self ans = *this;
		idx--;
		return ans;
	};
	inline self operator+( mwSize n ) {
		return self( *obj, idx+n);
	};
	inline self operator-( mwSize n) {
		return self( *obj, idx-n);
	};
	/*inline mwSize operator+( self b ) {
		return idx+b.idx;
	};*/
	inline mwSize operator-( self b) {
		return idx-b.idx;
	};
	inline void operator+=( mwSize n ) {
		idx += n;
	};
	inline void operator-=( mwSize n ) {
		idx -= n;
	};
	inline bool operator==(const self& other) {
		return (other.obj==obj)&& (other.idx==idx);
	};
	inline bool operator!= (const self& other) const {
		return (other.obj!=obj)|| (other.idx!=idx);
	};
	reference operator*() {
		return assignmenthelper<base>(*obj, idx);
	};
	reference operator[](mwSize n) {
		return assignmenthelper<base>(*obj, idx+n);
	};
	inline bool operator>( self other) {
		return idx > other.idx;
	}
	inline bool operator>=( self other) {
		return idx >= other.idx;
	}
	inline bool operator<( self other) {
		return idx < other.idx;
	}
	inline bool operator<( mwSize n) {
		return idx < n;
	}
	inline bool operator>( mwSize n) {
		return idx > n;
	}
	inline bool operator<=( self other) {
		return idx <= other.idx;
	}

};



template<typename base> inline bool operator<(mwSize n, iteratorhelper<base> b ) {
	return n<b.idx;
};

template<typename T1, typename T2> struct twoElStruct {
	T1 v1; T2 v2;
	twoElStruct(T1 v1_, T2 v2_): v1(v1_), v2(v2_) {};
};

template<typename T1, typename T2, typename int_t> class linked2vec   { 
private:
	T1 v1;
	T2 v2;
	mwSize len;
public:
	typedef linked2vec<T1, T2, int_t> self;
	typedef typename std::random_access_iterator_tag iterator_category;
	typedef twoElStruct< typename iterator_traits< T1 >::value_type, typename iterator_traits< T2 >::value_type> value_type;

	typedef ptrdiff_t difference_type;
    typedef difference_type distance_type;
    typedef int * pointer;
    typedef assignmenthelper<self> reference;
    typedef iteratorhelper<self> iterator;
	typedef int_t int_T;

	linked2vec( T1 v1_, T2 v2_, mwSize len_) {
		v1 = v1_; v2 = v2_;len = len_;
	}
	iterator begin(){
		return iterator(*this, 0);
	}
	iterator end(){
		return iterator(*this, len);
	}
	inline void swap(int_t i1, int_t i2) {
		std::swap( v1[i1], v1[i2] );
		std::swap( v2[i1], v2[i2] );
	}
	inline value_type read(int_t i) {
		return value_type( v1[i], v2[i]);
	};
	inline void write(  int_t i, const value_type newval) {
		v1[i] =  newval.v1;
		v2[i] =  newval.v2;
	};
	inline value_type operator[](int_t i) const {
		return read(i);
	};
	inline reference operator[](int_t i) {
		return reference(*this, i);
	};
};

template< typename T> bool sparseIdxSorter( T el1, T el2 ) {
	return (el1.v2 < el2.v2);
}
#endif


#ifdef ANISOTROPICVOXELS
    #define VOXELSPACINGARGS , const imgtype * voxelspacing
    #define VOXELSPACING2(dim) voxelspacing[ dim ]
#else
    #define VOXELSPACINGARGS 
    #define VOXELSPACING2(dim) 1
#endif
#ifdef EXPLICITSPARSEHESSIAN
	#define HESSIANARGS , double * hess_r, mwSize * hess_i,  mwSize * hess_jc, mwSize numHessEl , double sqrtHessModifScaleFact, 
#else 
	#define HESSIANARGS , double * hess_r
#endif
#ifdef GRADFILTER // defined in TotalVariation_Vector_gradfilt_fc.cpp
    #define GRADFILTERARGS , double * gradfilt, int gradfiltlen, int gradfiltoffset
#else
    #define GRADFILTERARGS 
#endif
template <typename TVtype, typename imgtype, typename distMeasureType, bool skipborder >
inline TVtype TV_core( const imgtype * __restrict source, imgtype * __restrict grad, 
					   const mwSize * __restrict size, const ptrdiff_t * __restrict stepSrc,
                       TVtype offset  VOXELSPACINGARGS , int ndims ,
					   distMeasureType & dist  HESSIANARGS GRADFILTERARGS ,  tempspace tmp)
{
#ifdef ONLYFORWARDGRADIENT // defined in TotalVariation_Vector_fc.cpp and TotalVariation_Vector_gradfilt_fc.cpp
    const bool incudebackgrad = false;
#else
    const bool incudebackgrad = true;
#endif
#ifndef GRADFILTER    
    const int gradfiltlen = 2;
    const int gradfiltoffset = 0;
#endif
    TVtype sqrtHessModifScaleFact;
#ifdef EXPLICITSPARSEHESSIAN
	const bool doExplicitHessian = (hess_r!=NULL);
	const bool doHessMulPrepare = false;
    sqrtHessModifScaleFact *= -.5;
#else
	const bool doExplicitHessian = false;
	const bool doHessMulPrepare = (hess_r!=NULL);
	mwSize * hess_i;
	mwSize * hess_jc;
	mwSize numHessEl ;
#endif 
	TVtype TV = 0;
	TVtype offset2 = offset * offset;

	mwSize * cumSpatialSize = tmp.get<mwSize>( ndims );
	TVtype ** gradcache = tmp.get<TVtype *>( ndims );
	cumSpatialSize[0] = 1;
	mwSize nel = 1;
	for (int dim = 0; dim < ndims; dim++ ) {
		nel *= size[dim];
	}
	for (int dim = 1; dim < ndims; dim++ ) {
//         if (skipborder) {
//             cumSpatialSize[dim] = cumSpatialSize[dim-1] * (size[dim]-1);// if skipborder : in iteration only last element in each dimension is skipped; first is not added but has to be computed for backward gradient.
//         } else {
            cumSpatialSize[dim] = cumSpatialSize[dim-1] * size[dim];
//         }
		if (incudebackgrad) {
		gradcache[dim] = tmp.get<TVtype>( cumSpatialSize[dim-1] );
		for (int k =0 ; k<cumSpatialSize[dim-1]; k++ ){
			gradcache[dim][ k ] = 0;
		}
		}
	}
	imgtype * g = tmp.get<imgtype>( size[0] );
	//imgtype * dg = NULL;
	imgtype * G_g = NULL;
	imgtype * G_dg = NULL;
	mwSize * G_steps = NULL;
#ifdef ANISOTROPICVOXELS
	int * G_dim = NULL;
#endif
	imgtype * Hg;
	mwSize * offsetH;
	mwSize numHg  = 0;
	mwSize curhessidx =0;
	if (grad!=NULL) {
		//dg = tmp.get<imgtype>( size[0] );
		G_g = tmp.get<imgtype>( ( doExplicitHessian ? 2* (ndims-1) :1)*size[0]);
		G_dg= tmp.get<imgtype>( ( doExplicitHessian ? 2* (ndims-1) :1)*size[0]);
		if ( doExplicitHessian ) {
			G_steps = tmp.get<mwSize >( 2* ndims );
			if (dist.constantHessian) {
				numHg = dist.hessianNorm( NULL, Hg, offsetH, tmp);
			}
#ifdef ANISOTROPICVOXELS
			G_dim = tmp.get<int>( 2* ndims );
#endif
		}
	}
	ptrdiff_t * pos = tmp.get<ptrdiff_t>( ndims );
	ptrdiff_t * curGradCachePos = tmp.get<ptrdiff_t>( ndims );
	for ( int dim =0; dim<ndims; dim++ ) {
		pos[dim]=0;
		curGradCachePos[dim]=0;
	};
	mwSize curpos = 0;
	for (int k = 0; k < cumSpatialSize[ndims-1]; k++ ) {
		TVtype TVk = offset2;
		for ( int dim =1; dim<ndims; dim++ ) {
			if (incudebackgrad) {
				TVk += gradcache[dim][ curGradCachePos[ dim ] ];
			}
			if ( (pos[dim]+gradfiltoffset>=0) &&  pos[ dim ]+gradfiltoffset+gradfiltlen <= ptrdiff_t( size[ dim ] ) ) { // check if we can do forward derivative.
#ifdef GRADFILTER 
				for (int d = 0; d<size[0]; d++ ){
                    g[d] = 0;
                    for (int i = 0; i < gradfiltlen;++i) {
                        g[d] += gradfilt[i] * source[ curpos + d + (gradfiltoffset+i)*stepSrc[ dim ] ];
                    }
				}
#else
				for (int d = 0; d<size[0]; d++ ){
					g[d] = source[ curpos +d ] - source[ curpos + d + stepSrc[ dim ] ];
				}
#endif                
				TVtype nrm2 = dist.computeNorm( g ) * VOXELSPACING2( dim );

				if (incudebackgrad) {
					gradcache[dim][ curGradCachePos[ dim ] ] = nrm2;
				}
				TVk += nrm2;
			}
			curGradCachePos[ dim ]++;
		}
        bool doadd = true;
        if (skipborder) {
            for ( int dim =1; dim<ndims; dim++ ) {
                if (pos[dim]==0) { // only need to test pos[dim]==0, since size[dim]-1 is skipped in loop control.
                    doadd = false;
                    break;
                };
            }
        }
        if (doadd) {
		TV += (sqrt(TVk)-offset);

		// compute gradient (if needed):
		if (grad!=NULL) {
			// d TV / d source_i = sum_kj d TV / d TVk * d TVk / d nrm2_j * d nrm2_j / d g_j * d g_j / d source_i
			//
			// d TV / d TVk = .5/sqrt(TVk) 
			// d TVk / d nrm2_j = 1;    if j selects one of the gradients involving k, 0 otherwise.
			// d nrm2_j / d g_j =  dist.derivativeNorm( g_j ) *  VOXELSPACING2( dim )  (== dg_j)
			// d g_j/ d source_i =  1;  if i == k 
			//				     = -1;  if i is neighbor of k
			//					 = 0 otherwise.
			//
			// => d TV / d source_i = .5/sqrt(TVk) *  dg * (1 or -1)
			TVtype rTVk = .5/sqrt(TVk);
			if (doHessMulPrepare) {
				hess_r[k] = rTVk ;
			}
			int numgrads = 0;
			imgtype * g_ = G_g;
			imgtype * dg_ = G_dg;
			for ( int dim =1; dim<ndims; dim++ ) {

				if (incudebackgrad && (pos[dim]+gradfiltoffset-1>=0) &&  (pos[ dim ]+gradfiltoffset+gradfiltlen-1 <= ptrdiff_t( size[ dim ] ) ) ) { // check if we can do backward derivative.
#ifdef GRADFILTER
                    for (int d = 0; d<size[0]; d++ ) {
                        g_[d] = 0;
                        for (int i = 0; i < gradfiltlen; ++i ) {
                            g_[d] += gradfilt[i] * source[ curpos + d + (gradfiltoffset+i-1)*stepSrc[ dim ] ];
                        }
                    }
#else
					for (int d = 0; d<size[0]; d++ ){
						g_[d] = source[ curpos + d ] - source[ curpos +d  - stepSrc[ dim ] ];
					}
#endif                    
					dist.derivativeNorm( g_ , dg_ ) ;
					for (int d = 0; d<size[0]; d++ ){
						dg_[d] *= VOXELSPACING2( dim );
#ifdef GRADFILTER
                        imgtype dg_d = dg_[d] * rTVk;
                        for (int i = 0; i < gradfiltlen;++i) {
                            grad[ curpos + d+ (gradfiltoffset+i-1)*stepSrc[ dim ] ] += gradfilt[i] * dg_d ;
    
                        }
#else
						grad[ curpos + d                 ] += dg_[d] * rTVk ;
						grad[ curpos + d -stepSrc[ dim ] ] -= dg_[d] * rTVk ;
#endif
					}
					if ( doExplicitHessian ) {
						G_steps[numgrads] = -stepSrc[ dim ];
#ifdef ANISOTROPICVOXELS
						G_dim[numgrads] = dim;
#endif
						numgrads++;g_ += size[0]; dg_ += size[0];
					}
				}
				if ( (pos[dim]+gradfiltoffset>=0) &&  (pos[ dim ]+gradfiltoffset+gradfiltlen <= ptrdiff_t( size[ dim ] ) ) ) { // check if we can do forward derivative.
#ifdef GRADFILTER
                    for (int d = 0; d<size[0]; d++ ){
                        g_[d] = 0;
                        for (int i = 0; i < gradfiltlen; ++i ) {
                            g_[d] += gradfilt[i] * source[ curpos + d + (gradfiltoffset+i)*stepSrc[ dim ] ];
                        }
                    }
#else                    
					for (int d = 0; d<size[0]; d++ ){
						g_[d] = source[ curpos + d ] - source[ curpos + d + stepSrc[ dim ] ];
					}
#endif
					dist.derivativeNorm( g_ , dg_) ;
					for (int d = 0; d<size[0]; d++ ){
						dg_[d] *= VOXELSPACING2( dim );
#ifdef GRADFILTER
                        imgtype dg_d = dg_[d] * rTVk;
                        for (int i = 0; i < gradfiltlen; ++i ) {
                            grad[ curpos + d+ (gradfiltoffset+i)*stepSrc[ dim ] ] += gradfilt[i] * dg_d ;
                        }
#else                        
						grad[ curpos + d                 ] += dg_[d] * rTVk ;
						grad[ curpos + d +stepSrc[ dim ] ] -= dg_[d] * rTVk ;
#endif                        
					}
					if ( doExplicitHessian ) {
						G_steps[numgrads] = stepSrc[ dim ];
#ifdef ANISOTROPICVOXELS
						G_dim[numgrads] = dim;
#endif
						numgrads++;g_ += size[0]; dg_ += size[0];
					}
				}
			}
			// Hessian & Multiply vector x with
			// d2 TV / d source_i d source_j =  sum_klm   d2 TV / d TVk d TVk * (d TVk / d nrm2_l * d nrm2_l / d g_l * d g_l/ d source_i) * (d TVk / d nrm2_m * d nrm2_m / d g_m * d g_m/ d source_j) 
			//									         + d TV / d TVk        *  d TVk / d nrm2_l * d2 nrm2_l/d g_l d g_l' * (d g_l/ d source_i) * (d g_l'/ d source_j )  

			// d2 TV / d TVk d TVk = -.25/(TVk * sqrt(TVk))
			// d2 nrm2/d g  d g  y    = dist.HessianNormMul( g , y)*  VOXELSPACING2( dim )  ( == Hg * y )
			// s_li =  1  if  k==i, -1 if k is neighbor of i , 0 otherwise
			//
			// => d2 TV / d source_i d source_j =  sum_klm  -.25/(TVk * sqrt(TVk)) * dg_l * s_li  *  dg_m' * s_mj 
			//									  +sum_kl    .5/sqrt(TVk)    *   s_li * Hg * s_lj 
			//
			// => sum_j d2 TV / d source_i d source_j * x_j = sum_kl ( s_li * dg_l *  (-.25/(TVk * sqrt(TVk)))  * sum_mj dg_m' * s_mj * x_j )
			//												  +sum_kl ( s_li * .5/sqrt(TVk)  * Hg * sum_j s_lj * x_j )

			if ( doExplicitHessian ) {
#ifdef EXPLICITSPARSEHESSIAN                
#ifdef GRADFILTER
                #error "cannot create explicit hessian for arbitrary gradient filter"
#endif
#endif        
				// store hessian explicitly (as sparse matrix).
				mwSize kkbase = curhessidx;
				mwSize hstep = size[0]*size[0];
				TVtype r2TVk= sqrtHessModifScaleFact *rTVk/TVk;
				//curhessidx += (2*numgrads+1)*hstep;
				// 'allocate' part of kk, ki, jk. Allocate once and re-use, since each of ki and jk is used up to 2*ndims times, and kk might be used (2*ndims)^2 times.
				for (int i=0; i< 2*numgrads+1 ; i++ ) {
					for (int i_i = 0; i_i < size[0]; i_i++ ) {
						for (int j_i = 0 ; j_i<size[0]; j_i++ ) {
							ptrdiff_t ii = curpos + (  ((i>0)&(i<numgrads+1)) ? G_steps[i-1]          : 0) + i_i;
							ptrdiff_t jj = curpos + (  ( i>numgrads         ) ? G_steps[i-numgrads-1] : 0) + j_i;
							hess_i[ curhessidx ] =  jj + nel* ii ;
							hess_r[ curhessidx ]=0; //However, it should already be, so we shouldnt need to set it.
							curhessidx++;
						}
					}
				}
				// store dg * dg' part of hessian.
				for (int i=0; i< numgrads ; i++ ) {
					imgtype * dg_i = G_dg + i *  size[0];
					ptrdiff_t iibase;
					for (int j =0 ; j<numgrads ; j++ ) {
						if (j==i) {
							iibase = curhessidx;
						}
						imgtype * dg_j = G_dg + j *  size[0];
						for (int i_i = 0; i_i < size[0]; i_i++ ) {
							for (int j_i = 0 ; j_i<size[0]; j_i++ ) {
								TVtype temp = r2TVk * dg_i[i_i] * dg_j[j_i] ;
								hess_r[ kkbase + j_i + size[0] * i_i ] += temp;
								hess_r[ kkbase + (i+1         )*hstep + j_i + size[0] * i_i ] -= temp;
								hess_r[ kkbase + (j+1+numgrads)*hstep + j_i + size[0] * i_i ] -= temp;

								ptrdiff_t ii = curpos + G_steps[i] + i_i ;
								ptrdiff_t jj = curpos + G_steps[j] + j_i;
								
								hess_i[curhessidx] = jj + nel * ii ;
								hess_r[curhessidx] = temp;
								curhessidx++;
							}
						}
					}
					if (!dist.constantHessian) {
						tempspace tmpcopy = tmp;
						numHg = dist.hessianNorm( G_g + i * size[0], Hg, offsetH, tmpcopy );
					}

					// store Hg part of hessian:
					for (int i_i = 0; i_i < numHg; i_i++ ) {
						TVtype temp = rTVk * Hg[i_i] * VOXELSPACING2( G_dim[i] ) ;
						hess_r[ kkbase +                         offsetH[ i_i] ] += temp;
						hess_r[ kkbase + (i+1         )*hstep +  offsetH[ i_i] ] -= temp;
						hess_r[ kkbase + (i+1+numgrads)*hstep +  offsetH[ i_i] ] -= temp;
						hess_r[ iibase +                         offsetH[ i_i] ] += temp;
					}

				}
	
			}

		}
        } // end if doadd

		// Increment position, both linear index as well as voxel coordinate 
		// reset curGradCachePos for the elements that need it.
		{
		int dim =1;
		bool dorepeat;
		do {
			pos[dim]++;
			curpos += stepSrc[ dim ];
			curGradCachePos[ dim ] =0;
            if (skipborder && (pos[dim]==(size[dim]-1))) {
                pos[dim]++;
                curpos += stepSrc[ dim ];
                k += cumSpatialSize[ dim - 1 ];
            }
			dorepeat = ( pos[dim] >= ptrdiff_t( size[dim] ) );
			if (dorepeat) {
				if (incudebackgrad) {
					for (int k2 =0 ; k2<cumSpatialSize[dim-1]; k2++ ) {
						gradcache[dim][k2] = 0;
					}
				}
				curpos -= pos[dim] * stepSrc[ dim ];
				pos[dim]=0;
				dim++;
				dorepeat = dim<ndims;
			}
		} while (dorepeat);
		};
	};

	if (doExplicitHessian) {
		if (curhessidx >= numHessEl) {
			mexErrMsgTxt("Sorry, I've written out of the bounds of the sparse hessian matrix. Memory is corrupted.");
		}
		using std::sort;
		hess_r[curhessidx] = 0;
		hess_i[curhessidx] = nel*nel; // set to one after last column, to fill the last Jc
		linked2vec< TVtype *, mwSize *, mwSize > hessv( hess_r, hess_i, curhessidx );

		sort( hessv.begin(), hessv.end() , sparseIdxSorter< typename linked2vec< TVtype *, mwSize *, mwSize >::value_type > );
		hessv.begin();
		mwSize stidx = 0;
		mwSize rdidx = 1;
		mwSize curcolidx = 1;
		double rnel = 1.0/( (double) nel);
		while (rdidx<=curhessidx) {
			if (hess_i[rdidx-1]>hess_i[rdidx]) { // DEBUG: check sorting.
				mexErrMsgTxt("Bad sorting.");
			}
			if (hess_i[rdidx]!=hess_i[stidx]) {
				// translate hess_i to just i index:
				mwSize j = hess_i[stidx] / nel;			// column index of last store.
				mwSize i = hess_i[stidx] - nel * j;		// row index of last store.
				hess_i[stidx] = i;
				while (curcolidx<=j) {
					hess_jc[curcolidx] = stidx;				// set column counts.
					curcolidx++;
				}
				stidx++;
				hess_i[stidx] = hess_i[rdidx];
				hess_r[stidx] = hess_r[rdidx];
			} else {
				hess_r[stidx] += hess_r[rdidx];
			}
			rdidx++;
		}
		hess_jc[curcolidx] = stidx; 
	}

#ifdef ADDGUARDS 

	tmp.checkGuard( cumSpatialSize, ndims);
    if (incudebackgrad) {
        tmp.checkGuard( gradcache,ndims );
        for (int dim = 1; dim < ndims; dim++ ) {
            tmp.checkGuard( gradcache[dim] , cumSpatialSize[dim-1] );
        }
    }
	tmp.checkGuard( g , size[0] );
	if (grad!=NULL) {
		tmp.checkGuard( G_g , (doExplicitHessian ? 2* (ndims-1) :1)*size[0] );
		tmp.checkGuard( G_dg, ( doExplicitHessian ? 2* (ndims-1) :1)*size[0]);
		if ( doExplicitHessian ) {
			tmp.checkGuard( G_steps , 2* ndims );
			tmp.checkGuard( Hg, numHg );
			tmp.checkGuard( offsetH, numHg );
#ifdef ANISOTROPICVOXELS
			tmp.checkGuard( G_dim , 2* ndims );
#endif
		}
	}
	tmp.checkGuard(  pos , ndims );
	tmp.checkGuard( curGradCachePos , ndims );
#endif 

	return TV;
};


template <typename TVtype, typename imgtype, typename distMeasureType >
inline void TV_hessianMul( const imgtype * __restrict x      ,  imgtype * __restrict Hx, 
						   const imgtype * __restrict source ,  const TVtype * __restrict TVk, 
					       const mwSize * __restrict size    ,  const ptrdiff_t * __restrict stepSrc
                           VOXELSPACINGARGS                  ,  int ndims ,
                           imgtype sqrtHessModifScaleFact    ,  // scale factor for the adjustment of the hessian due to the sqrt operation.
                                                                // For fully accurate hessian it should be 1
                                                                // however optimization convergence might be faster when it's lower. 
					       distMeasureType & dist               GRADFILTERARGS ,  tempspace tmp)
{
#ifdef ONLYFORWARDGRADIENT
    const bool incudebackgrad = false;
#else
    const bool incudebackgrad = true;
#endif	
	mwSize * cumSpatialSize = tmp.get<mwSize>( ndims );
	cumSpatialSize[0] = 1;
	mwSize nel = 1;
	for (int dim = 0; dim < ndims; dim++ ) {
		nel *= size[dim];
	}
	for (int dim = 1; dim < ndims; dim++ ) {
		cumSpatialSize[dim] = cumSpatialSize[dim-1] * size[dim];
	}
	
	imgtype * G_gs = tmp.get<imgtype>( 2* (ndims-1) *size[0]);
	imgtype * G_dg = tmp.get<imgtype>( 2* (ndims-1) *size[0]);
	imgtype * G_gx = tmp.get<imgtype>( size[0] );
	imgtype * Hx_i = tmp.get<imgtype>( size[0] );
	mwSize * G_steps = tmp.get<mwSize >( 2* ndims );
	
#ifdef ANISOTROPICVOXELS
	int * G_dim = tmp.get<int>( 2* ndims );
#endif

	mwSize * pos = tmp.get<mwSize>( ndims );
	for ( int dim =0; dim<ndims; dim++ ) {
		pos[dim]=0;
	};

	mwSize curpos = 0;
    sqrtHessModifScaleFact *= -2;
	for (int k = 0; k < cumSpatialSize[ndims-1]; k++ ) {
		
		TVtype rTVk = TVk[k];		// = .5/sqrt(TVk);
		TVtype r2TVk= sqrtHessModifScaleFact * rTVk * rTVk * rTVk;  // = -.25/( sqrt(TVk)*TVk );

		int numgrads = 0;
		for ( int dim =1; dim<ndims; dim++ ) {
			if (incudebackgrad && ( pos[ dim ] > 0 )) { // check if we can do backward derivative.
				G_steps[numgrads] = -stepSrc[ dim ];
#ifdef ANISOTROPICVOXELS
				G_dim[numgrads] = dim;
#endif
				numgrads++;
			}
			if ( pos[ dim ] < size[ dim ]-1 ) { // check if we can do forward derivative.
				G_steps[numgrads] = stepSrc[ dim ];
#ifdef ANISOTROPICVOXELS
				G_dim[numgrads] = dim;
#endif
				numgrads++;
			}
		}
		// Hessian & Multiply vector x with
		// d2 TV / d source_i d source_j =  sum_klm   d2 TV / d TVk d TVk * (d TVk / d nrm2_l * d nrm2_l / d g_l * d g_l/ d source_i) * (d TVk / d nrm2_m * d nrm2_m / d g_m * d g_m/ d source_j) 
		//									         + d TV / d TVk        *  d TVk / d nrm2_l * d2 nrm2_l/d g_l d g_l' * (d g_l/ d source_i) * (d g_l'/ d source_j )  

		// d2 TV / d TVk d TVk = -.25/(TVk * sqrt(TVk))
		// d2 nrm2/d g  d g  y    = dist.HessianNormMul( g , y)*  VOXELSPACING2( dim )  ( == Hg * y )
		// s_li =  1  if  k==i, -1 if k is neighbor of i , 0 otherwise
        // if GRADFILTER: s_li = 
		//
		// => d2 TV / d source_i d source_j =  sum_klm  -.25/(TVk * sqrt(TVk)) * dg_l * s_li  *  dg_m' * s_mj 
		//									  +sum_kl    .5/sqrt(TVk)    *   s_li * Hg * s_lj 
		//
		// => sum_j d2 TV / d source_i d source_j * x_j = sum_kl ( s_li * dg_l *  (-.25/(TVk * sqrt(TVk)))  * sum_mj dg_m' * s_mj * x_j )
		//												  +sum_kl ( s_li * .5/sqrt(TVk)  * Hg * sum_j s_lj * x_j )

		TVtype dg = 0;
		// apply dg * dg' part of hessian.
		for (int i=0; i< numgrads ; i++ ) {
			imgtype * gs_i = G_gs + i *  size[0];
			imgtype * dg_i = G_dg + i *  size[0];
			imgtype * gx_i = G_gx;// + i *  size[0];
			for (int d = 0; d<size[0]; d++ ){
				gs_i[d] = source[ curpos + d ] - source[ curpos +d + G_steps[i] ];
				gx_i[d] = x[      curpos + d ] - x[      curpos +d + G_steps[i] ];
			}
			dist.derivativeNorm( gs_i , dg_i ) ;

			// intermezzo Hessian of distance * gradient_TV:
			dist.hessianMul( gs_i , gx_i, Hx_i ) ;
			for (int d = 0; d<size[0]; d++ ){
				TVtype temp = Hx_i[d] * rTVk * VOXELSPACING2( G_dim[i] );
				Hx[curpos + d              ] += temp;
				Hx[curpos + d + G_steps[i] ] -= temp; 
			}	

			// dg_m' * s_mj * x_j
			TVtype dgm = 0;
			for (int d = 0; d<size[0]; d++ ){
				dgm += dg_i[d] * gx_i[d];
			}
			dg += dgm * VOXELSPACING2( G_dim[i] ) ;
		}
		dg *= r2TVk;

		for (int j =0 ; j<numgrads ; j++ ) {
			imgtype * dg_j = G_dg + j *  size[0];
			TVtype dgl = dg * VOXELSPACING2( G_dim[j] ); 
			for (int d = 0; d<size[0]; d++ ){
				TVtype temp = dgl * dg_j[d];
				Hx[curpos + d              ] += temp;
				Hx[curpos + d + G_steps[j] ] -= temp; 
			}	
		}




		// Increment position, both linear index as well as voxel coordinate 
		// reset curGradCachePos for the elements that need it.
		{
		int dim =1;
		bool dorepeat;
		do {
			pos[dim]++;
			curpos += stepSrc[ dim ];
			dorepeat = ( pos[dim] >= size[dim] );
			if (dorepeat) {
				curpos -= pos[dim] * stepSrc[ dim ];
				pos[dim]=0;
				dim++;
				dorepeat = dim<ndims;
			}
		} while (dorepeat);
		};
	};

	
}


#undef VOXELSPACINGARGS 
#undef VOXELSPACING2
#undef HESSIANARGS
#undef GRADFILTERARGS

#ifndef ANISOTROPICVOXELS    
#define ANISOTROPICVOXELS
#include __FILE__
#undef ANISOTROPICVOXELS

#ifndef SKIPBORDER
#define SKIPBORDER false
#endif


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
/*  [TV, grad] = TotalVariation_Vector_c(img, spacing, weights, offset)

	TV = sum_i sqrt(  offset^2 + num_n  (x_i-x_n)'* weights * (x_i-x_n) ) - offset.
	% where n = all neigbors of element i.
*/
 	bool returnSparse = true;// set to false to return components of sparse hessian, instead of a sparse matrix. for DEBUG. 

	// * Check for proper number of arguments. *
#ifdef GRADFILTER
    const int hessinpidx = 6;
	if ((nrhs!=6)&&(nrhs!=8)) {
		mexErrMsgTxt("six or eight inputs required:  [TV, grad, hess] = TotalVariation_Vector_c(img, spacing, weights, offset , gradfilt, gradfiltoffset, [,H, x]).");
	} else 
#else
    const int hessinpidx = 4;
	if ((nrhs!=4)&&(nrhs!=6)) {
		mexErrMsgTxt("Four or six inputs required:  [TV, grad, hess] = TotalVariation_Vector_c(img, spacing, weights, offset [,H, x]).");
	} else 
#endif
    if ((nlhs!=1)&&(nlhs!=2)&&(nlhs!=3)&& returnSparse)  {
		mexErrMsgTxt("One, two, or three outputs required.");
	}
	bool hessmul = (nrhs==(hessinpidx+2));
	bool computeGrad = (nlhs>1);
	bool computeHess = (nlhs>2);
#ifdef EXPLICITSPARSEHESSIAN
	bool computeExplicitHess = computeHess;
	#define HESSIANARGSs , hess_pr, hess_ir, hess_jc, numHessEl, sqrtHessModifScaleFact
#else
	bool computeExplicitHess = false;
	#define HESSIANARGSs , hess_pr
#endif

	// * parse inputs * /
	const mxArray * img     = prhs[0];
	const mxArray * spacing = prhs[1];
	const mxArray * weights = prhs[2];
	const mxArray * mxoffset  = prhs[3];
#ifdef GRADFILTER
	const mxArray * gradfilter  = prhs[4];
	const mxArray * gradfilteroffset  = prhs[5];
    #define GRADFILTERARGSs , gradfilterp, gradfilterlen, gradfilteroffseti
#else
    #define GRADFILTERARGSs 
#endif

	//img:
	int ndims = static_cast<int>( mxGetNumberOfDimensions( img ) ); // number of dimension should be representable by int
	const mwSize * size = mxGetDimensions(img);
	mxClassID imgclass = mxGetClassID(img);
    
	if ((imgclass!=mxDOUBLE_CLASS) & (imgclass!= mxSINGLE_CLASS) ) {
		mexErrMsgTxt("Image input should be single or double. (Other types can easily be added)");
	}
    if (mxIsComplex(img)) {
        mexErrMsgTxt("Image input should non complex");
    }
	double * source = mxGetPr(img);

	// spacing:
    double * voxspacing = NULL;
	if (!mxIsEmpty(spacing)!=0) {
		if ( (imgclass!= mxGetClassID(spacing)) || (mxGetNumberOfElements(spacing)!=ndims)) {
			mexErrMsgTxt("Reciprocal Voxel Spacing should be of same class as image and have ndims elements.");
		}
		voxspacing = mxGetPr( spacing );
    }
    
	// weights:
	const mwSize * sizeweight = mxGetDimensions(weights);
	if ( (mxGetNumberOfDimensions( weights )!=2) || (sizeweight[0]!=size[0]) || ((sizeweight[1]!=size[0]) && (sizeweight[1]!=1)) || (imgclass!= mxGetClassID(weights))) {
		mexErrMsgTxt("weights should be column vector or square matrix with size(img,1)  (i.e. the vector length) and data type of weights should be equal to data type of image (single/double).");
	}
	double * scale = mxGetPr( weights );
	bool isvecw = sizeweight[1] == 1;


	// offset
	if (mxGetNumberOfElements(mxoffset)!=1) { 
		mexErrMsgTxt("Offset should be scalar.");
	}
	double offset = mxGetScalar(mxoffset);

#ifdef GRADFILTER
    // gradilter
	const mwSize * sizegradfilter = mxGetDimensions(gradfilter);
	if ( (mxGetNumberOfDimensions( gradfilter )!=2) || (sizegradfilter[0]!=1) || (imgclass!= mxGetClassID(gradfilter))) {
		mexErrMsgTxt("gradilter should be row vector and data type should be equal to data type of image (single/double).");
	}
	double * gradfilterp = mxGetPr( gradfilter );
	int gradfilterlen = static_cast<int>( sizegradfilter[1] ) ; // gradfilter should be short, certainly not larger than maximum integer.
    
    // gradfilteroffset
	if (mxGetNumberOfElements(gradfilteroffset)!=1) { 
		mexErrMsgTxt("gradfilteroffset should be scalar.");
	}
	int gradfilteroffseti = static_cast<int>( mxGetScalar(gradfilteroffset) ) ; // should be smaller (in magnitude) than gradfilterlength
#endif
	// Calculate temporary memory required size:  (follows tmp requests in TV_core)
    mwSize ntempimg = 0;
	mwSize ntempInt = ndims;		//cumSpatialSize 
	mwSize ntempTVt = 0;		// gradcache

	double * pH;
	double * px;
	if (hessmul ) {
		// read hessian multiplication specific inputs:
		const mxArray * H = prhs[4];
		const mxArray * x  = prhs[5];
		mwSize nvox = 1;
		for (int dim = 1 ; dim< ndims; dim++ ) {
			nvox *= size[dim];
		}
		//H:
		if (ndims != mxGetNumberOfDimensions( H ) || nvox!=mxGetNumberOfElements( H ) || nvox*size[0] != mxGetNumberOfElements( x ) ) {
			mexErrMsgTxt("H should have same number of dimensions as img and numel(H) should be number of voxels should, and numel(x) should be numel(img)");
		}
		if ((imgclass!=mxGetClassID(H)) || (imgclass!= mxGetClassID(x)) ) {
			mexErrMsgTxt("H and x should be of same type as image.");
		}
		if (mxIsComplex(H)|| mxIsComplex(x)) {
			mexErrMsgTxt("H and x should not be complex.");
		}
		pH = mxGetPr(H);
		px = mxGetPr(x);


		// proceed with temporary memory calculations"
		ntempimg += size[0]*2*(1+ 2*(ndims-1)) ; //G_gs & G_dg & G_gx
		ntempInt += 2*2*ndims;				// G_steps & G_dim
		ntempInt += ndims;					// pos
	} else {
		ntempTVt += ndims;		// gradcache
		mwSize cmsz=1;
		for (int dim = 1; dim<ndims; dim++ ) {
			ntempTVt += cmsz;			// gradcache[dim]
			cmsz *= size[dim];
		};
		ntempimg += size[0];			// g
		if (computeGrad) {
			if (computeExplicitHess) {
				ntempimg += size[0]*2*2*(ndims-1) ; //G_g & G_dg
				ntempInt += 2*ndims ;
				if (isvecw ) {
					ntempTVt += size[0];		// for Hg, hessian
					ntempInt += size[0];		// index for hessian
				} else {
					ntempTVt += size[0]*size[0];// hessian of distance
					ntempInt += size[0]*size[0];// index of hessian of distance.
				}
			} else {
				ntempimg += size[0]*2;			//G_g & G_dg
			}
			ntempInt += 2*ndims;				// G_dim
		}
		ntempInt += ndims*2;	// pos & curGradCachePos
	}
	ntempInt += ndims; // stepDim;

	mwSize n_tempspace_bytes = mxGetElementSize(img) * ntempimg + sizeof(mwSize) * ntempInt + sizeof(double) * ntempTVt;
#ifdef ADDGUARDS 
	n_tempspace_bytes +=1000; // to account for extra bytes used by guards.
#endif
    char * tempspaceptr = (char *) mxMalloc( n_tempspace_bytes );
	tempspace tmp = tempspace(tempspaceptr, n_tempspace_bytes);

   
    ptrdiff_t * stepDim = tmp.get<ptrdiff_t>( ndims ); 
    stepDim[0] = 1;
    for (int dim = 1; dim<ndims; dim++) {
        stepDim[dim] = stepDim[dim-1] * size[dim-1];
    }
    int sz0 = static_cast<int>( size[0] ) ; // vectorlength (of image) should be rather short, certainly less than maximum integer (could check for that though).


	// create output:
    double sqrtHessModifScaleFact = 1;
	if (hessmul ) {
		// perform hessian multiplication.
		const mxArray * x  = prhs[5];
		plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions( x ), mxGetDimensions(x), imgclass, mxREAL); // create Hx
		double * Hx = mxGetPr( plhs[0] );

		if (imgclass==mxDOUBLE_CLASS) {
			if (isvecw) {
				dist_w_vecw< double, double* > dist( scale, sz0 );
				if (voxspacing==NULL) {
					TV_hessianMul( px, Hx, source, pH, size, stepDim,             ndims, sqrtHessModifScaleFact, dist GRADFILTERARGSs, tmp);
				} else {
					TV_hessianMul( px, Hx, source, pH, size, stepDim, voxspacing, ndims, sqrtHessModifScaleFact, dist GRADFILTERARGSs, tmp);
				}
			} else {
				dist_w_matw< double, double * > dist( scale, sz0 );
				if (voxspacing==NULL) {
					TV_hessianMul( px, Hx, source, pH, size, stepDim,	         ndims, sqrtHessModifScaleFact, dist GRADFILTERARGSs, tmp);
				} else {
					TV_hessianMul( px, Hx, source, pH, size, stepDim, voxspacing, ndims, sqrtHessModifScaleFact, dist GRADFILTERARGSs, tmp);
				}
			}
		} /*else if (imgclass==mxSINGLE_CLASS) {
			if (isvecw) {
				dist_w_vecw< double, float* > dist( (float *) scale, sz0 );
				if (voxspacing==NULL) {
					TV = TV_core( (float*) source, (float*) grad, size, stepDim, offset, ndims, dist, hess_pr, hess_ir, hess_jc, numHessEl, tmp);
				} else {
					TV = TV_core( (float*) source, (float*) grad, size, stepDim, offset, (float *) voxspacing, ndims, dist, hess_pr, hess_ir, hess_jc, numHessEl, tmp);
				}
			} else {
				dist_w_matw< double, float* > dist( (float *) scale, sz0 );
				if (voxspacing==NULL) {
					TV = TV_core( (float*) source, (float*) grad, size, stepDim, offset, ndims, dist, hess_pr, hess_ir, hess_jc, numHessEl, tmp);
				} else {
					TV = TV_core( (float*) source, (float*) grad, size, stepDim, offset, (float *) voxspacing, ndims, dist, hess_pr, hess_ir, hess_jc, numHessEl, tmp);
				}
			}
		}*/ else 
			mexErrMsgTxt("Unsupported image type.");


	} else {
		double * grad = NULL;
		if (computeGrad) {
			plhs[1] = mxCreateNumericArray(ndims, size, imgclass, mxREAL); // create grad
			grad = mxGetPr(plhs[1]);
		}
		mwSize numHessEl = 0;
		mwSize numel = 0;
		double * hess_pr = NULL;
		mwIndex * hess_ir = NULL;
		mwIndex * hess_jc = NULL;
		if (computeHess) {
			if (computeExplicitHess) {
				mwSize cumsz = 1;
				for (int dim = 1; dim<ndims; dim++) {
					cumsz  *= size[dim];
				}
				numel = cumsz * size[0];
				numHessEl = cumsz ;
				for (int dim = 1; dim<ndims; dim++) {
					numHessEl += 4 * cumsz/size[dim] * (size[dim]-1);
				}
				numHessEl += cumsz * (2*ndims)*(2*ndims);
				numHessEl *= size[0] *size[0];
				if (returnSparse) {
					// return as sparse
					plhs[2] = mxCreateSparse(numel, numel, numHessEl, mxREAL); // create sparse hessian
					hess_pr = mxGetPr(plhs[2]);
					hess_ir = mxGetIr(plhs[2]);
					hess_jc = mxGetJc(plhs[2]);
				} else {
					// return components of sparse.
					plhs[2] = mxCreateDoubleMatrix(numHessEl, 1, mxREAL); // create sparse hessian
					hess_pr = mxGetPr(plhs[2]);
					plhs[3] = mxCreateNumericMatrix(numHessEl, 1, mxINT64_CLASS, mxREAL); // create sparse hessian
					hess_ir = (mwIndex *) mxGetPr(plhs[3]);
					plhs[4] = mxCreateNumericMatrix(numHessEl, 1, mxINT64_CLASS, mxREAL); // create sparse hessian
					hess_jc = (mwIndex *) mxGetPr(plhs[4]);
				}
			} else {
				tempspace tmpcopy = tmp;
				mwSize * sizeH = tmpcopy.get<mwSize>( ndims ); 
				sizeH[0] = 1;
				for (int dim = 1 ; dim < ndims; dim++ ) {
					sizeH[dim] = size[dim];
				}
				plhs[2] = mxCreateNumericArray(ndims, sizeH, imgclass , mxREAL); 
				hess_pr = mxGetPr(plhs[2]);
			}
		}
		double TV =0;

		// Calculate TV:
// template <typename TVtype, typename imgtype, typename distMeasureType, bool skipborder=false >
// inline TVtype TV_core( const imgtype * __restrict source, imgtype * __restrict grad, 
// 					   const mwSize * __restrict size, const ptrdiff_t * __restrict stepSrc,
//                        TVtype offset  VOXELSPACINGARGS , int ndims ,
// 					   distMeasureType & dist  HESSIANARGS GRADFILTERARGS ,  tempspace tmp)
        
        
		if (imgclass==mxDOUBLE_CLASS) {
            typedef double TVtype;
            typedef double imgtype;
			if (isvecw) {
                typedef dist_w_vecw< double, double* > distMeasureType;
				distMeasureType dist( scale, sz0 );
                if (voxspacing==NULL) {
                    TV = TV_core<TVtype, imgtype, distMeasureType, SKIPBORDER>( source, grad, size, stepDim, offset,             ndims, dist HESSIANARGSs GRADFILTERARGSs, tmp);
                } else {
                    TV = TV_core<TVtype, imgtype, distMeasureType, SKIPBORDER>( source, grad, size, stepDim, offset, voxspacing, ndims, dist HESSIANARGSs GRADFILTERARGSs, tmp);
                }
			} else {
                typedef dist_w_matw< double, double * > distMeasureType;
				distMeasureType dist( scale, sz0 );
                if (voxspacing==NULL) {
                    TV = TV_core<TVtype, imgtype, distMeasureType, SKIPBORDER>( source, grad, size, stepDim, offset,             ndims, dist HESSIANARGSs GRADFILTERARGSs, tmp);
                } else {
                    TV = TV_core<TVtype, imgtype, distMeasureType, SKIPBORDER>( source, grad, size, stepDim, offset, voxspacing, ndims, dist HESSIANARGSs GRADFILTERARGSs, tmp);
                }
			}
		} /*else if (imgclass==mxSINGLE_CLASS) {
			if (isvecw) {
				dist_w_vecw< double, float* > dist( (float *) scale, sz0 );
				if (voxspacing==NULL) {
					TV = TV_core( (float*) source, (float*) grad, size, stepDim, offset, ndims, dist, hess_pr, hess_ir, hess_jc, numHessEl, tmp);
				} else {
					TV = TV_core( (float*) source, (float*) grad, size, stepDim, offset, (float *) voxspacing, ndims, dist, hess_pr, hess_ir, hess_jc, numHessEl, tmp);
				}
			} else {
				dist_w_matw< double, float* > dist( (float *) scale, sz0 );
				if (voxspacing==NULL) {
					TV = TV_core( (float*) source, (float*) grad, size, stepDim, offset, ndims, dist, hess_pr, hess_ir, hess_jc, numHessEl, tmp);
				} else {
					TV = TV_core( (float*) source, (float*) grad, size, stepDim, offset, (float *) voxspacing, ndims, dist, hess_pr, hess_ir, hess_jc, numHessEl, tmp);
				}
			}
		}*/ else 
			mexErrMsgTxt("Unsupported image type.");

		if ( returnSparse && computeExplicitHess ) {
			mwIndex nzmax = hess_jc[numel];
			mxSetNzmax(plhs[2], nzmax); 
			mxSetPr(plhs[2], (double*) mxRealloc(hess_pr, nzmax*sizeof(double)));
			mxSetIr(plhs[2], (mwIndex*) mxRealloc(hess_ir, nzmax*sizeof(mwIndex)));
		}
	    plhs[0] = mxCreateDoubleScalar(TV);
	}
#ifdef ADDGUARDS 
	tmp.checkGuard( stepDim , ndims ); 
#endif
	mxFree(tempspaceptr);

}
#endif

