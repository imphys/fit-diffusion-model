#ifndef EMM_WRAP_CPP // only include once.
#define EMM_WRAP_CPP

#include "emmintrin.h"
#include "xmmintrin.h"
//#define INCLUDE_SSE4
#define INCLUDE_SSE3

#ifdef INCLUDE_SSE4
#include <smmintrin.h>
#endif

#ifdef _MSC_VER
	// Compiling with MS Visual studio compiler:
	#include <intrin.h>
	// Convert MS macros/definitions:
	typedef __int32 int32_t;
	typedef __int64 int64_t; 
	#ifdef _M_X64 
		#define __x86_64__ _M_X64
	#endif
	#define ALIGN16 __declspec(align(16))
#else // add check for compiling with gcc.
	// Compiling with gcc
	#include <pmmintrin.h>
	#include <stdint.h>
	#define ALIGN16 __attribute__ ((aligned (16)))
	#define _MM_FROUND_TO_NEAREST_INT 0x0
#endif

#include <stdarg.h>
#define BAD mexErrMsgTxt("Currently unsupported");
#define BADARG mexErrMsgTxt("Wrong argument supplied.");

#define ALIGN_VEC ALIGN16

/* This file provides a templated wrapper around the xmm intrinsic functions.
 * Currently mainly single/double precision operations are defined, and some
 * int32 and int64 operations.
 * Currently (mainly) vectors of length 2 and 4 are supported.
 * Call vec<datatype, ignored>.naturalLen() to get the natural(/minimum) vector length for datatype.
 *
 * All standard operators are defined ( +, -, *, /, +=, -=, *=, /= )
 * Both with another vec (with the same type & length) as well as with scalars of the correct type 
 *
 * Type specification:  vec< data_type, vector_length >
 * Create new        :  vec< data_type, vector_length >( pointer_to_datatype [, stride_in_#elements] )
 *                      vec< data_type, vector_length >( variable_of_datatype )
 * Type conversion (double <-> single) is available.
 *
 * A few special member routines are provided as well. 
 * - A (very unsafe) scale routine ( v = v * 2^i, with vec<integer_type, len> i)
 *
 * And finally, some important functions that operate on vector objects:
 * - min, max
 * - round  : round each element to nearest integer
 * - accum  : accumulate values in vector and return scalar of datatype
 * - sqr    : return the element wise square of a vector.
 *
 * Created by Dirk Poot, Erasmus MC, 
 * Last modified 20-9-2010
 */

template <typename T, int vlen>  union vec
{
	typedef T value_type;
	enum {length = vlen};

    //Internal data fields:
	__m128i xmmi[sizeof(T)*vlen/16];
	__m128 xmm[sizeof(T)*vlen/16]; 
	__m128d xmmd[sizeof(T)*vlen/16]; 

    const static int naturalLen() {
        return (16/sizeof(T));
    }

    // Main external constructors:
	explicit vec (const T *v);
	template <typename int_t> explicit vec (const T *v, int_t stride);
//	vec (const T *v, int stride);
	explicit vec (const T v);
	static vec zero();
	explicit vec () {}; //  default constructor, no initialization.

    
    // Internal constructors:
	explicit vec(__m128 in1) {
		if (vlen * sizeof(T) !=1 *16) 
			BAD;
		xmm[0] = in1;
	}
	explicit vec(__m128 in1, __m128 in2) {
		if (vlen * sizeof(T) !=2 *16) 
			BAD;
		xmm[0] = in1;
		xmm[1] = in2;
	}
	explicit vec(__m128i in1) {
		if (vlen * sizeof(T) !=1 *16) 
			BAD;
		xmmi[0] = in1;
	}
	explicit vec(__m128i in1, __m128i in2) {
		if (vlen * sizeof(T) !=2 *16) 
			BAD;
		xmmi[0] = in1;
		xmmi[1] = in2;
	}
	explicit vec(__m128d in1) {
		if (vlen * sizeof(T) !=1 *16) 
			BAD;
		xmmd[0] =  in1;
	}
	explicit vec(__m128d in1, __m128d in2) {
		if (vlen * sizeof(T) !=2 *16) 
			BAD;
		xmmd[0] =  in1;
		xmmd[1] =  in2;
	}
    explicit vec(__m128 in[], const int arraylen) {
		if (vlen * sizeof(T) !=arraylen*16) 
			BADARG;
		xmm[0] = in[0];
		if (vlen * sizeof(T)>1* 16) {
		xmm[1] = in[1];
		if (vlen * sizeof(T)>2* 32) {
		xmm[2] = in[2];
		if (vlen * sizeof(T)>3* 32) {
		xmm[3] = in[3];
		if (vlen * sizeof(T)>4* 32) {
			BAD;
		}}}}
	}

	template < typename ToType > inline vec< ToType, (vlen*(sizeof (T))/(sizeof (ToType))) > reinterpret() 
	{
		return vec<ToType, vlen*(sizeof (T))/(sizeof (ToType)) >( xmm , vlen * (sizeof (T)) /16);
	}

    // Type conversion routine:
	template <typename From, int vlenFrom> explicit vec(const vec<From, vlenFrom> &v);
	
    // Store to memory routines:
	inline void store( T *v);  // store unaligned.
    typedef ALIGN16 T Ta;
	inline void storea( Ta * v); // store aligned.
	template <typename int_t> inline void store( T *v, int_t stride);
    
    
	template <typename From, int vlenFrom> inline void scale_unsafe(const vec<From, vlenFrom> &v);

	inline vec operator* (const vec<T, vlen> &v) const;
	inline vec operator+ (const vec<T, vlen> &v) const;
	inline vec operator- (const vec<T, vlen> &v) const;
	inline vec operator/ (const vec<T, vlen> &v) const;
	inline void operator*= (const vec<T, vlen> &v);
	inline void operator+= (const vec<T, vlen> &v);
	inline void operator-= (const vec<T, vlen> &v);
	inline void operator/= (const vec<T, vlen> &v);

	// scalar versions:
	inline vec operator* (const T &v) const;
	inline vec operator+ (const T &v) const;
	inline vec operator- (const T &v) const;
	inline vec operator/ (const T &v) const;
	inline void operator*= (const T &v);
	inline void operator+= (const T &v);
	inline void operator-= (const T &v);
	inline void operator/= (const T &v);

    // comparison to bitmask:
    inline vec operator> (const vec &v) const;
    inline vec operator>= (const vec &v) const;
    inline vec operator== (const vec &v) const;
    inline vec operator<= (const vec &v) const;
    inline vec operator< (const vec &v) const;
    inline vec operator!= (const vec &v) const;
    
	// bit wise operators:
    inline vec operator>> (const int shiftcount) const; 
    inline vec operator<< (const int shiftcount) const; 
    inline vec operator| (const vec<T, vlen> &v) const; 
    inline vec operator& (const vec<T, vlen> &v) const; 


/*	// other operations:
	static inline vec max(const vec<T, vlen> &a, const vec<T, vlen> &b);
	static inline vec min(const vec<T, vlen> &a, const vec<T, vlen> &b);
	static inline T accum(const vec<T, vlen> &v);*/
	inline vec rep(int idx) const; // Replicate vec[idx]; use only for compile time constant idx.
	inline vec insert( const vec &v, const vec &mask );
};

// Below the functions are defined in the following way. 
// Largest sections: 
//  - Load functions
//  - Store functions
//  - Operators
//  - Other functions (min/max ...)
//
//  For those applicable, vec - scalar  versions are stratisfied within type stratification.
//  Each of these is stratisfied on type (double/float ...)



// Load functions:
/* Unfortunately MSVC does (currently) not optimize the generic load function:
template <typename T, int vlen>  vec<T, vlen>::vec(const T *v) { 
  if( sizeof(T)==2) {
      for (int i=0; i<vlen/8; i++) {
        xmm[i] = _mm_loadu_ps((float *) (v + i*8)); 
      }
  } else if( sizeof(T)==4) {
      for (int i=0; i<vlen/4; i++) {
        xmm[i] = _mm_loadu_ps((float *) (v + i*4)); 
      }
  } else if (sizeof(T)==8) {
      for (int i=0; i<vlen/2; i++) {
        xmmd[i] = _mm_loadu_pd((double *) (v + i*2)); 
      }
  } else {BAD};
}*/

// load vector, (double/float, length 2/4/8)
template<>  vec<double, 2>::vec(const double *v) { 
	if (__alignof(v)>=16) {
		xmmd[0] = _mm_load_pd(v ); 
	} else { 
		xmmd[0] = _mm_loadu_pd(v ); 
	}
}

template <> template<typename int_t>  vec<double, 2>::vec(const double *v, int_t stride) { 
      xmmd[0] = _mm_load_sd(v ); 
      xmmd[0] = _mm_loadh_pd(xmmd[0], v+stride ); 
}
template<>  vec<double, 4>::vec(const double *v) { 
	if (__alignof(v)>=16) {
      xmmd[0] = _mm_load_pd(v ); 
      xmmd[1] = _mm_load_pd(v + 2); 
	} else {
      xmmd[0] = _mm_loadu_pd(v ); 
      xmmd[1] = _mm_loadu_pd(v + 2); 
	}
}
template<> template<typename int_t> vec<double, 4>::vec(const double *v, int_t stride) { 
      xmmd[0] = _mm_load_sd(v ); 
      xmmd[0] = _mm_loadh_pd(xmmd[0], v+stride ); 
      xmmd[1] = _mm_load_sd(v + 2*stride); 
      xmmd[1] = _mm_loadh_pd(xmmd[1], v+3*stride ); 
}

template<>  vec<float, 4>::vec(const float *v) { 
      xmm[0] = _mm_loadu_ps(v ); 
}
template<> template<typename int_t>  vec<float, 4>::vec(const float *v, int_t stride) { 
      xmm[0] = _mm_set_ps( v[0], v[stride], v[2*stride], v[3*stride]); 
}
template<>  vec<float, 8>::vec(const float *v) { 
      xmm[0] = _mm_loadu_ps(v ); 
      xmm[1] = _mm_loadu_ps(v + 4); 
}
template<>  vec<int64_t, 2>::vec(const int64_t *v) {
      xmm[0] = _mm_loadu_ps((const float *) v );
}

template<>  vec<double, 2>::vec(const double v) { 
      xmmd[0] = _mm_load_sd(&v ); 
	  xmmd[0] = _mm_unpacklo_pd(xmmd[0], xmmd[0]);
}
template<>  vec<double, 4>::vec(const double v) { 
      xmmd[0] = _mm_load_sd(&v ); 
	  xmmd[0] = _mm_unpacklo_pd(xmmd[0], xmmd[0]);
      xmmd[1] = xmmd[0]; 
}
template<>  vec<float, 4>::vec(const float v) { 
      xmm[0] = _mm_load_ss(&v ); 
	  xmm[0] = _mm_shuffle_ps(xmm[0],xmm[0],_MM_SHUFFLE(0, 0, 0, 0));
	  //xmm[0] = _mm_unpacklo_ps(xmm[0],xmm[0]); 
	  //xmm[0] = _mm_unpacklo_ps(xmm[0],xmm[0]); 
}
template<>  vec<float, 8>::vec(const float v) { 
      xmm[0] = _mm_load_ss(&v ); 
	  xmm[0] = _mm_shuffle_ps(xmm[0],xmm[0],_MM_SHUFFLE(0, 0, 0, 0));
	  xmm[1] = xmm[0];
}

//create as zero vector:

//template <> static vec<double, 4> vec<double, 4>::zero () { 
template <> vec<double, 4> vec<double, 4>::zero () { 
	return vec<double,4>(_mm_setzero_pd(), _mm_setzero_pd()); 
}
//template <> static vec<float, 4> vec<float, 4>::zero () { 
template <> vec<float, 4> vec<float, 4>::zero () { 
	return vec<float,4>(_mm_setzero_ps()); 
}

// Store functions
/*Unfortunately MSVC does (currently) not optimize the generic store function:
template <typename T, int vlen> void vec<T, vlen>::store(T *v) { 
	if( sizeof(T)==2) {
	  for (int i=0; i<vlen/8; i++) {
		  _mm_storeu_ps((float *) (v + i*8), xmm[i]);
	  }
	} else if( sizeof(T)==4) {
	  for (int i=0; i<vlen/4; i++) {
	   _mm_storeu_ps((float *) (v + i*4), xmm[i]);
	  }
	} else if (sizeof(T)==8) {
	  for (int i=0; i<vlen/2; i++) {
	   _mm_storeu_pd((double *) (v + i*2), xmmd[i]);
	  }
	} else {BAD};
}*/

template <> void vec<double, 2>::store(double *v) { 
	_mm_storeu_pd( (v    ), xmmd[0]);
}
template <> void vec<double, 2>::storea(double *v) {
	_mm_store_pd( (v    ), xmmd[0]);
}
template<> template<typename int_t>  void vec<double, 2>::store(double *v, int_t stride) { 
	_mm_store_sd(  (v    ), xmmd[0]);
	_mm_storeh_pd( (v +   stride), xmmd[0]);
}
template <> void vec<double, 4>::store(double *v) { 
	_mm_storeu_pd( (v    ), xmmd[0]);
	_mm_storeu_pd( (v + 2), xmmd[1]);
}
template <> void vec<double, 4>::storea(double *v) {
	_mm_store_pd( (v    ), xmmd[0]);
	_mm_store_pd( (v + 2), xmmd[1]);
}
template<> template<typename int_t> void vec<double, 4>::store(double *v, int_t stride) { 
	_mm_store_sd(  (v    ), xmmd[0]);
	_mm_storeh_pd( (v +   stride), xmmd[0]);
	_mm_store_sd(  (v + 2*stride), xmmd[1]);
	_mm_storeh_pd( (v + 3*stride), xmmd[1]);
}

template <> void vec<float, 4>::store(float *v) { 
	_mm_storeu_ps( (v    ), xmm[0]);
}
template <> void vec<float, 4>::storea(float *v) {
	_mm_store_ps( (v    ), xmm[0]);
}
template<> template<typename int_t>  void vec<float, 4>::store(float *v, int_t stride) { 
	_mm_store_ss( (v           ), xmm[0]);
    xmm[0] = _mm_shuffle_ps(xmm[0],xmm[0],_MM_SHUFFLE(0, 3, 2, 1));    
	_mm_store_ss( (v +   stride), xmm[0]);
    xmm[0] = _mm_shuffle_ps(xmm[0],xmm[0],_MM_SHUFFLE(0, 3, 2, 1));    
	_mm_store_ss( (v + 2*stride), xmm[0]);
    xmm[0] = _mm_shuffle_ps(xmm[0],xmm[0],_MM_SHUFFLE(0, 3, 2, 1));    
	_mm_store_ss( (v + 3*stride), xmm[0]);
    xmm[0] = _mm_shuffle_ps(xmm[0],xmm[0],_MM_SHUFFLE(0, 3, 2, 1)); // leave in original state
}
template <> void vec<float, 8>::store(float *v) { 
	_mm_storeu_ps( (v    ), xmm[0]);
	_mm_storeu_ps( (v + 4), xmm[1]);
}
template <> void vec<int64_t, 2>::store(int64_t *v) { 
	_mm_storeu_pd( (double *) (v    ), xmmd[0]);
}
template <> void vec<int64_t, 2>::storea(int64_t *v) {
	_mm_store_si128( (__m128i *) (v    ), xmmi[0]);
}

// Type conversion constructors:
/*template<> template<> vec<double, 4>::vec(const vec<int32_t, 4> &v) { 
	xmmd[0] = _mm_cvtepi32_pd(v.xmmi[0]); 
	xmmi[1] = _mm_unpackhi_epi64(v.xmmi[0],v.xmmi[0]);
    xmmd[1] = _mm_cvtepi32_pd(xmmi[1]); 
}
template<> template<> inline vec<int32_t, 4>::vec(const vec<double, 4> &v) { 
	xmmi[0] =     _mm_cvtpd_epi32(v.xmmd[0]); 
	__m128i tmp = _mm_cvtpd_epi32(v.xmmd[1]); 
	xmmi[0] = _mm_unpacklo_epi64(xmmi[0], tmp);
}*/
template<> template<> inline vec<int32_t, 4>::vec(const vec<float, 4> &v) { 
	xmmi[0] = _mm_cvtps_epi32(v.xmm[0]); 
}
#ifdef INCLUDE_SSE4
template<> template<> inline vec<int64_t, 2>::vec(const vec<double, 2> &v) { 
	xmmi[0] = _mm_cvtpd_epi32(v.xmmd[0]); 
	xmmi[0] = _mm_cvtepi32_epi64(xmmi[0]);
}
#endif
template<> template<> inline vec<double, 2>::vec(const vec<int64_t, 2> &v) { 
	//xmmi[0] = _mm_unpacklo_epi32(v.xmmi[0],v.xmmi[0]);
	xmmi[0] = _mm_shuffle_epi32(v.xmmi[0],_MM_SHUFFLE(2, 0, 2, 0));
	xmmd[0] = _mm_cvtepi32_pd(xmmi[0]); 
}

template<> template<> inline vec<float, 4>::vec(const vec<double, 4> &v) { 
	xmm[0] = _mm_cvtpd_ps(v.xmmd[0]);
    __m128 tmp = _mm_cvtpd_ps(v.xmmd[1]);
	xmm[0] = _mm_shuffle_ps(xmm[0], tmp, _MM_SHUFFLE(1, 0, 1, 0)); 
}

// Other members:
template <> template<> void vec<double, 2>::scale_unsafe(const vec<int64_t, 2> &v) { 
    // Computes self = self * 2^v
	// No INF or NAN checks, sign (and exponent) of value corrupted when exponent overflows.
	__m128i tmp = _mm_slli_epi64(v.xmmi[0],52);
	xmmi[0] = _mm_add_epi64(tmp, xmmi[0]);
}
template <> template<> void vec<float, 4>::scale_unsafe(const vec<int32_t, 4> &v) { 
    // Computes self = self * 2^v
	// No INF or NAN checks, sign (and exponent) of value corrupted when exponent overflows.
	__m128i tmp = _mm_slli_epi32(v.xmmi[0],23);
	xmmi[0] = _mm_add_epi32(tmp, xmmi[0]);
}

// Operators (Generic):
/*template <typename T, int vlen>  vec<T,vlen> vec<T, vlen>::operator* (const vec<T,vlen> &v) const	{  
	BAD;return v;
}
template <typename T, int vlen>  vec<T,vlen> vec<T, vlen>::operator+ (const vec<T, vlen> &v) const {
	BAD;return v;
}
template <typename T, int vlen>  vec<T,vlen> vec<T, vlen>::operator- (const vec<T, vlen> &v) const {
	BAD;return v;
}
template <typename T, int vlen>  vec<T,vlen> vec<T, vlen>::operator/ (const vec<T, vlen> &v) const {
	BAD;return v;
}*/

// Operators, Specialized versions (double precision, length 4):
template <> vec<double, 4> vec<double, 4>::operator* (const vec<double,4> &v) const	{ 
	return vec<double,4>(_mm_mul_pd(xmmd[0], v.xmmd[0]), _mm_mul_pd(xmmd[1], v.xmmd[1])); 
}
template <> vec<double, 4> vec<double, 4>::operator+ (const vec<double,4> &v) const	{ 
	return vec<double,4>(_mm_add_pd(xmmd[0], v.xmmd[0]), _mm_add_pd(xmmd[1], v.xmmd[1])); 
}
template <> vec<double, 4> vec<double, 4>::operator- (const vec<double,4> &v) const	{ 
	return vec<double,4>(_mm_sub_pd(xmmd[0], v.xmmd[0]), _mm_sub_pd(xmmd[1], v.xmmd[1])); 
}
template <> vec<double, 4> vec<double, 4>::operator/ (const vec<double,4> &v) const	{ 
	return vec<double,4>(_mm_div_pd(xmmd[0], v.xmmd[0]), _mm_div_pd(xmmd[1], v.xmmd[1])); 
}
template <> inline void vec<double, 4>::operator*= (const vec<double, 4> &v) { 
	xmmd[0] = _mm_mul_pd(xmmd[0], v.xmmd[0]); 
	xmmd[1] = _mm_mul_pd(xmmd[1], v.xmmd[1]); 
}
template <> inline void vec<double, 4>::operator+= (const vec<double, 4> &v) { 
	xmmd[0] = _mm_add_pd(xmmd[0], v.xmmd[0]); 
	xmmd[1] = _mm_add_pd(xmmd[1], v.xmmd[1]); 
}
template <> inline void vec<double, 4>::operator-= (const vec<double, 4> &v) { 
	xmmd[0] = _mm_sub_pd(xmmd[0], v.xmmd[0]); 
	xmmd[1] = _mm_sub_pd(xmmd[1], v.xmmd[1]); 
}
template <> inline void vec<double, 4>::operator/= (const vec<double, 4> &v) {
	xmmd[0] = _mm_div_pd(xmmd[0], v.xmmd[0]); 
	xmmd[1] = _mm_div_pd(xmmd[1], v.xmmd[1]); 
}

// Operators, Specialized versions (double precision, length 2):
template <> vec<double, 2> vec<double, 2>::operator* (const vec<double,2> &v) const	{ 
	return vec<double,2>(_mm_mul_pd(xmmd[0], v.xmmd[0])); 
}
template <> vec<double, 2> vec<double, 2>::operator+ (const vec<double,2> &v) const	{ 
	return vec<double,2>(_mm_add_pd(xmmd[0], v.xmmd[0])); 
}
template <> vec<double, 2> vec<double, 2>::operator- (const vec<double,2> &v) const	{ 
	return vec<double,2>(_mm_sub_pd(xmmd[0], v.xmmd[0])); 
}
template <> vec<double, 2> vec<double, 2>::operator/ (const vec<double,2> &v) const	{ 
	return vec<double,2>(_mm_div_pd(xmmd[0], v.xmmd[0])); 
}
template <> inline void vec<double, 2>::operator*= (const vec<double, 2> &v) { 
	xmmd[0] = _mm_mul_pd(xmmd[0], v.xmmd[0]); 
}
template <> inline void vec<double, 2>::operator+= (const vec<double, 2> &v) { 
	xmmd[0] = _mm_add_pd(xmmd[0], v.xmmd[0]); 
}
template <> inline void vec<double, 2>::operator-= (const vec<double, 2> &v) { 
	xmmd[0] = _mm_sub_pd(xmmd[0], v.xmmd[0]); 
}
template <> inline void vec<double, 2>::operator/= (const vec<double, 2> &v) {
	xmmd[0] = _mm_div_pd(xmmd[0], v.xmmd[0]); 
}
template <> vec<double, 2> vec<double, 2>::operator> (const vec<double,2> &v) const	{ 
	return vec<double,2>(_mm_cmpgt_pd(xmmd[0], v.xmmd[0])); 
}
template <> vec<double, 2> vec<double, 2>::operator>= (const vec<double,2> &v) const	{ 
	return vec<double,2>(_mm_cmpge_pd(xmmd[0], v.xmmd[0])); 
}
template <> vec<double, 2> vec<double, 2>::operator== (const vec<double,2> &v) const	{ 
	return vec<double,2>(_mm_cmpeq_pd(xmmd[0], v.xmmd[0])); 
}
template <> vec<double, 2> vec<double, 2>::operator<= (const vec<double,2> &v) const	{ 
	return vec<double,2>(_mm_cmple_pd(xmmd[0], v.xmmd[0])); 
}
template <> vec<double, 2> vec<double, 2>::operator< (const vec<double,2> &v) const	{ 
	return vec<double,2>(_mm_cmplt_pd(xmmd[0], v.xmmd[0])); 
}
template <> vec<double, 2> vec<double, 2>::operator!= (const vec<double,2> &v) const	{ 
	return vec<double,2>(_mm_cmpneq_pd(xmmd[0], v.xmmd[0])); 
}

// Operators, Specialized versions (single precision, length 8):
template <> vec<float, 8> vec<float, 8>::operator* (const vec<float,8> &v) const	{ 
	return vec<float,8>(_mm_mul_ps(xmm[0], v.xmm[0]), _mm_mul_ps(xmm[1], v.xmm[1])); 
}
template <> vec<float, 8> vec<float, 8>::operator+ (const vec<float,8> &v) const	{ 
	return vec<float,8>(_mm_add_ps(xmm[0], v.xmm[0]), _mm_add_ps(xmm[1], v.xmm[1])); 
}
template <> vec<float, 8> vec<float, 8>::operator- (const vec<float,8> &v) const	{ 
	return vec<float,8>(_mm_sub_ps(xmm[0], v.xmm[0]), _mm_sub_ps(xmm[1], v.xmm[1])); 
}
template <> vec<float, 8> vec<float, 8>::operator/ (const vec<float,8> &v) const	{ 
	return vec<float,8>(_mm_div_ps(xmm[0], v.xmm[0]), _mm_div_ps(xmm[1], v.xmm[1])); 
}
template <> inline void vec<float, 8>::operator*= (const vec<float, 8> &v) { 
	xmm[0] = _mm_mul_ps(xmm[0], v.xmm[0]); 
	xmm[1] = _mm_mul_ps(xmm[1], v.xmm[1]); 
}
template <> inline void vec<float, 8>::operator+= (const vec<float, 8> &v) { 
	xmm[0] = _mm_add_ps(xmm[0], v.xmm[0]); 
	xmm[1] = _mm_add_ps(xmm[1], v.xmm[1]); 
}
template <> inline void vec<float, 8>::operator-= (const vec<float, 8> &v) { 
	xmm[0] = _mm_sub_ps(xmm[0], v.xmm[0]); 
	xmm[1] = _mm_sub_ps(xmm[1], v.xmm[1]); 
}
template <> inline void vec<float, 8>::operator/= (const vec<float, 8> &v) {
	xmm[0] = _mm_div_ps(xmm[0], v.xmm[0]); 
	xmm[1] = _mm_div_ps(xmm[1], v.xmm[1]); 
}

// Operators, Specialized versions (single precision, length 4):
template <> vec<float, 4> vec<float, 4>::operator* (const vec<float,4> &v) const	{ 
	return vec<float,4>(_mm_mul_ps(xmm[0], v.xmm[0])); 
}
template <> vec<float, 4> vec<float, 4>::operator+ (const vec<float,4> &v) const	{ 
	return vec<float,4>(_mm_add_ps(xmm[0], v.xmm[0])); 
}
template <> vec<float, 4> vec<float, 4>::operator- (const vec<float,4> &v) const	{ 
	return vec<float,4>(_mm_sub_ps(xmm[0], v.xmm[0])); 
}
template <> vec<float, 4> vec<float, 4>::operator/ (const vec<float,4> &v) const	{ 
	return vec<float,4>(_mm_div_ps(xmm[0], v.xmm[0])); 
}
template <> inline void vec<float, 4>::operator*= (const vec<float, 4> &v) { 
	xmm[0] = _mm_mul_ps(xmm[0], v.xmm[0]); 
}
template <> inline void vec<float, 4>::operator+= (const vec<float, 4> &v) { 
	xmm[0] = _mm_add_ps(xmm[0], v.xmm[0]); 
}
template <> inline void vec<float, 4>::operator-= (const vec<float, 4> &v) { 
	xmm[0] = _mm_sub_ps(xmm[0], v.xmm[0]); 
}
template <> inline void vec<float, 4>::operator/= (const vec<float, 4> &v) {
	xmm[0] = _mm_div_ps(xmm[0], v.xmm[0]); 
}


//  Operators, scalar versions (double, length 2):
template <> vec<double, 2> vec<double, 2>::operator* (const double &v) const	{ 
	__m128d tmp = _mm_load_sd(&v);
	tmp = _mm_unpacklo_pd(tmp, tmp);
	return vec<double,2>(_mm_mul_pd(xmmd[0], tmp)); 
}
template <> vec<double, 2> vec<double, 2>::operator+ (const double &v) const	{ 
	__m128d tmp = _mm_load_sd(&v);
	tmp = _mm_unpacklo_pd(tmp, tmp);
	return vec<double,2>(_mm_add_pd(xmmd[0], tmp)); 
}
template <> vec<double, 2> vec<double, 2>::operator- (const double &v) const	{ 
	__m128d tmp = _mm_load_sd(&v);
	tmp = _mm_unpacklo_pd(tmp, tmp);
	return vec<double,2>(_mm_sub_pd(xmmd[0], tmp)); 
}
template <> vec<double, 2> vec<double, 2>::operator/ (const double &v) const	{ 
	__m128d tmp = _mm_load_sd(&v);
	tmp = _mm_unpacklo_pd(tmp, tmp);
	return vec<double,2>(_mm_div_pd(xmmd[0], tmp)); 
}
template <> inline void vec<double, 2>::operator*= (const double &v) { 
	__m128d tmp = _mm_load_sd(&v);
	tmp = _mm_unpacklo_pd(tmp, tmp);
	xmmd[0] = _mm_mul_pd(xmmd[0], tmp); 
}
template <> inline void vec<double, 2>::operator+= (const double &v) { 
	__m128d tmp = _mm_load_sd(&v);
	tmp = _mm_unpacklo_pd(tmp, tmp);
	xmmd[0] = _mm_add_pd(xmmd[0], tmp); 
}
template <> inline void vec<double, 2>::operator-= (const double &v) { 
	__m128d tmp = _mm_load_sd(&v);
	tmp = _mm_unpacklo_pd(tmp, tmp);
	xmmd[0] = _mm_sub_pd(xmmd[0], tmp); 
}
template <> inline void vec<double, 2>::operator/= (const double &v) {
	__m128d tmp = _mm_load_sd(&v);
	tmp = _mm_unpacklo_pd(tmp, tmp);
	xmmd[0] = _mm_div_pd(xmmd[0], tmp); 
}

//  Operators, scalar versions (double, length 4):
template <> vec<double, 4> vec<double, 4>::operator* (const double &v) const	{ 
	__m128d tmp = _mm_load_sd(&v);
	tmp = _mm_unpacklo_pd(tmp, tmp);
	return vec<double,4>(_mm_mul_pd(xmmd[0], tmp), _mm_mul_pd(xmmd[1], tmp)); 
}
template <> vec<double, 4> vec<double, 4>::operator+ (const double &v) const	{ 
	__m128d tmp = _mm_load_sd(&v);
	tmp = _mm_unpacklo_pd(tmp, tmp);
	return vec<double,4>(_mm_add_pd(xmmd[0], tmp), _mm_add_pd(xmmd[1], tmp)); 
}
template <> vec<double, 4> vec<double, 4>::operator- (const double &v) const	{ 
	__m128d tmp = _mm_load_sd(&v);
	tmp = _mm_unpacklo_pd(tmp, tmp);
	return vec<double,4>(_mm_sub_pd(xmmd[0], tmp), _mm_sub_pd(xmmd[1], tmp)); 
}
template <> vec<double, 4> vec<double, 4>::operator/ (const double &v) const	{ 
	__m128d tmp = _mm_load_sd(&v);
	tmp = _mm_unpacklo_pd(tmp, tmp);
	return vec<double,4>(_mm_div_pd(xmmd[0], tmp), _mm_div_pd(xmmd[1], tmp)); 
}
template <> inline void vec<double, 4>::operator*= (const double &v) { 
	__m128d tmp = _mm_load_sd(&v);
	tmp = _mm_unpacklo_pd(tmp, tmp);
	xmmd[0] = _mm_mul_pd(xmmd[0], tmp); 
	xmmd[1] = _mm_mul_pd(xmmd[1], tmp); 
}
template <> inline void vec<double, 4>::operator+= (const double &v) { 
	__m128d tmp = _mm_load_sd(&v);
	tmp = _mm_unpacklo_pd(tmp, tmp);
	xmmd[0] = _mm_add_pd(xmmd[0], tmp); 
	xmmd[1] = _mm_add_pd(xmmd[1], tmp); 
}
template <> inline void vec<double, 4>::operator-= (const double &v) { 
	__m128d tmp = _mm_load_sd(&v);
	tmp = _mm_unpacklo_pd(tmp, tmp);
	xmmd[0] = _mm_sub_pd(xmmd[0], tmp); 
	xmmd[1] = _mm_sub_pd(xmmd[1], tmp); 
}
template <> inline void vec<double, 4>::operator/= (const double &v) {
	__m128d tmp = _mm_load_sd(&v);
	tmp = _mm_unpacklo_pd(tmp, tmp);
	xmmd[0] = _mm_div_pd(xmmd[0], tmp); 
	xmmd[1] = _mm_div_pd(xmmd[1], tmp); 
}

//  Operators, scalar versions of operators (float, length 8):
template <> vec<float, 8> vec<float, 8>::operator* (const float &v) const	{ 
	__m128 tmp = _mm_load_ss(&v);
	tmp = _mm_shuffle_ps(tmp,tmp,_MM_SHUFFLE(0, 0, 0, 0));
	return vec<float,8>(_mm_mul_ps(xmm[0], tmp), _mm_mul_ps(xmm[1], tmp)); 
}
template <> vec<float, 8> vec<float, 8>::operator+ (const float &v) const	{ 
	__m128 tmp = _mm_load_ss(&v);
	tmp = _mm_shuffle_ps(tmp,tmp,_MM_SHUFFLE(0, 0, 0, 0));
	return vec<float,8>(_mm_add_ps(xmm[0], tmp), _mm_add_ps(xmm[1], tmp)); 
}
template <> vec<float, 8> vec<float, 8>::operator- (const float &v) const	{ 
	__m128 tmp = _mm_load_ss(&v);
	tmp = _mm_shuffle_ps(tmp,tmp,_MM_SHUFFLE(0, 0, 0, 0));
	return vec<float,8>(_mm_sub_ps(xmm[0], tmp), _mm_sub_ps(xmm[1], tmp)); 
}
template <> vec<float, 8> vec<float, 8>::operator/ (const float &v) const	{ 
	__m128 tmp = _mm_load_ss(&v);
	tmp = _mm_shuffle_ps(tmp,tmp,_MM_SHUFFLE(0, 0, 0, 0));
	return vec<float,8>(_mm_div_ps(xmm[0], tmp), _mm_div_ps(xmm[1], tmp)); 
}
template <> inline void vec<float, 8>::operator*= (const float &v) { 
	__m128 tmp = _mm_load_ss(&v);
	tmp = _mm_shuffle_ps(tmp,tmp,_MM_SHUFFLE(0, 0, 0, 0));
	xmm[0] = _mm_mul_ps(xmm[0], tmp); 
	xmm[1] = _mm_mul_ps(xmm[1], tmp); 
}
template <> inline void vec<float, 8>::operator+= (const float &v) { 
	__m128 tmp = _mm_load_ss(&v);
	tmp = _mm_shuffle_ps(tmp,tmp,_MM_SHUFFLE(0, 0, 0, 0));
	xmm[0] = _mm_add_ps(xmm[0], tmp); 
	xmm[1] = _mm_add_ps(xmm[1], tmp); 
}
template <> inline void vec<float, 8>::operator-= (const float &v) { 
	__m128 tmp = _mm_load_ss(&v);
	tmp = _mm_shuffle_ps(tmp,tmp,_MM_SHUFFLE(0, 0, 0, 0));
	xmm[0] = _mm_sub_ps(xmm[0], tmp); 
	xmm[1] = _mm_sub_ps(xmm[1], tmp); 
}
template <> inline void vec<float, 8>::operator/= (const float &v) {
	__m128 tmp = _mm_load_ss(&v);
	tmp = _mm_shuffle_ps(tmp,tmp,_MM_SHUFFLE(0, 0, 0, 0));
	xmm[0] = _mm_div_ps(xmm[0], tmp); 
	xmm[1] = _mm_div_ps(xmm[1], tmp); 
}

// Operators,  scalar versions (float, length 4):
template <> vec<float, 4> vec<float, 4>::operator* (const float &v) const	{ 
	__m128 tmp = _mm_load_ss(&v);
	tmp = _mm_shuffle_ps(tmp,tmp,_MM_SHUFFLE(0, 0, 0, 0));
	return vec<float,4>(_mm_mul_ps(xmm[0], tmp)); 
}
template <> vec<float, 4> vec<float, 4>::operator+ (const float &v) const	{ 
	__m128 tmp = _mm_load_ss(&v);
	tmp = _mm_shuffle_ps(tmp,tmp,_MM_SHUFFLE(0, 0, 0, 0));
	return vec<float,4>(_mm_add_ps(xmm[0], tmp)); 
}
template <> vec<float, 4> vec<float, 4>::operator- (const float &v) const	{ 
	__m128 tmp = _mm_load_ss(&v);
	tmp = _mm_shuffle_ps(tmp,tmp,_MM_SHUFFLE(0, 0, 0, 0));
	return vec<float,4>(_mm_sub_ps(xmm[0], tmp)); 
}
template <> vec<float, 4> vec<float, 4>::operator/ (const float &v) const	{ 
	__m128 tmp = _mm_load_ss(&v);
	tmp = _mm_shuffle_ps(tmp,tmp,_MM_SHUFFLE(0, 0, 0, 0));
	return vec<float,4>(_mm_div_ps(xmm[0], tmp)); 
}
template <> inline void vec<float, 4>::operator*= (const float &v) { 
	__m128 tmp = _mm_load_ss(&v);
	tmp = _mm_shuffle_ps(tmp,tmp,_MM_SHUFFLE(0, 0, 0, 0));
	xmm[0] = _mm_mul_ps(xmm[0], tmp); 
}
template <> inline void vec<float, 4>::operator+= (const float &v) { 
	__m128 tmp = _mm_load_ss(&v);
	tmp = _mm_shuffle_ps(tmp,tmp,_MM_SHUFFLE(0, 0, 0, 0));
	xmm[0] = _mm_add_ps(xmm[0], tmp); 
}
template <> inline void vec<float, 4>::operator-= (const float &v) { 
	__m128 tmp = _mm_load_ss(&v);
	tmp = _mm_shuffle_ps(tmp,tmp,_MM_SHUFFLE(0, 0, 0, 0));
	xmm[0] = _mm_sub_ps(xmm[0], tmp); 
}
template <> inline void vec<float, 4>::operator/= (const float &v) {
	__m128 tmp = _mm_load_ss(&v);
	tmp = _mm_shuffle_ps(tmp,tmp,_MM_SHUFFLE(0, 0, 0, 0));
	xmm[0] = _mm_div_ps(xmm[0], tmp); 
}

// Operators, Specialized versions (int64_t precision, length 2):
template <> vec<int64_t, 2> vec<int64_t, 2>::operator+ (const vec<int64_t, 2> &v) const	{ 
	return vec<int64_t, 2>(  _mm_add_epi64(xmmi[0], v.xmmi[0]) ); 
}
template <> vec<int64_t, 2> vec<int64_t, 2>::operator- (const vec<int64_t, 2> &v) const	{ 
	return vec<int64_t, 2>(  _mm_sub_epi64(xmmi[0], v.xmmi[0]) ); 
}
int64_t ALIGN16 negation64[2] = {0x8000000000000000,0x8000000000000000};
template <> vec<int64_t, 2> vec<int64_t, 2>::operator>> (const int shiftcount) const	{ 
	// shift right not supported for signed int64 values. (ARGH!!, why???)
	__m128i tmp = _mm_load_si128( (__m128i *) &negation64[0] );
	//tmp.m128i_u64[0] = 0x8000000000000000;
	//tmp.m128i_u64[1] = 0x8000000000000000;
	return vec<int64_t, 2>( _mm_sub_epi64( _mm_srli_epi64( _mm_add_epi64( xmmi[0], tmp), shiftcount), _mm_srli_epi64( tmp, shiftcount) ) ); 
}
template <> vec<int64_t, 2> vec<int64_t, 2>::operator<< (const int shiftcount) const	{ 
	return vec<int64_t, 2>( _mm_slli_epi64(xmmi[0], shiftcount) ); 
}
template <> vec<int64_t, 2> vec<int64_t, 2>::operator| (const vec<int64_t, 2> &v) const	{ 
	return vec<int64_t, 2>( _mm_or_si128(xmmi[0], v.xmmi[0]) ); 
}
template <> vec<int64_t, 2> vec<int64_t, 2>::operator& (const vec<int64_t, 2> &v) const	{ 
	return vec<int64_t, 2>( _mm_and_si128(xmmi[0], v.xmmi[0]) ); 
}
#ifdef INCLUDE_SSE4
template <> vec<int64_t, 2> vec<int64_t, 2>::operator>(const vec<int64_t, 2> &v) const	{ 
	return vec<int64_t, 2>( _mm_cmpgt_epi64(xmmi[0], v.xmmi[0]) ); 
}
/* Somehow lt comparison of 64 bit integers is not supported (which idiot decided that?)
// obviously !((a>b) | (a==b)) is the same, but much slower.
template <> vec<int64_t, 2> vec<int64_t, 2>::operator<(const vec<int64_t, 2> &v) const	{ 
	return vec<int64_t, 2>( _mm_cmplt_epi64(xmmi[0], v.xmmi[0]) ); 
}*/
template <> vec<int64_t, 2> vec<int64_t, 2>::operator==(const vec<int64_t, 2> &v) const	{ 
	return vec<int64_t, 2>( _mm_cmpeq_epi64(xmmi[0], v.xmmi[0]) ); 
}
#endif

// repeat element
template <> vec<double, 4> vec<double, 4>::rep (int idx) const	{ 
	__m128d tmp = ((idx & 2) ? xmmd[1] : xmmd[0]);
	tmp = ( (idx &1 ) ? _mm_unpackhi_pd(tmp, tmp) : _mm_unpacklo_pd(tmp, tmp));
	return vec<double,4>(tmp, tmp); 
}
template <> vec<float, 4> vec<float, 4>::rep (int idx) const	{ 
	__m128 tmp = ((idx & 2) ? _mm_unpackhi_ps(xmm[0], xmm[0]) : _mm_unpacklo_ps(xmm[0], xmm[0]));
	tmp = ( (idx &1 ) ? _mm_unpackhi_ps(tmp, tmp) : _mm_unpacklo_ps(tmp, tmp));
	return vec<float,4>(tmp); 
}


// other functions (min/max ..) (double, length 2)
inline vec<double, 2> max(const vec<double, 2> &a, const vec<double, 2> &b){
	return vec<double,2>(_mm_max_pd(a.xmmd[0], b.xmmd[0])); 
}
inline vec<double, 2> min(const vec<double, 2> &a, const vec<double, 2> &b){
	return vec<double,2>(_mm_min_pd(a.xmmd[0], b.xmmd[0])); 
}
#ifdef INCLUDE_SSE3
inline double accum(const vec<double, 2> &v){
	__m128d tmp = _mm_hadd_pd(v.xmmd[0],v.xmmd[0]);
	double tmpd;
	_mm_store_sd(&tmpd,tmp);
	return tmpd; 
}
#endif
#ifdef INCLUDE_SSE4
inline vec<double, 2> round(const vec<double, 2> &v){
	return vec<double,2>(_mm_round_pd(v.xmmd[0], _MM_FROUND_TO_NEAREST_INT )); 
}
#endif
// other functions (min/max ..) (double, length 4)
inline vec<double, 4> max(const vec<double, 4> &a, const vec<double, 4> &b){
	return vec<double,4>(_mm_max_pd(a.xmmd[0], b.xmmd[0]), _mm_max_pd(a.xmmd[1], b.xmmd[1])); 
}
inline vec<double, 4> min(const vec<double, 4> &a, const vec<double, 4> &b){
	return vec<double,4>(_mm_min_pd(a.xmmd[0], b.xmmd[0]), _mm_min_pd(a.xmmd[1], b.xmmd[1])); 
}
#ifdef INCLUDE_SSE3
inline double accum(const vec<double, 4> &v){
	__m128d tmp = _mm_add_pd(v.xmmd[0], v.xmmd[1]);
	tmp = _mm_hadd_pd(tmp,tmp);
	double tmpd;
	_mm_store_sd(&tmpd,tmp);
	return tmpd; 
}
#endif
#ifdef INCLUDE_SSE4
inline vec<double, 4> round(const vec<double, 4> &v){
	return vec<double,4>(_mm_round_pd(v.xmmd[0], _MM_FROUND_TO_NEAREST_INT ), _mm_round_pd(v.xmmd[0], _MM_FROUND_TO_NEAREST_INT )); 
}
#endif
// other functions (min/max ..) (float, length 8)
inline vec<float, 8> max(const vec<float, 8> &a, const vec<float, 8> &b){
	return vec<float,8>(_mm_max_ps(a.xmm[0], b.xmm[0]), _mm_max_ps(a.xmm[1], b.xmm[1])); 
}
inline vec<float, 8> min(const vec<float, 8> &a, const vec<float, 8> &b){
	return vec<float,8>(_mm_min_ps(a.xmm[0], b.xmm[0]), _mm_min_ps(a.xmm[1], b.xmm[1])); 
}
#ifdef INCLUDE_SSE3
inline float accum(const vec<float, 8> &v){
	__m128 tmp = _mm_add_ps(v.xmm[0], v.xmm[1]);
	tmp = _mm_hadd_ps(tmp,tmp);
	tmp = _mm_hadd_ps(tmp,tmp);
	float tmpd;
	_mm_store_ss(&tmpd,tmp);
	return tmpd; 
}
#endif
// other functions (min/max ..) (float, length 4)
inline vec<float, 4> max(const vec<float, 4> &a, const vec<float, 4> &b){
	return vec<float,4>(_mm_max_ps(a.xmm[0], b.xmm[0])); 
}
inline vec<float, 4> min(const vec<float, 4> &a, const vec<float, 4> &b){
	return vec<float,4>(_mm_min_ps(a.xmm[0], b.xmm[0])); 
}
#ifdef INCLUDE_SSE3
inline float accum(const vec<float, 4> &v){
	__m128 tmp = _mm_hadd_ps(v.xmm[0], v.xmm[0]);
	tmp = _mm_hadd_ps(tmp,tmp);
	float tmpd;
	_mm_store_ss(&tmpd,tmp);
	return tmpd; 
}
#endif
inline vec<float, 4> abs(const vec<float, 4> &a){
    __m128i tmp = _mm_set1_epi32( (0x7FFFFFFF) );
	return vec<float,4>(_mm_and_ps(a.xmm[0], _mm_castsi128_ps(tmp))); 
}
#ifdef __x86_64__ 
inline vec<double, 2> abs(const vec<double, 2> &a){
#ifdef INCLUDE_SSE4
	__m128i tmp;
	tmp =  _mm_srl_epi64( _mm_cmpeq_epi64(tmp,tmp) , 1);
#else
    __m128i tmp = _mm_set1_epi64x( (0x7FFFFFFFFFFFFFFF) );
#endif
	return vec<double,2>(_mm_and_pd(a.xmmd[0], _mm_castsi128_pd(tmp))); 
}
#endif
inline vec<double, 2> operator^(const vec<double, 2> &a, const vec<double, 2> &b){
	return vec<double,2>(_mm_xor_pd(a.xmmd[0], b.xmmd[0])); 
}
template <> vec<double, 2> vec<double, 2>::operator& (const vec<double, 2> &v) const	{ 
	return vec<double,2>(_mm_and_pd(xmmd[0], v.xmmd[0])); 
}
template <> vec<double, 2> vec<double, 2>::operator| (const vec<double, 2> &v) const	{ 
	return vec<double,2>(_mm_or_pd(xmmd[0], v.xmmd[0])); 
}
inline vec<double, 2> andnot(const vec<double, 2> &a, const vec<double, 2> &b){
	return vec<double,2>(_mm_andnot_pd(a.xmmd[0], b.xmmd[0])); 
}

template <typename T, int vlen> vec<T, vlen> vec<T, vlen>::insert( const vec<T, vlen> &b, const vec<T, vlen> &mask ) { 
	return andnot(mask, *this ) | (b & mask);
}


template <typename T>  
inline void conditional_swap(T &a, T &b, const T condition){
/* for all bits in condition that are 1 (i.e. when condition is true)
   the corresponding bits of a and b are swapped. 
   Example:
   	  conditional_swap(a,b, a>b) 
   performs for all elements in a and b: if a[i]>b[i] then swap(a[i],b[i]) 
   so after evaluation of this function a[i]<=b[i]
   */
     
    T tmp = ((a^b)& condition);
    // if true (i.e. all bits set)  : tmp = a xor b
    //   			           else : tmp = 0
    a = (a^tmp);
    b = (b^tmp);
}

template <typename T>  inline T sqr(T a) {
/* Sqr: compute square. For vector arguments, the element wise square is computed.
    Created by Dirk Poot, Erasmus MC
   */
   return a*a;  
}

#ifdef INCLUDE_SSE4
inline vec<float, 4> round(const vec<float, 4> &v){
	return vec<float,4>(_mm_round_ps(v.xmm[0], _MM_FROUND_TO_NEAREST_INT )); 
}
// other functions (min/max ..) (int64, length 2)
inline vec<int64_t, 2> max_bad(const vec<int64_t, 2> &a, const vec<int64_t, 2> &b){
    // _bad : using 32 bit maximum function
	return vec<int64_t,2>(_mm_max_epi32(a.xmmi[0], b.xmmi[0])); 
}
inline vec<int64_t, 2> min_bad(const vec<int64_t, 2> &a, const vec<int64_t, 2> &b){
    // _bad : using 32 bit maximum function
	return vec<int64_t,2>(_mm_min_epi32(a.xmmi[0], b.xmmi[0])); 
}

// other functions (min/max ..) (int32, length 4)
inline vec<int32_t, 4> max(const vec<int32_t, 4> &a, const vec<int32_t, 4> &b){
	return vec<int32_t,4>(_mm_max_epi32(a.xmmi[0], b.xmmi[0])); 
}
inline vec<int32_t, 4> min(const vec<int32_t, 4> &a, const vec<int32_t, 4> &b){
	return vec<int32_t,4>(_mm_min_epi32(a.xmmi[0], b.xmmi[0])); 
}
#endif


#undef BAD

/* ORIGINAL:

struct vec4
{
  __m128 xmm;

  vec4 (__m128 v) : xmm (v) {}

  vec4 (float v) { xmm = _mm_set1_ps(v); }

  vec4 (float x, float y, float z, float w)
  { xmm = _mm_set_ps(w,z,y,x); }

  vec4 (const float *v) { xmm = _mm_load_ps(v); }

  vec4 operator* (const vec4 &v) const
  { return vec4(_mm_mul_ps(xmm, v.xmm)); }

  vec4 operator+ (const vec4 &v) const
  { return vec4(_mm_add_ps(xmm, v.xmm)); }

  vec4 operator- (const vec4 &v) const
  { return vec4(_mm_sub_ps(xmm, v.xmm)); }

  vec4 operator/ (const vec4 &v) const
  { return vec4(_mm_div_ps(xmm, v.xmm)); }

  void operator*= (const vec4 &v)
  { xmm = _mm_mul_ps(xmm, v.xmm); }

  void operator+= (const vec4 &v)
  { xmm = _mm_add_ps(xmm, v.xmm); }

  void operator-= (const vec4 &v)
  { xmm = _mm_sub_ps(xmm, v.xmm); }

  void operator/= (const vec4 &v)
  { xmm = _mm_div_ps(xmm, v.xmm); }

  void operator>> (float *v)
  { _mm_store_ps(v, xmm); }

};
*/
#endif