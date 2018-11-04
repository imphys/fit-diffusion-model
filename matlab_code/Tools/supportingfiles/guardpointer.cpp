#ifndef GUARDEDPOINTER_CPP
#define GUARDEDPOINTER_CPP
/*  guardpointer.cpp
  Code that allows range checking on pointers, to avoid out of bounds accesses.
  When enabled (compile time), an error is thrown when out of bound elements are accessed.

  Enable bound checking by defining 'CHECKPOINTERBOUNDS' (#define CHECKPOINTERBOUNDS)

  Has 2 modi, - checks pointers you already have  (single function call upon creation when CHECKPOINTERBOUNDS not defined)
		      - check mxArrays  (convenience method that avoids even the single function call)
  Useage:
  - for pointers:
      typedef typename guard_pointer_type< data_pointer_type >::type  guarded_data_pointer_type;
	  guarded_data_pointer_type  guarded_pointer( first, end, currentpos);

	  currentpos can be ommitted, defaults to first;

  - for mxArrays call the macro:
	  mxGetPtr( data_pointer_type, varname , sourcename);
	  or : mxGetComplexPtr( data_pointer_type, varname, sourename)
 
    This does define a pointer type:
      varname_pointer_type;  (if CHECKPOINTERBOUNDS is not defined: data_pointer_type == varname_type )
    and (effectively) performs
	  varname_type varname = mxGetData( sourcename );
      or: varname_type varname_R = mxGetData( sourcename );
	      varname_type varname_I = mxGetImagData( sourcename );

	Example: 
      mxGetPtr( const double * , test_ptr , prhs[0] );

	  which defines a variable 'test_ptr' with a type called 'test_ptr_type' which is initialized 
	  to point to the first element of prhs[0], which is assumed to be an array of (non complex) double precision values.


  Use iterator_traits and templates to make your programming independed of the actual type of varname_pointer_type;

  When CHECKPOINTERBOUNDS is defined the pointer type performs the range checking (costs a bit of computation time)
  otherwise a standard pointer is returned (so absolutely no computation cost is added)

  This routine only provides the range checking; it does not check whether the mx arrays are of the specified type.

  Created by Dirk Poot, Erasmus MC, 4-3-2013
*/

#ifdef EMM_VEC_CPP 
#error Import guardpointer.cpp prior to emm_vec.cpp
#endif

#ifdef CHECKPOINTERBOUNDS
#pragma message("Guarding pointers that are created with guard_pointer or mxGet(Complex)Ptr .")
#define GUARDPOINTER_DEFINETYPE( DATATYPE, VARNAME ) typedef guard_pointer< DATATYPE > VARNAME ## _type
#define GUARDPOINTER_GETPOINTER( DATATYPE, SOURCENAME, MXFUN ) ( DATATYPE MXFUN( SOURCENAME ), DATATYPE ( (char *) MXFUN( SOURCENAME ) + mxGetNumberOfElements( SOURCENAME ) * mxGetElementSize( SOURCENAME ) ) )
// check guard class
#include <iterator>

	using std::iterator_traits;
template <typename T> class guard_pointer  {
public:
	typedef guard_pointer<T> self;
	typedef self type;
	typedef typename iterator_traits< T >::value_type        value_type;
	typedef typename iterator_traits< T >::iterator_category iterator_category;
	typedef typename iterator_traits< T >::difference_type   difference_type ;
	typedef typename iterator_traits< T >::pointer           pointer;
	typedef typename iterator_traits< T >::reference         reference;
private:
	T curpos;
	T first;
	T end;
	void check( T value ) {
		if (!((value>=first) && (value<end))) {
			mexErrMsgTxt("Out of range access detected in guard_pointer");
		}
	}
	friend class guard_pointer;
public:
	guard_pointer(T first_, T end_ ): curpos(first_), first(first_), end(end_) {;};
	//guard_pointer(T & first_, T & end_ , T & curpos_): curpos(curpos_),first(first_),end(end_) {;};
	guard_pointer(T  first_, T  end_ , T curpos_): curpos(curpos_),first(first_),end(end_) {;};
	// Typecasting pointers:
	template < typename otherT > guard_pointer( guard_pointer<otherT> const & other ): curpos(other.curpos),first(other.first),end(other.end) {;};
	inline bool operator!=(const self & b) const {
		return (curpos != b.curpos);
	};
	inline value_type operator*() const {
		check(curpos);
		return *curpos;
	};
	inline reference operator*() {
		check(curpos);
		return *curpos;
	};
	inline self operator+(const difference_type b) const {
		return self(first, end, curpos+b) ;
	}
	inline self operator-(const difference_type b) const {
		return self(first, end, curpos-b) ;
	}
	inline self* operator++() {
		++curpos;
		return this;
	};
	inline self operator++(int) {
		self tmp(first, end, curpos);
		curpos++;
		return tmp;
	};
	inline self* operator--() {
		--curpos;
		return this;
	};
	inline self operator--(int) {
		self tmp(first, end,curpos);
		curpos--;
		return tmp;
	};
	inline self* operator+=( difference_type cnt) {
		curpos+=cnt;
		return this;
	};
	inline self* operator-=( difference_type cnt) {
		curpos-=cnt;
		return this;
	};
	inline value_type operator[]( difference_type i ) const {
		check(curpos+i);
		return curpos[i];
	};
	inline reference operator[]( difference_type i) {
		check(curpos+i);
		return curpos[i];
	}
};

template <typename T> class guard_pointer_type  {
public:
	typedef guard_pointer<T> type;
};


#else
template <typename T> class guard_pointer_type  {
public:
	typedef T type;
};
template<typename T> inline T guard_pointer( T first, T end ) {
	return first;
};
template<typename T> inline T guard_pointer( T first, T end , T curpos) {
	return curpos;
};

#define GUARDPOINTER_DEFINETYPE( DATATYPE, VARNAME ) typedef DATATYPE  VARNAME ## _type
#define GUARDPOINTER_GETPOINTER( DATATYPE, SOURCENAME, MXFUN ) ( DATATYPE  MXFUN( SOURCENAME ) )

#endif

#define mxGetPtr( DATATYPE, VARNAME, SOURCENAME)        GUARDPOINTER_DEFINETYPE( DATATYPE, VARNAME );\
												        VARNAME ## _type VARNAME GUARDPOINTER_GETPOINTER( (DATATYPE), SOURCENAME, mxGetData)
#define mxGetPtri( DATATYPE, VARNAME, SOURCENAME, INTERMEDTYPE) GUARDPOINTER_DEFINETYPE( DATATYPE, VARNAME );\
												        VARNAME ## _type VARNAME GUARDPOINTER_GETPOINTER( (DATATYPE) (INTERMEDTYPE), SOURCENAME,  mxGetData)
#define mxGetComplexPtr( DATATYPE, VARNAME, SOURCENAME) GUARDPOINTER_DEFINETYPE( DATATYPE, VARNAME );\
	                                                    VARNAME ## _type VARNAME ## _R GUARDPOINTER_GETPOINTER( (DATATYPE), SOURCENAME, mxGetData);\
	                                                    VARNAME ## _type VARNAME ## _I GUARDPOINTER_GETPOINTER( (DATATYPE), SOURCENAME, mxGetImagData)
#define mxGetComplexPtri( DATATYPE, VARNAME, SOURCENAME, INTERMEDTYPE) GUARDPOINTER_DEFINETYPE( DATATYPE, VARNAME );\
	                                                    VARNAME ## _type VARNAME ## _R GUARDPOINTER_GETPOINTER( (DATATYPE) (INTERMEDTYPE), SOURCENAME,  mxGetData);\
	                                                    VARNAME ## _type VARNAME ## _I GUARDPOINTER_GETPOINTER( (DATATYPE) (INTERMEDTYPE), SOURCENAME, mxGetImagData)

#endif