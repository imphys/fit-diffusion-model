#ifndef LINEREADWRITERS
#define LINEREADWRITERS
/* lineReadWriters.cpp
 * 
 * Helper classes to read elements from 'lines'.
 * Each class is templated over image type and indexing type.
 * Three methods:
 *  one creator function, takes arguments applicable for each type.
 *  - value = read( i )     : reads 'element' i into 'value'
 *  - write( i , value)     : writes 'value' to 'element' i 
 *  - accumulate(i , value) : adds value to element i
 * 
 * Naming: 
 * lineType[2][c]
 * [2] : indicates 2D 'lines'. 'all' values in second dim are placed in a vector, 
 *	      creator function takes 'step_sdim' as argument and vlen as extra template argument
 *       valuetype = vec<imgtype, vlen>
 * [c] : indicates complex values. Creator function takes imR & imI instead of im.
 *		 type returned by read and required by write = complex<valuetype>
 *
 * An extra lineTypeAccum type is provided, which accumulates values instead of just writing.
 * This class is templated over the other lineTypes.
 *
 * You might want to use step_iterators.cpp instead
 *       
 */

#include <iterator>
#include <complex>
using std::complex;
/*
template<typename base> class lineTypeAccum : public base {
public:
	inline void write( typename base::int_t i, typename base::valuetype newval) {
		base::write( i, read( i ) + newval );
	}
};
*/


template < typename base > class assignmenthelper {
public:
	typedef typename base::int_T int_t ;
	typedef typename base::value_type value_type;
private:
	base & obj;
	int_t idx;

public:
	assignmenthelper( base & obj_, int_t idx_) : obj(obj_), idx(idx_) {};
	inline void operator=(value_type newval) {
		obj.write(idx,newval);
	}
	inline void operator+=(value_type newval) {
		obj.accumulate(idx,newval);
	}
	/*inline value_type operator*() {
		return obj.read(idx);
	}*/
	inline operator value_type() { 
		return obj.read(idx); 
	};
};

template < typename base > class assignmenthelper0 {
public:
	typedef typename base::value_type value_type;
private:
	base obj;
public:
	assignmenthelper0( base * obj_) : obj(*obj_){};
	assignmenthelper0( base & obj_) : obj(obj_){};
	inline void operator=(value_type newval) {
		obj.write(newval);
	}
	inline void operator+=(value_type newval) {
		obj.accumulate(newval);
	}
	inline void operator-=(value_type newval) {
		obj.accumulate(-newval);
	}
	/*inline value_type operator*() {
		return obj.read();
	}*/
	inline operator value_type() { 
		return obj.read(); 
	};
};
// base class that implements the common stuff (but doesnt store anything).
// derived needs to:
//  -  store the pointers (& other required info)
//  -  implement the following methods:
//    inline void step( int_t i );			// inline to optimize for i==1 or i==-1
//    value_type read();					// for value = *it, value = it[k], 
//    void write( const value_type value)		// for *it = newval, it[k]=newval,
//    bool operator==( const self testvalue)			// for == and !=
//   optionally:  void accumulate( const value_type )   // for *it += newval, it[k] += newval, *it-=newval, it[k]-=newval.
//			    (a default implementation using read and write is provided)
//  - be copy-able (i.e. fast to copy)
template < typename derived , typename valueT = typename derived::value_type , typename intT = typename derived::int_t > class line_type_base : public std::iterator< std::random_access_iterator_tag , valueT , intT > { 
public:
	typedef derived self;
	typedef valueT value_type;
	//typedef typename std::iterator_traits< typename derived::imgtype >::value_type imgValType;
    typedef self pointer;
    typedef assignmenthelper0<self> reference;
	friend class assignmenthelper0<self>;

	typedef intT int_t ;
public:
	line_type_base() {};
/*	inline bool operator==(const self & b) const {
		return (curindex == b.curindex);
	};*/
	inline bool operator!=(const self & b) const {
		return !(static_cast<derived*>(this) == b);
	};
	inline const value_type operator*() const {
		return static_cast<derived*>(this)->read();
	};
	inline reference operator*() {
		return reference(static_cast<derived*>(this));
	};
	inline derived operator+(const int_t b) const {
		derived temp( *(static_cast< const derived * >(this)) );
		temp.step(b);
		return temp;
	}
	inline derived operator - (const  int_t b) const {
		derived temp( static_cast<const derived * const>(this) );
		temp.step(-b);
		return temp;
	}
	inline derived * operator++() {
		(static_cast<derived*>(this))->step( 1 );
		return static_cast<derived*>(this);
	};
	inline derived operator++(int) {
		derived tmp = static_cast<derived*>(this);
		static_cast<derived*>(this)->step( 1 ) ;
		return tmp;
	};
	inline derived & operator--() {
		return (static_cast<derived*>(this))->step(-1);
	};
	inline derived operator--(int) {
		derived tmp = static_cast<derived*>(this);
		static_cast<derived*>(this)->step( -1 ) ;
		return tmp;
	};
	inline derived * operator+=(int_t cnt) {
		return (static_cast<derived*>(this))->step( cnt ) ;
	};
	inline derived * operator-=(int_t cnt) {
		return static_cast<derived*>(this)->step( -cnt ) ;
	};
	inline value_type operator[]( int_t i ) const {
		derived tmp( static_cast<derived*>(this) );
		tmp.step(i);
		return tmp.read();
	};
	inline reference operator[](int_t i) {
		derived tmp( *(static_cast<derived*>(this)) );
		tmp.step(i);
		return reference(tmp);
	}
private:
	void accumulate( const value_type value) { // default accumulate implementation.
		static_cast<derived*>(this)->write( static_cast<derived*>(this)->read() + value );
	}
};

template<typename imgtype, typename int_t> class lineType0  { 
private:
	imgtype* im;                                          
	int_t   length_;
public:
	typedef imgtype value_type;
	typedef value_type * iterator;
	typedef int_t int_T;

	lineType0( imgtype * im_ , int_t len) {
		im = im_;
		length_ = len;
	}
	inline imgtype read(int_t i) {
		return im[i];
	};
	inline void write(  int_t i, const imgtype newval) {
		im[ i ] = newval;
	};
	inline void accumulate(  int_t i, const imgtype addval) {
		im[ i ] += addval;
	};

	inline imgtype * operator[](int_t i) {
		return &im[i];
	};

	// value_type handles all assignment of values

	imgtype * begin() { // iterator begin.
		return im;
	};
	imgtype * end() { // iterator end.
		return im+length_;
	};
};


template<typename imgtype, typename valueT= typename std::iterator_traits< imgtype >::value_type , typename intT = ptrdiff_t> class lineType : public line_type_base< lineType<imgtype, valueT, intT> , valueT, intT>  { 
public:
	typedef line_type_base< lineType<imgtype, valueT, intT> , valueT, intT> parent;
	typedef typename parent::int_t int_t;
private:
	imgtype im;                                                      
	int_t step_Tdim;
public:
	typedef valueT value_type;
	lineType() {
		im = NULL;
		step_Tdim = 0;
	};
	lineType( const lineType * init ): im(init->im), step_Tdim(init->step_Tdim) {;};
	lineType( imgtype im_, int_t step_Tdim_ ) : im(im_), step_Tdim(step_Tdim_) {;};
	inline void step( int_t i) {
		im += i*step_Tdim;
	};
	inline value_type read() {
		return *im;
	};
	inline void write(value_type newval) {
		*im =newval;
	};
	inline void accumulate(value_type newval) {
		*im +=newval;
	};
/*	inline value_type read(int_t i) {
		return (valuetype) im[i*step_Tdim];
	};
	inline void write(  int_t i, const  value_type  newval) {
		im[ i * step_Tdim ] = (imgValType) newval;
	};
	inline void accumulate(  int_t i, const value_type  addval) {
		im[ i * step_Tdim ] += (imgValType) addval;
	};*/
};

template<typename imgtype, typename valuetype, typename int_t> class lineTypeL : lineType<imgtype, valuetype, int_t> { 
public:
	typedef lineType<imgtype, valuetype, int_t> parent;
private:
	int_t   length_;
public:
	lineTypeL( imgtype* im_, int_t step_Tdim_ , int_t length) : parent(im_, step_Tdim_) {
		length_ = length;
	}
};


template<typename imgtype, typename valueT, typename intT> class lineTypec : public line_type_base< lineTypec<imgtype, valueT, intT> , complex<valueT>, intT>   { 
private:
	imgtype imR;
	imgtype imI;
	intT step_Tdim;
public:
	typedef lineTypec<imgtype, valueT, intT> self;
	typedef typename std::iterator_traits< imgtype >::value_type imgValType;
/*	typedef typename std::random_access_iterator_tag iterator_category;*/
	typedef complex<valueT> value_type;
/*	typedef ptrdiff_t difference_type;
    typedef difference_type distance_type;
    typedef int * pointer;
    typedef assignmenthelper<self> reference;
	typedef intT int_T;*/
    lineTypec(){};// imR( NULL ), imI( NULL), step_Tdim(0) {;}; // default constructor, initialize to null.
	lineTypec( const lineTypec & init): imR(init.imR), imI(init.imI), step_Tdim(init.step_Tdim) {;};
	lineTypec( imgtype imR_, imgtype imI_, intT step_Tdim_ ) : imR(imR_), imI(imI_), step_Tdim(step_Tdim_) {};

	inline self * step( intT i ) {
		imR += i*step_Tdim;
		imI += i*step_Tdim;
		return this;
	}
	inline value_type read() {
		return complex<valueT>( (valueT) *imR, (valueT) *imI);
	};
	inline void write( const value_type newval) {
		*imR =  (imgValType) newval.real();
		*imI =  (imgValType) newval.imag();
	};
	/*bool operator==( const self testval) {
		return imR==testval.imR; // assume if Real part is equal, all pointers are equal.
	}*/

/*	inline complex<valueT> read(intT i) {
		return complex<valueT>( (valueT) imR[i*step_Tdim], (valueT) imI[i*step_Tdim]);
	};
	inline void write(  intT i, const complex<valueT> newval) {
		imR[i*step_Tdim] =  (imgValType) newval.real();
		imI[i*step_Tdim] =  (imgValType) newval.imag();
	};
	inline void accumulate(  intT i, const complex<valueT> newval) {
		imR[i*step_Tdim] +=  (imgValType) newval.real();
		imI[i*step_Tdim] +=  (imgValType) newval.imag();
	};
	inline value_type operator[](intT i) const {
		return read(i);
	};
	inline reference operator[](intT i) {
		return reference(*this, i);
	};*/
};

template<typename imgtype, typename valuetype, typename int_t> class lineType2  { 
private:
	imgtype im;                                                      
	int_t step_Tdim;
	int_t step_sdim;
public:
	typedef valuetype value_type;

	lineType2( imgtype im_, int_t step_Tdim_, int_t step_sdim_) {
		im = im_; step_Tdim = step_Tdim_; step_sdim =step_sdim_;
	}
	inline valuetype read(int_t i) {
		return valuetype( im + i*step_Tdim , step_sdim);
	};
	inline void write(  int_t i, valuetype  newval) {
		newval.store( im  + i*step_Tdim , step_sdim);
	};
	inline void accumulate(  int_t i, valuetype  newval) {
		newval += read( i );
		newval.store( im  + i*step_Tdim , step_sdim);
	};
};

template<typename imgtype, typename valuetype, typename int_t> class lineType2c  { 
private:
	imgtype imR;
	imgtype imI;
	int_t step_Tdim;
	int_t step_sdim;
public:
	typedef complex<valuetype> value_type;

	lineType2c( imgtype imR_, imgtype* imI_, int_t step_Tdim_ , int_t step_sdim_) {
		imR = imR_; imI = imI_; step_Tdim = step_Tdim_; step_sdim =step_sdim_;
	}
	inline complex<valuetype> read(int_t i) {
		return complex<valuetype>( valuetype( imR + i*step_Tdim , step_sdim), valuetype( imI + i*step_Tdim , step_sdim));
	};
	inline void write(  int_t i, const complex<valuetype> newval) {
		newval.real().store( imR  + i*step_Tdim , step_sdim);
		newval.imag().store( imI  + i*step_Tdim , step_sdim);
	};
	inline void accumulate(  int_t i, const complex<valuetype> newval) {
		newval += read( i );
		newval.real().store( imR  + i*step_Tdim , step_sdim);
		newval.imag().store( imI  + i*step_Tdim , step_sdim);
	};
};

#endif