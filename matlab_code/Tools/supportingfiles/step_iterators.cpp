#ifndef STEP_ITERATORS
#define STEP_ITERATORS
/* step_iterators.cpp
 * 
 * Helper classes to read elements from 'lines'.
 * Each class is templated over image type and indexing type.
 * Three methods:
 *  one creator function, takes arguments applicable for each type.
 *  
 * 
 * Naming: 
 * step_iterator[c]< pointertype, itegertype >
 * [c] : indicates complex values. Creator function takes imR & imI instead of im.
 *		 type returned by read and required by write = complex< valuetype >
 *
 * Created by Dirk Poot, Erasmus MC, 23-9-2013
 */

#include <iterator>
#include <complex>
using std::complex;


template < typename pointerType > class step_iterator_assignmenthelper {
public:
	typedef typename std::iterator_traits< pointerType >::value_type value_type;
private:
	pointerType obj;
public:
	step_iterator_assignmenthelper( pointerType obj_ ) : obj(obj_){};
	inline void operator=(value_type newval) {
		*obj = newval;
	}
	inline void operator+=(value_type newval) {
		*obj += newval;
	}
	inline void operator-=(value_type newval) {
		*obj -= newval;
	}
	inline value_type operator+( value_type val ) {
		return *obj + val;
	}
	inline value_type operator-( value_type val ) {
		return *obj - val;
	}
	inline value_type operator*( value_type val ) {
		return *obj * val;
	}
	inline value_type operator/( value_type val ) {
		return *obj / val;
	}
	template < typename otherT> inline void operator*=( otherT val ) {
		*obj *= val;
	}
	template < typename otherT > value_type operator*( otherT val ) const {
		return (*obj) * val;
	}
	inline operator value_type() { 
		return *obj; 
	};
};


template < typename pointerType, typename intT = int> class step_iterator : public std::iterator< std::random_access_iterator_tag , typename std::iterator_traits< pointerType >::value_type , intT > { 
public:
	typedef step_iterator<pointerType, intT> self;
	typedef typename std::iterator_traits< pointerType >::value_type value_type;
    typedef self pointer;
    typedef step_iterator_assignmenthelper<pointerType> reference;
	typedef CHAR * charptr;
protected:
    pointerType ptr;
    intT step;
public:
	step_iterator(pointerType ptr_, intT step_) : ptr(ptr_), step(step_) {};
	inline bool operator==(const self & b) const {
		return (ptr == b.ptr);
	};
	inline bool operator!=(const self & b) const {
		return (ptr != b.ptr);
	};
	inline value_type operator*() const {
		return *ptr;
	};
	inline reference operator*() {
		return reference( ptr );
	};
	inline self operator+(const intT b) const {
		return self( ptr + b * step, step );
	}
	inline self operator-(const intT b) const {
		return self( ptr - b* step, step );
	}
	inline self & operator++() {
		ptr += step;
		return *this;
	};
	inline self operator++(int) {
		self tmp(ptr, step) ;
		ptr += step;
		return tmp;
	};
	inline self & operator--() {
        ptr -= step;
		return *this;
	};
	inline self operator--(int) {
		self tmp(ptr, step);
		ptr -= step;
		return tmp;
	};
	inline self & operator+=(intT cnt) {
        ptr += step*cnt;
		return *this ;
	};
	inline self & operator-=(intT cnt) {
        ptr -= step*cnt;
		return *this;
	};
	inline value_type operator[]( intT i ) const {
		return ptr[ step*i ] ;
	};
	inline reference operator[](intT i) {
		return reference( ptr + step*i );
	};
	inline operator charptr() {
		return (charptr) ptr;
	}
};

template < typename pointerType, typename intT = int> class step_iterator2D : public step_iterator< pointerType, intT > {
public:
	typedef pointerType  value_type;    
    typedef step_iterator< pointerType, intT > parent ;
	step_iterator2D(pointerType ptr_, intT step_) : parent (ptr_,step_) {};
    inline value_type operator*() {
		return parent::ptr;
	};
    inline value_type operator[]( intT i ) {
		return parent::ptr+ parent::step*i ;
	};
};

template < typename pointerType > class complex_pointer : public std::iterator< std::random_access_iterator_tag , complex< typename std::iterator_traits< pointerType >::value_type > , typename std::iterator_traits< pointerType >::difference_type >  {
public:
    typedef typename std::iterator< std::random_access_iterator_tag , complex< typename std::iterator_traits< pointerType >::value_type > , typename std::iterator_traits< pointerType >::difference_type > parent; 
    typedef complex_pointer< pointerType > self;
    typedef typename std::iterator_traits< pointerType >::difference_type difference_type;
    typedef difference_type intT;
    typedef typename parent::value_type value_type;
    typedef typename parent::reference reference;
private:
    pointerType re;
    difference_type offset2Im;
    complex_pointer( pointerType re_, difference_type offset2Im_) : re(re_), offset2Im(offset2Im_) {};
public:
    complex_pointer( pointerType re_, pointerType im) : re(re_), offset2Im(im-re) {};
    
    inline bool operator==(const self & b) const {
		return (re == b.re) && (offset2Im==b.offset2Im);
	};
	inline bool operator!=(const self & b) const {
		return (re != b.re) || (offset2Im!=b.offset2Im);
	};
	inline value_type operator*() const {
		return value_type( *re, *( re+offset2Im ) );
	};
	/*inline reference operator*() {
		return reference( self );
	};*/
	inline self operator+(const intT b) const {
		return self( re + b , offset2Im );
	};
	inline self operator-(const intT b) const {
		return self( re - b, offset2Im );
	};
	inline self & operator++() {
		++re;
		return this;
	};
	inline self operator++(int) {
		self tmp( re, offset2Im) ;
		++re;
		return tmp;
	};
	inline self & operator--() {
        --re;
		return this;
	};
	inline self operator--(int) {
		self tmp(re, offset2Im);
		--re;
		return tmp;
	};
	inline self & operator+=(intT cnt) {
        re += cnt;
		return this ;
	};
	inline self & operator-=(intT cnt) {
        re -= cnt;
		return this;
	};
	inline value_type operator[]( intT i ) const {
		return value_type( re[ i ] , (re + offset2Im)[i]);
	};
	inline reference operator[](intT i) {
		return reference( self(re + i, offset2Im) );
	};
};



#endif