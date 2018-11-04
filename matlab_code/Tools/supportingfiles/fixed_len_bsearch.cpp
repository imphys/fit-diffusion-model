#ifndef FIXEDLENBSEARCH
#define FIXEDLENBSEARCH
/* fixed_len_bsearch.cpp
 * A binary search with a fixed length of the array. 
 * the search tree is fully unroled, so no data dependend jumps are made. 
 * Only use when performance is critical and the length of the search array is known at compile time.
 *
 * Useage:
 *
 * bin = binsearch<N, iteratorT>( array, value ) 
 * locates the bin such that:
 *  array[bin] <= value  < array[bin+1]  if  bin>=0 && bin+1<N
 * and: bin=0   if value <  array[0]
 *      bin=N-1 if value >= array[N-1]
 * Array should be defined for elements 0...N-1 and sorted with respect to '<='
 * Always uses exactly ceil(log2(N)) comparisons.
 *
 * Although this is slightly more than the minimal required average of log2(N), 
 * the fact that this is completely unroled and does not contain jumps (when compiler optimizations are turned on)
 * makes up for that. See test_unroled_binary_search.m for algorithmic correctness tests.
 */


using std::iterator_traits;

template<int N, int minindx, int maxindx, bool lastloop, typename itT> class binsearchi;

#define MULWITHBOOL
template < int N, int minindx, int maxindx, typename itT> class binsearchi< N, minindx, maxindx, true, itT > {
public : static inline int get( itT iterator, typename iterator_traits<itT>::value_type value, int curindex) {
	typedef typename iterator_traits<itT>::value_type value_type;
	// last loop; stepup==0, stepdown==1
	const int stepdown = 1;
	const int stepup   = 0;
#ifdef MULWITHBOOL
	return curindex -stepdown+ ( iterator[curindex] <= value )*(stepup+stepdown);
#else
	return curindex + ( iterator[curindex] <= value ? stepup : -stepdown);
#endif
}} ;

template<int N, int minindx, int maxindx, typename itT> class binsearchi<N, minindx, maxindx, false, itT>{ 
public: static inline int get( itT iterator, typename iterator_traits<itT>::value_type value, int curindex) {
	typedef typename iterator_traits<itT>::value_type value_type;
	const int stepdown = (minindx)/2;
	const int stepup   = (N-maxindx)/2;
	const bool lastloop = minindx-stepdown<=1;
#ifdef MULWITHBOOL
	curindex += -stepdown+( iterator[curindex] <= value )*(stepup+stepdown);
#else
	curindex += ( iterator[curindex] <= value ? stepup : -stepdown);
#endif
	return binsearchi<N,minindx-stepdown, maxindx+stepup, lastloop, itT>::get( iterator , value, curindex);
}
} ;

#ifdef MULWITHBOOL
#undef MULWITHBOOL
#endif
template<int N, typename itT> int binsearch( itT iterator, typename iterator_traits<itT>::value_type value) {
	typedef typename iterator_traits<itT>::value_type value_type;
	if (N>1) {
		const int initindx = (N+1)/2; 
		return binsearchi<N, initindx, initindx, (N<=2) , itT>::get( iterator, value, initindx );
	} else {
		return N-1;
	}
}
#endif