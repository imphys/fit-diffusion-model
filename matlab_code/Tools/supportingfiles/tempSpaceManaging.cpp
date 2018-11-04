#ifndef TEMPSPACEMANAGING
#define TEMPSPACEMANAGING
/*
 * This (tiny) helper file provides a simple class for temporary space managing.
 * Instances of tempspace should typically be passed by value, as this automagically 
 * clears the tempspace used by subroutines when they exit.
 * Pass by reference if you return variables that reside in tempspace. 
 *
 * - create:
 *   tempspace tmp( pointer_to_tempspace, number_of_bytes);
 *			you should allocate number_of_bytes bytes of memory, starting at pointer_to_tempspace
 *			Then this function creates the temp-space managing class tmp.
 *   getTempspaceSize tmp( reference_to_ptrdiff_t )
 *          Helper class that allows to compute the required temp space easily.
 *			Use in the same way as tempspace; e.g. in an overloaded function that is created to just 
 *		    compute the required space.
 *			the get method returns void. 
 * - get some space:
 *   T * var = tmp.get<T>( number_of_elements );
 *          allocates a pointer to an array of type T having number_of_elements elements and proper alignment.
 *          Errors if not sufficient space is available.
 * - Split tempspace in parts for 'OMP' multiprocessing:
 *   tmp_in_paralell = tmp.getParalellPart();
 *          Divides the remaining temporary space evenly over the threads.
 *			Call from each thread separately.
 *          Make sure you define 'USEOMP' before including this file.
 * - Check guards:
 *    tmp.checkGuard( var , number_of_elements );  
 *          if ADDGUARDS is defined, checks if the guards around var are still present, 
 *			and errors if not (=> out of bounds writing). If ADDGUARDS is not defined it does nothing. 
 *          var should be created by var = tmp.get<T>( number_of_elements );
 *			(NOTE: tempmanaging does not store info about the variables allocated, 
 *                 so this check cannot know with how many elements var was created.)
 * 
 * This file has two options (#define's):
 *    NOSPACELEFTCHECK  : (do-not-check-the-space-that-is-still-left)
 *						  if defined, the amount of space left is not checked (and invalid pointers are returned when to much is requested). 
 *						  Avoids a simple '<0' test per 'get'
 *                        So should only be used when you are absolutely sure enough tempspace has been allocated and 'get' is used extremely often.
 *    ADDGUARDS         : Adds a 4 byte 'guard' before and after each array returned by 'get' (take that into account when allocating pointer_to_tempspace)
 *                        call tmp.checkGuard( var , number_of_elements );
 *					      to check if the guard is still in place. 
 *   
 */
#ifdef _MSC_VER
	// Compiling with MS Visual studio compiler:
#define NOTALIASED __declspec(restrict)
#else
	// Compiling with other compiler: TODO: find way to declare output pointer is not aliased.
#define NOTALIASED 
#endif

#ifdef ADDGUARDS
typedef int guardT;
const guardT guardnumber = 0xA9876543 ; 
#pragma message("Adding guards in temporary space. Please call checkGuard at the end of any function that did use a part of the temporary space block.")
#endif

class tempspace {
private:
	char * baseptr;
	ptrdiff_t bytesleft;
public:
	tempspace( char * baseptr_, size_t numbytes) {
		baseptr = baseptr_;
		bytesleft = numbytes;
	};
	template <typename T, int align > NOTALIASED T * get(size_t numel) {
		ptrdiff_t numbytesrequested = numel * sizeof(T);
#ifdef ADDGUARDS
		numbytesrequested +=sizeof(guardT);
#endif
		T * ret = (T*) ( (ptrdiff_t) (baseptr + bytesleft - numbytesrequested) & (- (ptrdiff_t) align ) );
		bytesleft = ((char*) ret) - baseptr;
#ifdef ADDGUARDS
		bytesleft -= sizeof(guardT) ;
		*(((guardT *) ret)-1) = guardnumber; // store guard just before ret.
		*((guardT *) (ret + numel)) = guardnumber; // store guard just after ret.
#endif
#ifndef NOSPACELEFTCHECK
		if (bytesleft<0) {
			mexErrMsgTxt("Insufficient temporary space allocated.");
		}
#endif
		return ret;
	};
	template <typename T> inline NOTALIASED T * get( size_t numel) {
		return get<T,__alignof( T )>( numel );
	};

	tempspace getParalellPart() {
#ifdef USEOMP
		int thread_num = omp_get_thread_num();
		int num_threads = omp_get_num_threads();
		ptrdiff_t bytesperthread = bytesleft/num_threads;
		return tempspace( baseptr + thread_num * bytesperthread , bytesperthread );
#else
		return *this;
#endif
	}
	template <typename T> void checkGuard( T var , size_t numel ) {
#ifdef ADDGUARDS
		if ( (*(((guardT *) var)-1) != guardnumber) ||(*((guardT *) (var + numel)) != guardnumber)) {
			mexErrMsgTxt("A check determined that out-of-bounds writing has occured, please debug to find the offending instruction(s).");
		}
#endif
	}
};


class getTempspaceSize {
private:
	int nthreads;
	ptrdiff_t bytesused;
	ptrdiff_t & maxbytesused;
public:
	getTempspaceSize ( ptrdiff_t & maxbytesused_ ) : maxbytesused(maxbytesused_) , bytesused(maxbytesused_),nthreads(1) {};
	template <typename T> void get(size_t numel) {
		ptrdiff_t numbytesrequested = numel * sizeof(T) + (__alignof( T )-1); // worst case alignment requirement.
#ifdef ADDGUARDS
		numbytesrequested += 2 *sizeof(guardT);
#endif
		bytesused += numbytesrequested * nthreads;
		if (bytesused>maxbytesused) {
			maxbytesused = bytesused;
		}
	};
	tempspace getParalellPart() {
#ifdef USEOMP
		nthreads = omp_get_num_threads();
#else
		nthreads = 1;
#endif
	}
};


#endif