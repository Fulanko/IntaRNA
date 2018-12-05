
#ifndef INTARNA_GENERAL_H_
#define INTARNA_GENERAL_H_


//////////////  GENERAL CONFIGURE FLAGS  ////////////////

#ifdef HAVE_CONFIG_H
	#include "config.h"
#endif

#include "IntaRNA/intarna_config.h"

//! central flag whether or not debug mode is enabled
#define INTARNA_IN_DEBUG_MODE ((defined(_DEBUG)) || (!defined (NDEBUG)))

////////////////  CENTRAL LOGGING LIB  //////////////////

#include "easylogging++.h"

////////////////  GARBAGE COLLECTION  ///////////////////

#define  INTARNA_CLEANUP(pointer) if (pointer != NULL) {delete pointer; pointer=NULL;}

////////////////  ABORT NON-IMPLEMENTED  ////////////////

#include <stdexcept>

	#define INTARNA_NOT_IMPLEMENTED(message)	\
		throw std::runtime_error( \
					std::string("\nSTOP DUE TO MISSING IMPLEMENTATION : ") \
					+ message \
				);

////////////////  ABORT ON NULL POINTER  ////////////////

#include <stdexcept>

	#define INTARNA_CHECK_NOT_NULL(pointer,message) if (pointer == NULL) {	\
		throw std::runtime_error( \
					std::string("\nSTOP DUE TO NULL POINTER : ") \
					+ message \
				); }


////////////////  STRING GENERATION  ////////////////////

#include <boost/lexical_cast.hpp>

	#ifdef toString
		#error toString already defined
	#endif
	#define toString( x ) boost::lexical_cast<std::string>(x)


////////////////  NUMBER STRING LENGTH  ////////////////////

#include <cmath>

	#ifdef numStringLength
		#error numStringLength already defined
	#endif
	#define numStringLength( x ) ( (x==0) ? 1 : (1+std::floor(std::log10(x))) )


////////////////  GLOBAL TYPEDEFS  //////////////////////

#include <cmath>
#include <limits>

namespace IntaRNA {

	//! type for energy values in kcal/mol units for in-/output only
	typedef float E_kcal_type;

	//! type for energy values (energy + accessibility [ED]) (internally)
	typedef int E_type;
	const E_type E_INF = (std::numeric_limits<E_type>::max() / 2) + 1;
	const E_type E_MAX = E_INF / 2;

	//! type for probabilities, RT and Boltzmann values
	typedef float Z_type;
	const Z_type Z_INF = std::numeric_limits<Z_type>::infinity();

	//! type for temperature values (has to be compatible with Z_type for computations)
	typedef Z_type T_type;

} // namespace

#ifdef E_2_Ekcal
  #error E_2_Ekcal already defined
#endif
  //! convert internal energy type to energy value in kcal/mol units
#define E_2_Ekcal( e ) ( static_cast<E_kcal_type>(e) / 100.0 )

#ifdef Ekcal_2_E
  #error Ekcal_2_E already defined
#endif
  //! convert energy in kcal/mol units to internal energy type
#define Ekcal_2_E( e ) ( static_cast<E_type>(e * 100) )

#ifdef E_2_Z
  #error E_2_Z already defined
#endif
  //! convert E_type to Z_type
#define E_2_Z( e ) ( static_cast<Z_type>(e) / 100.0 )

#ifdef Z_2_E
  #error Z_2_E already defined
#endif
  //! convert Z_type to E_type
#define Z_2_E( e ) ( static_cast<E_type>(e * 100) )

#ifdef E_precisionEpsilon
	#error E_precisionEpsilon already defined
#endif
	//! the delta difference range to consider two energies equivalent
	//! only needed for floating point energy types:
    //! for these use sqrt(representable positive value closest to zero)
#define E_precisionEpsilon 0

#ifdef E_equal
	#error E_equal already defined
#endif
	//! check if two energies are equal
#define E_equal( e1, e2 ) ( e1 == e2 )

#ifdef E_isNotINF
	#error E_isNotINF already defined
#endif
	//! check if a given energy is NOT set to E_INF
#define E_isNotINF( e ) ( E_INF > e )

#ifdef E_isINF
	#error E_isINF already defined
#endif
	//! check if a given energy is set to E_INF
#define E_isINF( e ) ( E_INF <= e )



#ifdef Z_precisionEpsilon
	#error Z_precisionEpsilon already defined
#endif
	//! the delta difference range to consider two energies equivalent
    //! using sqrt(representable positive value closest to zero)
#define Z_precisionEpsilon std::sqrt(std::numeric_limits<Z_type>::min())

#ifdef Z_equal
	#error Z_equal already defined
#endif
	//! check if two energies are equal according to some epsilon
#define Z_equal( e1, e2 ) ( std::abs((e1)-(e2)) < Z_precisionEpsilon)
// another option from http://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
//#define Z_equal_ULP 2
//#define Z_equal( e1, e2 ) ( \
//	/* the machine epsilon has to be scaled to the magnitude of the values used */ \
//	/* and multiplied by the desired precision in ULPs (units in the last place) */ \
//	std::abs(e1-e2) < std::numeric_limits<T>::epsilon() * std::abs(e1+e2) * Z_equal_ULP \
//	/* unless the result is subnormal */ \
//	|| std::abs(e1-e2) < std::numeric_limits<T>::min() \
//)

#ifdef Z_isNotINF
	#error Z_isNotINF already defined
#endif
	//! check if a given energy is NOT set to Z_INF
#define Z_isNotINF( e ) ( std::numeric_limits<Z_type>::max() >= e )

#ifdef Z_isINF
	#error Z_isINF already defined
#endif
	//! check if a given energy is set to Z_INF
#define Z_isINF( e ) (  std::numeric_limits<Z_type>::max() < e )


////////////////  UTILITY FUNCTION  /////////////////////

#include <iostream>
#include <string>

namespace IntaRNA {

/**
 * Returns an open output stream for the given output name, i.e.
 *
 * - & std::cout : if outName == STDOUT
 * - & std::cerr : if outName == STDERR
 * - new std::fstream( outName ) : else if outName non-empty
 *
 * @param outName the name of the output to open. use STDOUT/STDERR for the
 *        respective output stream or otherwise a filename to be created.
 *
 * @return the open output stream, or NULL in error case
 *
 */
std::ostream *
newOutputStream( const std::string & outName );


/**
 * Deletes an output stream pointer returned by newOutputStream().
 * If it is a file stream, it is closed before deletion.
 * If it is pointing to std::cout or std::cerr, nothing is done.
 *
 * @param outStream the pointer to the output stream to delete
 */
void
deleteOutputStream( std::ostream * outStream );

} // namespace

#endif /* GENERAL_H_ */
