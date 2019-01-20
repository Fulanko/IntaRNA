
#ifndef INTARNA_PREDICTORMFEENS2DSEEDEXTENSION_H_
#define INTARNA_PREDICTORMFEENS2DSEEDEXTENSION_H_

#include "IntaRNA/PredictorMfe2d.h"
#include "IntaRNA/SeedHandlerIdxOffset.h"

namespace IntaRNA {

/**
 * Implements seed-based space-efficient interaction prediction.
 *
 * Note, for each seed start (i1,i2) only the mfe seed is considered for the
 * overall interaction computation instead of considering all possible seeds
 * starting at (i1,i2).
 *
 * @author Martin Mann
 *
 */
class PredictorMfeEns2dSeedExtension: public PredictorMfe2d {

protected:

	typedef float Z_type;
	const Z_type Z_INF = std::numeric_limits<Z_type>::infinity();

	//! the delta difference range to consider two energies equivalent
  //! using sqrt(representable positive value closest to zero)
	#define Z_precisionEpsilon 0.0001

	//! check if two energies are equal according to some epsilon
	#define Z_equal( e1, e2 ) ( std::abs((e1)-(e2)) < Z_precisionEpsilon )

	//! check if a given energy is NOT set to Z_INF
	#define Z_isNotINF( e ) ( std::numeric_limits<Z_type>::max() >= e )

	//! check if a given energy is set to Z_INF
	#define Z_isINF( e ) (  std::numeric_limits<Z_type>::max() < e )

	//! matrix type to hold the mfe energies for interaction site starts
	typedef boost::numeric::ublas::matrix<Z_type> E2dMatrix;

public:


	/**
	 * Constructs a predictor and stores the energy and output handler
	 *
	 * @param energy the interaction energy handler
	 * @param output the output handler to report mfe interactions to
	 * @param predTracker the prediction tracker to be used or NULL if no
	 *         tracking is to be done; if non-NULL, the tracker gets deleted
	 *         on this->destruction.
	 * @param seedHandler the seed handler to be used
	 */
	PredictorMfeEns2dSeedExtension(
			const InteractionEnergy & energy
			, OutputHandler & output
			, PredictionTracker * predTracker
			, SeedHandler * seedHandler );


	/**
	 * data cleanup
	 */
	virtual ~PredictorMfeEns2dSeedExtension();


	/**
	 * Computes the mfe for the given sequence ranges (i1-j1) in the first
	 * sequence and (i2-j2) in the second sequence and reports it to the output
	 * handler.
	 *
	 * Each considered interaction contains a seed according to the seed handler
	 * constraints.
	 *
	 * @param r1 the index range of the first sequence interacting with r2
	 * @param r2 the index range of the second sequence interacting with r1
	 * @param outConstraint constrains the interactions reported to the output handler
	 *
	 */
	virtual
	void
	predict( const IndexRange & r1 = IndexRange(0,RnaSequence::lastPos)
			, const IndexRange & r2 = IndexRange(0,RnaSequence::lastPos)
			, const OutputConstraint & outConstraint = OutputConstraint() );


protected:


	//! access to the interaction energy handler of the super class
	using PredictorMfe2d::energy;

	//! access to the output handler of the super class
	using PredictorMfe2d::output;

	//! energy of all interaction hybrids that end in position p (seq1) and
	//! q (seq2) and do not necessarily contain a seed interaction
	E2dMatrix hybridE_pq;

	//! the seed handler (with idx offset)
	SeedHandlerIdxOffset seedHandler;

	//! energy of all interaction hybrids that start in position p (seq1) and
	//! q (seq2)
	E2dMatrix hybridE_right;

protected:

	/**
	 * does nothing but to ignore the calls from fillHybridE()
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1
	 * @param energy ignored
	 * @param isHybridE ignored
	 */
	virtual
	void
	updateOptima( const size_t i1, const size_t j1
			, const size_t i2, const size_t j2
			, const E_type energy
			, const bool isHybridE );

	/**
	 * Computes all entries of the hybridE matrix for interactions ending in
	 * p=j1 and q=j2 and report all valid interactions to updateOptima()
	 *
	 * @param j1 end of the interaction within seq 1
	 * @param j2 end of the interaction within seq 2
	 * @param outConstraint constrains the interactions reported to the output handler
	 * @param i1init smallest value for i1
	 * @param i2init smallest value for i2
	 *
	 */
	void
	fillHybridE( const size_t j1, const size_t j2
				, const OutputConstraint & outConstraint
				, const size_t i1init, const size_t i2init
				);

	/**
	 * Computes all entries of the hybridE matrix for interactions starting in
	 * p=i1 and q=i2 and report all valid interactions to updateOptima()
	 *
	 * @param i1 end of the interaction within seq 1
	 * @param i2 end of the interaction within seq 2
	 * @param outConstraint constrains the interactions reported to the output handler
	 * @param j1init largest value for j1
	 * @param j2init largest value for j2
	 *
	 */
	void
	fillHybridE_right( const size_t i1, const size_t i2
				, const OutputConstraint & outConstraint
				, const size_t j1init, const size_t j2init
				);

	/**
	 * Fills a given interaction (boundaries given) with the according
	 * hybridizing base pairs using hybridE_seed.
	 * @param interaction IN/OUT the interaction to fill
	 * @param outConstraint constrains the interactions reported to the output handler
	 */
	virtual
	void
	traceBack( Interaction & interaction, const OutputConstraint & outConstraint  );

	/**
	 * Identifies the next best interaction with an energy equal to or higher
	 * than the given interaction. The new interaction will not overlap any
	 * index range stored in reportedInteractions.
	 *
	 * NOTE: this is not possible for this predictor (unless a full recomputation
	 * of the matrices is done). Thus, calling this method raises an exception.
	 *
	 * @param curBest ignored (see method comment)
	 */
	virtual
	void
	getNextBest( Interaction & curBest );

	void
	printMatrix( const E2dMatrix & matrix );

};

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

inline
void
PredictorMfeEns2dSeedExtension::
updateOptima( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const E_type energy
		, const bool isHybridE )
{
	// do nothing and ignore calls from fillHybridE()
}

//////////////////////////////////////////////////////////////////////////

inline
void
PredictorMfeEns2dSeedExtension::
printMatrix( const E2dMatrix & matrix )
{
	for (int i = 0; i < matrix.size1(); i++) {
		std::cout << "| ";
		for (int j = 0; j < matrix.size2(); j++) {
			std::cout << matrix(i, j) << " | ";
		}
		std::cout << std::endl;
	}
}

} // namespace

#endif /* INTARNA_PREDICTORMFEENS2DSEEDEXTENSION_H_ */
