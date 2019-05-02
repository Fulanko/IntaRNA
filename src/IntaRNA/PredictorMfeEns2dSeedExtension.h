
#ifndef INTARNA_PREDICTORMFEENS2DSEEDEXTENSION_H_
#define INTARNA_PREDICTORMFEENS2DSEEDEXTENSION_H_

#include "IntaRNA/PredictorMfeEns.h"
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
class PredictorMfeEns2dSeedExtension: public PredictorMfeEns {

protected:

	//! matrix type to hold the partition functions for interaction site starts
	typedef boost::numeric::ublas::matrix<Z_type> Z2dMatrix;

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
	using PredictorMfeEns::energy;

	//! access to the output handler of the super class
	using PredictorMfeEns::output;

	//! partition function of all interaction hybrids that start on the left side of the seed including E_init
	Z2dMatrix hybridZ_left;

	//! the seed handler (with idx offset)
	SeedHandlerIdxOffset seedHandler;

	//! partition function of all interaction hybrids that start on the right side of the seed excluding E_init
	Z2dMatrix hybridZ_right;

	//
	struct ZPartition {
		size_t i1;
		size_t j1;
		size_t i2;
		size_t j2;
		Z_type partZ;
	};
	std::unordered_map<size_t, ZPartition> Z_partitions;

protected:

	/**
	 * Updates the overall hybridization partition function
	 * and the mfe information.
	 *
	 * Note: if called multiple time for the same boundaries then the
	 * reported partition functions have to represent disjoint interaction sets!
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1
	 * @param partFunct the partition function of the interaction
	 * @param isHybridZ whether or not the given Z is only covering
	 *        hybridization energy terms (init+loops) or the total
	 *        interaction energy
	 */
	virtual
	void
	updateZ( const size_t i1, const size_t j1
				, const size_t i2, const size_t j2
				, const Z_type partFunct
				, const bool isHybridZ );

	/**
	 * Computes all entries of the hybridE matrix for interactions ending in
	 * p=j1 and q=j2 and report all valid interactions to updateOptima()
	 *
	 * @param j1 end of the interaction within seq 1
	 * @param j2 end of the interaction within seq 2
	 * @param outConstraint constrains the interactions reported to the output handler
	 *
	 */
	virtual
	void
	fillHybridZ_left( const size_t j1, const size_t j2
				, const OutputConstraint & outConstraint );

	/**
	 * Computes all entries of the hybridE matrix for interactions starting in
	 * i1 and i2 and report all valid interactions to updateOptima()
	 *
	 * Note: (i1,i2) have to be complementary (right-most base pair of seed)
	 *
	 * @param i1 end of the interaction within seq 1
	 * @param i2 end of the interaction within seq 2
	 * @param outConstraint constrains the interactions reported to the output handler
	 *
	 */
	virtual
	void
	fillHybridZ_right( const size_t i1, const size_t i2
				, const OutputConstraint & outConstraint );

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

	/**
	 * Returns the hybridization energy of the non overlapping part of seeds si and sj
	 *
	 * @param si1 the index of seed1 in the first sequence
	 * @param si2 the index of seed1 in the second sequence
	 * @param sj1 the index of seed2 in the first sequence
	 * @param sj2 the index of seed2 in the second sequence
	 */
	virtual
	E_type
	getNonOverlappingEnergy( const size_t si1, const size_t si2, const size_t sj1, const size_t sj2 );

	// debug function
	void
	printMatrix( const Z2dMatrix & matrix );

};

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

inline
void
PredictorMfeEns2dSeedExtension::
updateZ( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const Z_type partZ
		, const bool isHybridZ )
{
	// update overall hybridization partition function
	PredictorMfeEns::updateZ( i1, j1, i2, j2, partZ, isHybridZ );

	// store partial Z
	size_t maxLength = std::max(energy.getAccessibility1().getMaxLength(), energy.getAccessibility2().getMaxLength());
	size_t key = 0;
	key += i1;
	key += j1 * pow(maxLength, 1);
	key += i2 * pow(maxLength, 2);
	key += j2 * pow(maxLength, 3);
	// TODO: check if key overflow
	if ( Z_partitions.find(key) == Z_partitions.end() ) {
		// create new entry
		ZPartition zPartition;
		zPartition.i1 = i1;
		zPartition.j1 = j1;
		zPartition.i2 = i2;
		zPartition.j2 = j2;
		zPartition.partZ = partZ;
		Z_partitions[key] = zPartition;
	} else {
		// update entry
		ZPartition & zPartition = Z_partitions[key];
		zPartition.partZ += partZ;
	}
}

//////////////////////////////////////////////////////////////////////////

inline
void
PredictorMfeEns2dSeedExtension::
printMatrix( const Z2dMatrix & matrix )
{
	std::cout << "==============" << std::endl;
	for (int i = 0; i < matrix.size1(); i++) {
		std::cout << "| ";
		for (int j = 0; j < matrix.size2(); j++) {
			std::cout << matrix(i, j) << " | ";
		}
		std::cout << std::endl;
	}
	std::cout << "==============" << std::endl;
}

//////////////////////////////////////////////////////////////////////////


} // namespace

#endif /* INTARNA_PREDICTORMFEENS2DSEEDEXTENSION_H_ */
