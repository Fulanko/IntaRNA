
#ifndef INTARNA_PREDICTIONTRACKERBASEPAIRPROB_H_
#define INTARNA_PREDICTIONTRACKERBASEPAIRPROB_H_

#include "IntaRNA/PredictionTracker.h"
#include "IntaRNA/InteractionEnergy.h"
#include "IntaRNA/Interaction.h"
#include "IntaRNA/PredictorMfeEns.h"

#include <iostream>

#include <boost/algorithm/string.hpp>

namespace IntaRNA {

/**
 * Collects partition function parts Z(i,j), computes base-pair probabilities
 * and prints them as a dot plot
 *
 * The information is written to stream on destruction.
 */
class PredictionTrackerBasePairProb: public PredictionTracker
{

public:

	/**
	 * Constructs a PredictionTracker that collects probability information
	 * for a set of interaction spots by computing the Boltzmann probability
	 * of any interaction enclosing the positions.
	 *
	 * @param energy the energy function used for energy calculation
	 * @param spots the interaction spots to track for their probabilities
	 * @param outStreamName the stream name where the probability data
	 *        is to be written to. use STDOUT/STDERR for the respective stream.
	 *        Otherwise, an according file is created
	 */
	PredictionTrackerBasePairProb(
				const InteractionEnergy & energy
				, const std::string & outStreamName
			);

	/**
	 * Constructs a PredictionTracker that collects probability information
	 * for a set of interaction spots by computing the Boltzmann probability
	 * of any interaction enclosing the positions.
	 *
	 * @param energy the energy function used for energy calculation
	 * @param spots the interaction spots to track for their probabilities
	 * @param outStream the stream where the probability data
	 *        is to be written to.
	 */
	PredictionTrackerBasePairProb(
				const InteractionEnergy & energy
				, std::ostream & outStream
			);


	/**
	 * destruction: write the probabilities to stream.
	 */
	virtual ~PredictionTrackerBasePairProb();


	/**
	 * Updates the probability information for each Predictor.updateOptima() call.
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1
	 * @param energy the overall energy of the interaction site
	 */
	virtual
	void
	updateOptimumCalled( const size_t i1, const size_t j1
						, const size_t i2, const size_t j2
						, const E_type energy
						);

	/**
	 * Updates the probability information.
	 *
	 * @param predictor the predictor providing the probability information
	 */
	virtual
	void
	updateZ( PredictorMfeEns *predictor ) override;

protected:

	struct StructureProb {
		size_t i1;
		size_t j1;
		size_t i2;
		size_t j2;
		float prob;
	};

	//! energy handler used for predictions
	const InteractionEnergy & energy;

	//! the stream to write the probabilities to
	std::ostream * outStream;

	//! whether or not outStream is to be deleted on destruction
	const bool deleteOutStream;

	//! overall partition function of all reported interactions
	Z_type overallZ;

	//! map storing structure probabilities
	std::unordered_map<size_t, StructureProb> structureProbs;

	//! pointer to mfeEnsPredictor
	PredictorMfeEns* mfeEnsPredictor;

	/**
	 * Generates key for storing values in map
	 */
	virtual
	size_t
	generateMapKey( const size_t i1, const size_t j1
						, const size_t i2, const size_t j2 ) const;

	/**
	 * computes probability of a sub structure (i1, j1, i2, j2) given a
	 * structure (k1, l1, k2, l2) with probability structProb
	 */
	virtual
	float
	computeSubStructureProb( const size_t i1, const size_t j1
						, const size_t i2, const size_t j2
						, const size_t k1, const size_t l1
						, const size_t k2, const size_t l2, const float structProb );

	/**
	 * Recursively computes structure probabilities of sub structures given
	 * an initial structure (k1, l1, k2, l2) with probability structProb
	 */
	virtual
	float
	computeStructureProbRecursive( const size_t k1, const size_t l1
						, const size_t k2, const size_t l2, const float structProb );

};

//////////////////////////////////////////////////////////////////////

} // namespace

#endif /* INTARNA_PREDICTIONTRACKERBASEPAIRPROB_H_ */
