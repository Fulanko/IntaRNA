#ifndef INTARNA_PREDICTORMFESEEDONLY_H_
#define INTARNA_PREDICTORMFESEEDONLY_H_

#include "IntaRNA/PredictorMfe.h"
#include "IntaRNA/SeedHandlerIdxOffset.h"

namespace IntaRNA {

/**
 * Prediction of seed interaction only.
 *
 * @author Martin Mann
 *
 */
class PredictorMfeSeedOnly: public PredictorMfe {

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
	PredictorMfeSeedOnly(
			const InteractionEnergy & energy
			, OutputHandler & output
			, PredictionTracker * predTracker
			, SeedHandler * seedHandler );


	/**
	 * data cleanup
	 */
	virtual ~PredictorMfeSeedOnly();


	/**
	 * Computes the mfe seed for the given sequence ranges (i1-j1) in the first
	 * sequence and (i2-j2) in the second sequence and reports it to the output
	 * handler.
	 *
	 * @param r1 the index range of the first sequence interacting with r2
	 * @param r2 the index range of the second sequence interacting with r1
	 *
	 */
	virtual
	void
	predict( const IndexRange & r1 = IndexRange(0,RnaSequence::lastPos)
			, const IndexRange & r2 = IndexRange(0,RnaSequence::lastPos) );


protected:


	//! access to the interaction energy handler of the super class
	using PredictorMfe::energy;

	//! access to the output handler of the super class
	using PredictorMfe::output;

	//! access to the output handler of the super class
	using PredictorMfe::reportedInteractions;

	//! the seed handler (with idx offset)
	SeedHandlerIdxOffset seedHandler;

	//! last position of a seed of interest (in index-range of seedHandler)
	size_t seedLastPos1;
	//! last position of a seed of interest (in index-range of seedHandler)
	size_t seedLastPos2;

protected:

	/**
	 * Fills a given interaction (boundaries given) with the according
	 * hybridizing base pairs using hybridE_seed.
	 * @param interaction IN/OUT the interaction to fill
	 */
	virtual
	void
	traceBack( Interaction & interaction );

	/**
	 * Identifies the next best interaction with an energy equal to or higher
	 * than the given interaction. The new interaction will not overlap any
	 * index range stored in reportedInteractions.
	 *
	 * @param curBest the interaction to be replaced
	 */
	virtual
	void
	getNextBest( Interaction & curBest );

};

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

} // namespace

#endif /* INTARNA_PREDICTORMFESEEDONLY_H_ */
