
#include "IntaRNA/PredictorMfe2dHeuristicSeedExtension.h"

namespace IntaRNA {

//////////////////////////////////////////////////////////////////////////

PredictorMfe2dHeuristicSeedExtension::
PredictorMfe2dHeuristicSeedExtension(
		const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker
		, SeedHandler * seedHandlerInstance )
 :
	PredictorMfe2dSeedExtension(energy,output,predTracker,seedHandlerInstance)
{
	assert( seedHandler.getConstraint().getBasePairs() > 1 );
}

//////////////////////////////////////////////////////////////////////////

PredictorMfe2dHeuristicSeedExtension::
~PredictorMfe2dHeuristicSeedExtension()
{
}

//////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dHeuristicSeedExtension::
predict( const IndexRange & r1, const IndexRange & r2
		, const OutputConstraint & outConstraint )
{
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"predicting mfe interactions with seed in O(n^2) space and O(n^4) time..."; }
	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

	// suboptimal setup check
	if (outConstraint.reportMax>1 && outConstraint.reportOverlap != OutputConstraint::ReportOverlap::OVERLAP_BOTH) {
		throw std::runtime_error("PredictorMfe2dHeuristicSeedExtension : the enumeration of non-overlapping suboptimal interactions is not supported in this prediction mode");
	}

#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!(r1.isAscending() && r2.isAscending()) )
		throw std::runtime_error("PredictorMfe2d::predict("+toString(r1)+","+toString(r2)+") is not sane");
#endif

	// setup index offset
	energy.setOffset1(r1.from);
	energy.setOffset2(r2.from);
	seedHandler.setOffset1(r1.from);
	seedHandler.setOffset2(r2.from);

	const size_t interaction_size1 = std::min( energy.size1()
			, (r1.to==RnaSequence::lastPos?energy.size1()-1:r1.to)-r1.from+1 );
	const size_t interaction_size2 = std::min( energy.size2()
			, (r2.to==RnaSequence::lastPos?energy.size2()-1:r2.to)-r2.from+1 );

	// compute seed interactions for whole range
	// and check if any seed possible
	if (seedHandler.fillSeed( 0, interaction_size1-1, 0, interaction_size2-1 ) == 0) {
		// trigger empty interaction reporting
		initOptima(outConstraint);
		reportOptima(outConstraint);
		// stop computation
		return;
	}

	// initialize mfe interaction for updates
	initOptima( outConstraint );

	for (size_t si1 = 0; si1 <= interaction_size1-seedHandler.getConstraint().getBasePairs(); si1++) {
	for (size_t si2 = 0; si2 <= interaction_size2-seedHandler.getConstraint().getBasePairs(); si2++) {
		E_type seedE = seedHandler.getSeedE(si1, si2);
		// check if valid left seed base pair
		if (E_isINF(seedE))
		  continue;

		const size_t sl1 = seedHandler.getSeedLength1(si1, si2)-1;
		const size_t sl2 = seedHandler.getSeedLength2(si1, si2)-1;
		const size_t sj1 = si1+sl1;
		const size_t sj2 = si2+sl2;
		const size_t maxMatrixLen1 = energy.getAccessibility1().getMaxLength()-sl1;
		const size_t maxMatrixLen2 = energy.getAccessibility2().getMaxLength()-sl2;
		// check if seed fits into interaction range
		if (sj1 > interaction_size1 || sj2 > interaction_size2)
			continue;

		// init optimal right boundaries
		j1opt = sj1;
		j2opt = sj2;

		// update mfe for seed only
		PredictorMfe2d::updateOptima( si1,sj1,si2,sj2, energy.getE_init() + seedHandler.getSeedE(si1, si2), true );

		// ER
		hybridE_right.resize( std::min(interaction_size1-sj1, maxMatrixLen1), std::min(interaction_size2-sj2, maxMatrixLen2) );
		fillHybridE_right(sj1, sj2, outConstraint, si1, si2);

		// EL
		hybridE_left.resize( std::min(si1+1, maxMatrixLen1), std::min(si2+1, maxMatrixLen2) );
		fillHybridE_left(si1, si2, outConstraint);


	} // si2
	} // si1

	// report mfe interaction
	reportOptima( outConstraint );
}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dHeuristicSeedExtension::
fillHybridE_right( const size_t i1, const size_t i2
			, const OutputConstraint & outConstraint
			, const size_t si1, const size_t si2 )
{

	// global vars to avoid reallocation
	size_t j1,j2,k1,k2;
	//////////  FIRST ROUND : COMPUTE HYBRIDIZATION ENERGIES ONLY  ////////////

	// current minimal value
	E_type curMinE = E_INF;
	// iterate over all window starts j1 (seq1) and j2 (seq2)
	for (j1=i1; j1-i1 < hybridE_right.size1(); j1++) {
		// screen for right boundaries in seq2
		for (j2=i2; j2-i2 < hybridE_right.size2(); j2++) {

			// init current cell (0 if just left (i1,i2) base pair)
			hybridE_right(j1-i1,j2-i2) = i1==j1 && i2==j2 ? 0 : E_INF;

			// check if complementary
			if( i1<j1 && i2<j2 && energy.areComplementary(j1,j2) ) {
				curMinE = E_INF;

				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				for (k1=j1; k1-- > i1; ) {
					// ensure maximal loop length
					if (j1-k1 > energy.getMaxInternalLoopSize1()+1) break;
				for (k2=j2; k2-- > i2; ) {
					// ensure maximal loop length
					if (j2-k2 > energy.getMaxInternalLoopSize2()+1) break;
					// check if (k1,k2) are valid left boundary
					if ( E_isNotINF( hybridE_right(k1-i1,k2-i2) ) ) {
						curMinE = std::min( curMinE,
								(energy.getE_interLeft(k1,j1,k2,j2)
										+ hybridE_right(k1-i1,k2-i2) )
								);
					}
				}
				}

				// store value
				hybridE_right(j1-i1,j2-i2) = curMinE;
				// update mfe if needed
				updateOptima( si1,i1+j1,si2,i2+j2, hybridE_right(j1-i1,j2-i2) + energy.getE_init() + seedHandler.getSeedE(si1, si2), true, si1, si2 );
			}
		}
	}

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dHeuristicSeedExtension::
fillHybridE_left( const size_t i1, const size_t i2
			, const OutputConstraint & outConstraint )
{

	// global vars to avoid reallocation
	size_t j1,j2,k1,k2;
	//////////  FIRST ROUND : COMPUTE HYBRIDIZATION ENERGIES ONLY  ////////////

	// current minimal value
	E_type curMinE = E_INF;
	// iterate over all window starts j1 (seq1) and j2 (seq2)
	for (j1=i1; i1-j1 < hybridE_left.size1(); j1--) {
		// screen for right boundaries in seq2
		for (j2=i2; i2-j2 < hybridE_left.size2(); j2--) {
			// init current cell (0 if just left (i1,i2) base pair)
			hybridE_left(i1-j1,i2-j2) = i1==j1 && i2==j2 ? energy.getE_init() : E_INF;
			// check if complementary
			if( j1<i1 && j2<i2 && energy.areComplementary(j1,j2) ) {
				curMinE = E_INF;

				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				for (k1=j1; k1++ < i1; ) {
					// ensure maximal loop length
					if (k1-j1 > energy.getMaxInternalLoopSize1()+1) break;
				for (k2=j2; k2++ < i2; ) {
					// ensure maximal loop length
					if (k2-j2 > energy.getMaxInternalLoopSize2()+1) break;
					// check if (k1,k2) are valid left boundary
					if ( E_isNotINF( hybridE_left(i1-k1,i2-k2) ) ) {
						curMinE = std::min( curMinE,
								(energy.getE_interLeft(j1,k1,j2,k2)
										+ hybridE_left(i1-k1,i2-k2) )
								);
					}
				} // k2
			  } // k1

				// store value
				hybridE_left(i1-j1,i2-j2) = curMinE;

				// update mfe if needed
				const size_t sl1 = seedHandler.getSeedLength1(i1, i2)-1;
				const size_t sl2 = seedHandler.getSeedLength2(i1, i2)-1;
				const size_t sj1 = i1+sl1;
				const size_t sj2 = i2+sl2;
				PredictorMfe2d::updateOptima( i1-std::min(i1, hybridE_left.size1()),j1opt,i2-std::min(i2, hybridE_left.size2()),j2opt, hybridE_right(j1opt-sj1, j2opt-sj2) + hybridE_left(i1-j1,i2-j2) + seedHandler.getSeedE(j1, j2), true );

			}
		}
	}

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dHeuristicSeedExtension::
getNextBest( Interaction & curBest )
{
	INTARNA_NOT_IMPLEMENTED("PredictorMfe2dHeuristicSeedExtension::getNextBest() not implemented yet");
}

//////////////////////////////////////////////////////////////////////////


} // namespace
