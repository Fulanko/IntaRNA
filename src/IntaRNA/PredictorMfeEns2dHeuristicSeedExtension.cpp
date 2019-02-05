#include "IntaRNA/PredictorMfeEns2dHeuristicSeedExtension.h"

namespace IntaRNA {

//////////////////////////////////////////////////////////////////////////

PredictorMfeEns2dHeuristicSeedExtension::
PredictorMfeEns2dHeuristicSeedExtension(
		const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker
		, SeedHandler * seedHandlerInstance )
 :
	PredictorMfeEns2dSeedExtension(energy,output,predTracker,seedHandlerInstance)
, E_right_opt(E_INF)
, j1opt(0)
, j2opt(0)
{
	assert( seedHandler.getConstraint().getBasePairs() > 1 );
}

//////////////////////////////////////////////////////////////////////////

PredictorMfeEns2dHeuristicSeedExtension::
~PredictorMfeEns2dHeuristicSeedExtension()
{
}

//////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dHeuristicSeedExtension::
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
		throw std::runtime_error("PredictorMfeEns2dHeuristicSeedExtension : the enumeration of non-overlapping suboptimal interactions is not supported in this prediction mode");
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
		// check if seed fits into interaction range
		if (sj1 > interaction_size1 || sj2 > interaction_size2)
			continue;

		// init optimal right boundaries
		E_right_opt = E_INF;
		j1opt = sj1;
		j2opt = sj2;

		// ER (min size*() == 1 == right-most bp of seed)
		hybridZ_right.resize( std::min(r1.to - sj1 + 1, interaction_size1-sl1+1)
							, std::min(r2.to - sj2 + 1, interaction_size2-sl2+1) );
		fillHybridZ_right(sj1, sj2, outConstraint, sj1+hybridZ_right.size1()-1, sj2+hybridZ_right.size2()-1, si1, si2);

		// EL
		// has to be resized based on interaction_size* and sl* (index shift needed) and ranges
		hybridZ_left.resize( std::min( si1 - r1.from + 1, interaction_size1-sl1+1)
						,  std::min( si2 - r2.from + 1, interaction_size2-sl2+1) );
		fillHybridZ_left(si1, si2, outConstraint, 0, 0);

//		// ER
//		hybridZ_right.resize( interaction_size1-sj1, interaction_size2-sj2);
//		fillHybridZ_right(sj1, sj2, outConstraint, interaction_size1-1, interaction_size2-1, si1, si2);
//
//		// EL
//		hybridZ_left.resize( si1+1, si2+1 );
//		fillHybridZ_left(si1, si2, outConstraint, 0, 0);

	} // si2
	} // si1

	// report mfe interaction
	reportOptima( outConstraint );

}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dHeuristicSeedExtension::
fillHybridZ_right( const size_t i1, const size_t i2
			, const OutputConstraint & outConstraint
			, const size_t j1max, const size_t j2max
			, const size_t si1, const size_t si2 )
{

	// global vars to avoid reallocation
	size_t j1,j2,k1,k2;

	// current minimal value
	const Z_type seedZ = energy.getBoltzmannWeight(seedHandler.getSeedE(si1, si2));
	const Z_type initZ = energy.getBoltzmannWeight(energy.getE_init());
	// iterate over all window starts j1 (seq1) and j2 (seq2)
	for (j1=i1; j1 <= j1max; j1++ ) {
		// screen for right boundaries in seq2
		for (j2=i2; j2 <= j2max; j2++ ) {

			// init current cell (0 if just left (i1,i2) base pair)
			hybridZ_right(j1-i1,j2-i2) = i1==j1 && i2==j2 ? energy.getBoltzmannWeight(0.0) : 0.0;

			// check if complementary
			if( i1<j1 && i2<j2 && energy.areComplementary(j1,j2) ) {

				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				for (k1=j1; k1-- > i1; ) {
					// ensure maximal loop length
					if (j1-k1 > energy.getMaxInternalLoopSize1()+1) break;
					for (k2=j2; k2-- > i2; ) {
						// ensure maximal loop length
						if (j2-k2 > energy.getMaxInternalLoopSize2()+1) break;
						// check if (k1,k2) is valid left boundary
						if ( ! Z_equal(hybridZ_right(k1-i1,k2-i2), 0.0) ) {
							// store value
							hybridZ_right(j1-i1,j2-i2) += energy.getBoltzmannWeight(energy.getE_interLeft(k1,j1,k2,j2)) * hybridZ_right(k1-i1,k2-i2);
						}
					}
				}
			}

			// update optimum of right extension if needed
			updateOptRightZ( si1,j1,si2,j2, energy.getE(hybridZ_right(j1-i1,j2-i2) * initZ * seedZ) );

		}
	}

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dHeuristicSeedExtension::
fillHybridZ_left( const size_t j1, const size_t j2
			, const OutputConstraint & outConstraint
			, const size_t i1init, const size_t i2init )
{

	// global vars to avoid reallocation
	size_t i1,i2,k1,k2;

	INTARNA_NOT_IMPLEMENTED("rewrite analogously to non-heuristic variant")
//
//	// current minimal value
//	Z_type curZ = 0.0;
//	// iterate over all window starts j1 (seq1) and j2 (seq2)
//	for (i1=j1+1; i1-- > i1init; ) {
//		// screen for left boundaries in seq2
//		for (i2=j2+1; i2-- > i2init; ) {
//
//			// init current cell (0 if just left (i1,i2) base pair)
//			hybridZ_left(i1,i2) = i1==j1 && i2==j2 ? energy.getBoltzmannWeight(energy.getE_init()) : 0.0;
//
//			// check if complementary
//			if( i1<j1 && i2<j2 ) {
//				curZ = 0.0;
//
//				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
//				for (k1=j1-1; k1 > i1; k1--) {
//					// ensure maximal loop length
//					if (k1-i1 > energy.getMaxInternalLoopSize1()+1) break;
//				for (k2=j2-1; k2 > i2; k2--) {
//					// ensure maximal loop length
//					if (k2-i2 > energy.getMaxInternalLoopSize2()+1) break;
//					// check if (k1,k2) are valid left boundary
//					if ( hybridZ_left(k1,k2) != 0.0 ) {
//					  curZ += energy.getBoltzmannWeight(energy.getE_interLeft(i1,k1,i2,k2)) * hybridZ_left(k1,k2);
//					}
//				}
//				}
//
//				// store value
//				hybridZ_left(i1,i2) = curZ;
//			}
//
//			// update mfe if needed
//			const size_t sl1 = seedHandler.getSeedLength1(j1, j2)-1;
//			const size_t sl2 = seedHandler.getSeedLength2(j1, j2)-1;
//			const size_t sj1 = j1+sl1;
//			const size_t sj2 = j2+sl2;
//			PredictorMfe2d::updateOptima( i1,j1opt,i2,j2opt, energy.getE(hybridZ_right(j1opt-sj1, j2opt-sj2) * hybridZ_left(i1,i2) * energy.getBoltzmannWeight(seedHandler.getSeedE(j1, j2))), true );
//
//		}
//	}

}

////////////////////////////////////////////////////////////////////////////
//
//void
//PredictorMfeEns2dHeuristicSeedExtension::
//traceBack( Interaction & interaction, const OutputConstraint & outConstraint  )
//{
//	// check if something to trace
//	if (interaction.basePairs.size() < 2) {
//		return;
//	}
//
//#if INTARNA_IN_DEBUG_MODE
//	// sanity checks
//	if ( interaction.basePairs.size() != 2 ) {
//		throw std::runtime_error("PredictorMfeEns2dHeuristicSeedExtension::traceBack() : given interaction does not contain boundaries only");
//	}
//#endif
//
//	// check for single interaction
//	if (interaction.basePairs.at(0).first == interaction.basePairs.at(1).first) {
//		// delete second boundary (identical to first)
//		interaction.basePairs.resize(1);
//		// update done
//		return;
//	}
//
//#if INTARNA_IN_DEBUG_MODE
//	// sanity checks
//	if ( ! interaction.isValid() ) {
//		throw std::runtime_error("PredictorMfeEns2dHeuristicSeedExtension::traceBack() : given interaction not valid");
//	}
//#endif
//
//	// ensure sorting
//	interaction.sort();
//	// get indices in hybridE for boundary base pairs
//	size_t	i1 = energy.getIndex1(interaction.basePairs.at(0)),
//			j1 = energy.getIndex1(interaction.basePairs.at(1)),
//			i2 = energy.getIndex2(interaction.basePairs.at(0)),
//			j2 = energy.getIndex2(interaction.basePairs.at(1))
//			;
//
//#if INTARNA_IN_DEBUG_MODE
//	// check if intervals are larger enough to contain a seed
//	if (std::min(j1-i1,j2-i2)+1 < seedHandler.getConstraint().getBasePairs()) {
//		// no seed possible, abort computation
//		throw std::runtime_error("PredictorMfeEns2dHeuristicSeedExtension::traceBack() : given boundaries "+toString(interaction)+" can not hold a seed of "+toString(seedHandler.getConstraint().getBasePairs())+" base pairs");
//	}
//#endif
//
//PredictorMfeEns2dSeedExtension::traceBack(interaction, outConstraint);
//
//}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dHeuristicSeedExtension::
getNextBest( Interaction & curBest )
{
	INTARNA_NOT_IMPLEMENTED("PredictorMfeEns2dHeuristicSeedExtension::getNextBest() not implemented yet");
}

//////////////////////////////////////////////////////////////////////////


} // namespace
