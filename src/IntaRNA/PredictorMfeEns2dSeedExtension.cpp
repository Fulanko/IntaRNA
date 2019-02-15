
#include "IntaRNA/PredictorMfeEns2dSeedExtension.h"

namespace IntaRNA {

//////////////////////////////////////////////////////////////////////////

PredictorMfeEns2dSeedExtension::
PredictorMfeEns2dSeedExtension(
		const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker
		, SeedHandler * seedHandlerInstance )
 :
	PredictorMfe2d(energy,output,predTracker)
	, seedHandler(seedHandlerInstance)
	, hybridZ_right( 0,0 )
{
	assert( seedHandler.getConstraint().getBasePairs() > 1 );
}

//////////////////////////////////////////////////////////////////////////

PredictorMfeEns2dSeedExtension::
~PredictorMfeEns2dSeedExtension()
{
}

//////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dSeedExtension::
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
		throw std::runtime_error("PredictorMfeEns2dSeedExtension : the enumeration of non-overlapping suboptimal interactions is not supported in this prediction mode");
	}

#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!(r1.isAscending() && r2.isAscending()) )
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::predict("+toString(r1)+","+toString(r2)+") is not sane");
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
		const Z_type seedZ = energy.getBoltzmannWeight( seedE );

		const size_t sl1 = seedHandler.getSeedLength1(si1, si2)-1;
		const size_t sl2 = seedHandler.getSeedLength2(si1, si2)-1;
		const size_t sj1 = si1+sl1;
		const size_t sj2 = si2+sl2;
		const size_t maxMatrixLen1 = energy.getAccessibility1().getMaxLength()-sl1;
		const size_t maxMatrixLen2 = energy.getAccessibility2().getMaxLength()-sl2;
		// check if seed fits into interaction range
		if (sj1 > interaction_size1 || sj2 > interaction_size2)
			continue;

		// EL
		// has to be resized based on interaction_size* and sl* (index shift needed) and ranges
		hybridZ_left.resize( std::min(si1+1, maxMatrixLen1), std::min(si2+1, maxMatrixLen2) );
		fillHybridZ_left(si1, si2, outConstraint);

		// ER (min size*() == 1 == right-most bp of seed)
		hybridZ_right.resize( std::min(interaction_size1-sj1, maxMatrixLen1), std::min(interaction_size2-sj2, maxMatrixLen2) );
		fillHybridZ_right(sj1, sj2, outConstraint);

		// update Optimum for all boundary combinations
		for (size_t i1 = 0; i1<hybridZ_left.size1(); i1++) {
			// ensure max interaction length in seq 1
			for (size_t j1 = 0; j1 < hybridZ_right.size1() ; j1++) {
				// check interaction length
				if (sj1+j1-si1+i1 > energy.getAccessibility1().getMaxLength()) continue;
				for (size_t i2 = 0; i2< hybridZ_left.size2(); i2++) {
					// check complementarity of boundary
					if ( Z_equal(hybridZ_left(i1,i2), 0.0) ) continue;
					// ensure max interaction length in seq 2
					for (size_t j2 = 0; j2 < hybridZ_right.size2() ; j2++) {
						// check interaction length
						if (sj2+j2-si2+i2 > energy.getAccessibility2().getMaxLength()) continue;
						// check complementarity of boundary
						if (Z_equal(hybridZ_right(j1,j2),0.0)) continue;
						// compute full partition function
						Z_type fullZ = seedZ * hybridZ_left(i1,i2) * hybridZ_right(j1,j2);
						// update ensemble mfe
						if (! Z_equal(fullZ, 0.0)) {
							PredictorMfe::updateOptima( si1-i1, sj1+j1, si2-i2, sj2+j2, energy.getE(fullZ), true );
						}
					} // dj2
				} // di2
			} // dj1
		} // di1

	} // si2
	} // si1

	// report mfe interaction
	reportOptima( outConstraint );

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dSeedExtension::
fillHybridZ_left( const size_t j1orig, const size_t j2orig
			, const OutputConstraint & outConstraint )
{
#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!energy.areComplementary(j1orig,j2orig) )
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::fillHybridZ_left("+toString(j1orig)+","+toString(j2orig)+",..) are not complementary");
#endif


	// NOTE: indices not in global indexing but in local hybrid_pq matrix
	const size_t j1 = hybridZ_left.size1()-1, j2 = hybridZ_left.size2()-1;
	// compute shift to global sequence indexing
	const size_t shift1 = j1orig-j1, shift2 = j2orig-j2;
	// global vars to avoid reallocation
	size_t i1,i2,k1,k2;

	// iterate over all window starts i1 (seq1) and i2 (seq2)
	for (i1=j1; j1-i1 < hybridZ_left.size1(); i1-- ) {
		for (i2=j2; j2-i2 < hybridZ_left.size2(); i2-- ) {

			// init current cell (0 if not just right-most (j1,j2) base pair)
			hybridZ_left(j1-i1,j2-i2) = i1==j1 && i2==j2 ? energy.getBoltzmannWeight(energy.getE_init()) : 0.0;

			// check if complementary (use global sequence indexing)
			if( i1<j1 && i2<j2 && energy.areComplementary(shift1+i1,shift2+i2) ) {

				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				for (k1=i1; k1++ < j1; ) {
					// ensure maximal loop length
					if (k1-i1 > energy.getMaxInternalLoopSize1()+1) break;
					for (k2=i2; k2++ < j2; ) {
						// ensure maximal loop length
						if (k2-i2 > energy.getMaxInternalLoopSize2()+1) break;
						// check if (k1,k2) are valid left boundary
						if ( ! Z_equal(hybridZ_left(j1-k1,j2-k2), 0.0) ) {
							hybridZ_left(j1-i1,j2-i2) += energy.getBoltzmannWeight(energy.getE_interLeft(shift1+i1,shift1+k1,shift2+i2,shift2+k2)) * hybridZ_left(j1-k1,j2-k2);
						}
					} // k2
				} // k1
			}
		}
	}

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dSeedExtension::
fillHybridZ_right( const size_t i1, const size_t i2
			, const OutputConstraint & outConstraint )
{
#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!energy.areComplementary(i1,i2) )
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::fillHybridZ_right("+toString(i1)+","+toString(i2)+",..) are not complementary");
#endif


	// global vars to avoid reallocation
	size_t j1,j2,k1,k2;

	// iterate over all window ends j1 (seq1) and j2 (seq2)
	for (j1=i1; j1-i1 < hybridZ_right.size1(); j1++ ) {
		for (j2=i2; j2-i2 < hybridZ_right.size2(); j2++ ) {

			// init partition function for current cell -> (i1,i2) are complementary per definition
			hybridZ_right(j1-i1,j2-i2) = i1==j1 && i2==j2 ? energy.getBoltzmannWeight(0.0) : 0.0;

			// check if complementary free base pair
			if( i1<j1 && i2<j2 && energy.areComplementary(j1,j2) ) {

				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				for (k1=j1; k1-- > i1; ) {
					// ensure maximal loop length
					if (j1-k1 > energy.getMaxInternalLoopSize1()+1) break;
					for (k2=j2; k2-- > i2; ) {
						// ensure maximal loop length
						if (j2-k2 > energy.getMaxInternalLoopSize2()+1) break;
						// check if (k1,k2) are valid left boundary
						if ( ! Z_equal(hybridZ_right(k1-i1,k2-i2), 0.0) ) {
							// update partition function
							hybridZ_right(j1-i1,j2-i2) += energy.getBoltzmannWeight(energy.getE_interLeft(k1,j1,k2,j2)) * hybridZ_right(k1-i1,k2-i2);
						}
					} // k2
				} // k1
			}
		}
	}

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dSeedExtension::
traceBack( Interaction & interaction, const OutputConstraint & outConstraint  )
{
	// check if something to trace
	if (interaction.basePairs.size() < 2) {
		return;
	}

#if INTARNA_IN_DEBUG_MODE
	// sanity checks
	if ( interaction.basePairs.size() != 2 ) {
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::traceBack() : given interaction does not contain boundaries only");
	}
#endif

	// check for single interaction
	if (interaction.basePairs.at(0).first == interaction.basePairs.at(1).first) {
		// delete second boundary (identical to first)
		interaction.basePairs.resize(1);
		// update done
		return;
	}

#if INTARNA_IN_DEBUG_MODE
	// sanity checks
	if ( ! interaction.isValid() ) {
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::traceBack() : given interaction not valid");
	}
#endif

	// ensure sorting
	interaction.sort();
	// get indices in hybridE for boundary base pairs
	size_t	i1 = energy.getIndex1(interaction.basePairs.at(0)),
			j1 = energy.getIndex1(interaction.basePairs.at(1)),
			i2 = energy.getIndex2(interaction.basePairs.at(0)),
			j2 = energy.getIndex2(interaction.basePairs.at(1))
			;

#if INTARNA_IN_DEBUG_MODE
	// check if intervals are larger enough to contain a seed
	if (std::min(j1-i1,j2-i2)+1 < seedHandler.getConstraint().getBasePairs()) {
		// no seed possible, abort computation
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::traceBack() : given boundaries "+toString(interaction)+" can not hold a seed of "+toString(seedHandler.getConstraint().getBasePairs())+" base pairs");
	}
#endif

	E_type fullE = interaction.energy;

	for (size_t si1 = i1; si1 <= j1+1-seedHandler.getConstraint().getBasePairs(); si1++) {
		for (size_t si2 = i2; si2 <= j2+1-seedHandler.getConstraint().getBasePairs(); si2++) {
			E_type seedE = seedHandler.getSeedE(si1, si2);
			if ( E_isNotINF( seedE ) ) {

				const size_t sl1 = seedHandler.getSeedLength1(si1, si2)-1;
				const size_t sl2 = seedHandler.getSeedLength2(si1, si2)-1;
				const size_t sj1 = si1+sl1;
				const size_t sj2 = si2+sl2;
				const size_t maxMatrixLen1 = energy.getAccessibility1().getMaxLength()-sl1;
				const size_t maxMatrixLen2 = energy.getAccessibility2().getMaxLength()-sl2;

				// compute auxiliary matrices for seed check
				hybridZ_left.resize( std::min(si1+1, maxMatrixLen1), std::min(si2+1, maxMatrixLen2) );
				fillHybridZ_left( si1, si2, outConstraint );
				hybridZ_right.resize( std::min(j1-sj1+1, maxMatrixLen1), std::min(j2-sj2+1, maxMatrixLen2) );
				fillHybridZ_right( sj1, sj2, outConstraint );

				// check if we found the right seed for the interaction energy
				if ( E_equal( fullE,
						(energy.getE(i1, j1, i2, j2, energy.getE(energy.getBoltzmannWeight(seedE) * hybridZ_left( si1-i1, si2-i2 ) * hybridZ_right( j1-sj1, j2-sj2 ))))))
				{
					// found seed -> traceback seed base pairs
					if (si1 > i1 && si2 > i2) {
						interaction.basePairs.push_back( energy.getBasePair(si1,si2) );
					}
					// store seed information
					interaction.setSeedRange(
									energy.getBasePair(si1,si2),
									energy.getBasePair(sj1,sj2),
									energy.getE(si1,sj1,si2,sj2,seedE)+energy.getE_init());
					// trace seed base pairs
					seedHandler.traceBackSeed( interaction, si1, si2 );
					if (sj1 < j1 && sj2 < j2) {
						interaction.basePairs.push_back( energy.getBasePair(sj1,sj2) );
					}
					// sort output interaction
					interaction.sort();

					// stop searching for seeds
					return;
				}
			}

		} // si2
	} // si1

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dSeedExtension::
getNextBest( Interaction & curBest )
{
	INTARNA_NOT_IMPLEMENTED("PredictorMfeEns2dSeedExtension::getNextBest() not implemented yet");
}

//////////////////////////////////////////////////////////////////////////


} // namespace
