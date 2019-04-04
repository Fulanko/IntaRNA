
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
	PredictorMfeEns(energy,output,predTracker)
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

	size_t si1 = RnaSequence::lastPos, si2 = RnaSequence::lastPos;
	while( seedHandler.updateToNextSeed(si1,si2
			, 0, interaction_size1+1-seedHandler.getConstraint().getBasePairs()
			, 0, interaction_size2+1-seedHandler.getConstraint().getBasePairs()) )
	{
		E_type seedE = seedHandler.getSeedE(si1, si2);
		const Z_type seedZ = energy.getBoltzmannWeight( seedE );

		const size_t sl1 = seedHandler.getSeedLength1(si1, si2);
		const size_t sl2 = seedHandler.getSeedLength2(si1, si2);
		const size_t sj1 = si1+sl1-1;
		const size_t sj2 = si2+sl2-1;
		const size_t maxMatrixLen1 = energy.getAccessibility1().getMaxLength()-sl1;
		const size_t maxMatrixLen2 = energy.getAccessibility2().getMaxLength()-sl2;
		// check if seed fits into interaction range
		if (sj1 > interaction_size1 || sj2 > interaction_size2)
			continue;

		// EL
		hybridZ_left.resize( std::min(si1+1, maxMatrixLen1), std::min(si2+1, maxMatrixLen2) );
		fillHybridZ_left(si1, si2, outConstraint);

		//printMatrix(hybridZ_left);

		// ER
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
						// compute overall ensemble energy
						E_type fullE = seedE + energy.getE(hybridZ_left(i1,i2)) + energy.getE(hybridZ_right(j1,j2));
						// update ensemble mfe
						PredictorMfe::updateOptima( si1-i1, sj1+j1, si2-i2, sj2+j2, fullE, true );
						// update Z
						updateZ(si1-i1, sj1+j1, si2-i2, sj2+j2, fullE, true);
						// store partial Z
						size_t key = getHashKey(si1-i1, sj1+j1, si2-i2, sj2+j2);
						if ( Z_partitions.find(key) == Z_partitions.end() ) {
						  Z_partitions[key] = fullE;
						} else {
						  Z_partitions[key] += fullE;
						}
					} // dj2
				} // di2
			} // dj1
		} // di1

	} // si1 / si2

	std::cout << "Z: " << getHybridZ() << std::endl;
	for (std::unordered_map<size_t, Z_type >::const_iterator it = Z_partitions.begin(); it != Z_partitions.end(); ++it)
  {
    std::cout << it->first << " " << it->second << "\n";
  }

	// report mfe interaction
	reportOptima( outConstraint );

}

////////////////////////////////////////////////////////////////////////////

size_t
PredictorMfeEns2dSeedExtension::
getHashKey( const size_t i1, const size_t j1, const size_t i2, const size_t j2) {
	size_t maxLength = std::max(energy.getAccessibility1().getMaxLength(), energy.getAccessibility2().getMaxLength());
	size_t key = 0;
	key += i1 + pow(maxLength, 0);
	key += j1 + pow(maxLength, 1);
	key += i2 + pow(maxLength, 2);
	key += j2 + pow(maxLength, 3);
	return key;
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dSeedExtension::
fillHybridZ_left( const size_t j1, const size_t j2
			, const OutputConstraint & outConstraint )
{
#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!energy.areComplementary(j1,j2) )
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::fillHybridZ_left("+toString(j1)+","+toString(j2)+",..) are not complementary");
#endif

	// global vars to avoid reallocation
	size_t i1,i2,k1,k2;
	std::cout << "fill left" << std::endl;

	// position of last seen complementary seed
	size_t last_si1 = j1;
	size_t last_si2 = j2;

	// iterate over all window starts i1 (seq1) and i2 (seq2)
	for (i1=j1; j1-i1 < hybridZ_left.size1(); i1-- ) {
		for (i2=j2; j2-i2 < hybridZ_left.size2(); i2-- ) {
			// init current cell (0 if not just right-most (j1,j2) base pair)
			hybridZ_left(j1-i1,j2-i2) = i1==j1 && i2==j2 ? energy.getBoltzmannWeight(energy.getE_init()) : 0.0;

			// check if complementary (use global sequence indexing)
			if( i1<j1 && i2<j2 && energy.areComplementary(i1,i2) ) {

				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				for (k1=i1; k1++ < j1; ) {
					// ensure maximal loop length
					if (k1-i1 > energy.getMaxInternalLoopSize1()+1) break;
					for (k2=i2; k2++ < j2; ) {
						// ensure maximal loop length
						if (k2-i2 > energy.getMaxInternalLoopSize2()+1) break;
						// check if (k1,k2) are valid left boundary
						if ( ! Z_equal(hybridZ_left(j1-k1,j2-k2), 0.0) ) {
							hybridZ_left(j1-i1,j2-i2) += energy.getBoltzmannWeight(energy.getE_interLeft(i1,k1,i2,k2)) * hybridZ_left(j1-k1,j2-k2);
						}
					} // k2
				} // k1

				// subtract seeds in left matrix
				E_type seedE = seedHandler.getSeedE(i1, i2);
				if (E_isNotINF(seedE)) {
					std::cout << "Complementary SEED at " << i1 << ":" << i2 << "= " << seedE << std::endl;
					const size_t sl1 = seedHandler.getSeedLength1(i1, i2)-1;
					const size_t sl2 = seedHandler.getSeedLength2(i1, i2)-1;
					const size_t sj1 = std::min(i1+sl1, last_si1);
					const size_t sj2 = std::min(i2+sl2, last_si2);
					std::cout << hybridZ_left(j1-i1,j2-i2) << std::endl;
					std::cout << "sub: " << energy.getBoltzmannWeight(energy.getE_interLeft(i1,sj1,i2,sj2)) * hybridZ_left(j1-sj1,j2-sj2) << std::endl;
					hybridZ_left(j1-i1,j2-i2) -= energy.getBoltzmannWeight(energy.getE_interLeft(i1,sj1,i2,sj2)) * hybridZ_left(j1-sj1,j2-sj2);
					std::cout << hybridZ_left(j1-i1,j2-i2) << std::endl;
					last_si1 = i1;
					last_si2 = i2;
				}
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
	std::cout << "===========================================================\n";

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

	size_t si1 = RnaSequence::lastPos, si2 = RnaSequence::lastPos;
	while( seedHandler.updateToNextSeed(si1,si2
			, i1,j1+1-seedHandler.getConstraint().getBasePairs()
			, i2,j2+1-seedHandler.getConstraint().getBasePairs() ) )
	{
		E_type seedE = seedHandler.getSeedE(si1, si2);

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
			seedHandler.addSeeds( interaction );

			// stop searching for seeds
			return;
		}

	} // si1 / si2

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
