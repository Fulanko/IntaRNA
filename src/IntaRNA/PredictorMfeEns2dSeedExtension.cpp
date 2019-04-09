
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

	// TODO check if getMaxLength > sqrt3(size_t) -> error
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
		const size_t maxMatrixLen1 = energy.getAccessibility1().getMaxLength()-sl1+1;
		const size_t maxMatrixLen2 = energy.getAccessibility2().getMaxLength()-sl2+1;
		// check if seed fits into interaction range
		if (sj1 > interaction_size1 || sj2 > interaction_size2)
			continue;

		// EL
		hybridZ_left.resize( std::min(si1+1, maxMatrixLen1), std::min(si2+1, maxMatrixLen2) );
		fillHybridZ_left(si1, si2, outConstraint);

		// ER
		hybridZ_right.resize( std::min(interaction_size1-sj1, maxMatrixLen1), std::min(interaction_size2-sj2, maxMatrixLen2) );
		fillHybridZ_right(sj1, sj2, outConstraint);

		// update Optimum for all boundary combinations
		for (size_t i1 = 0; i1<hybridZ_left.size1(); i1++) {
			// ensure max interaction length in seq 1
			for (size_t j1 = 0; j1 < hybridZ_right.size1() ; j1++) {
				// check interaction length
				if (sj1+j1-si1+i1 >= energy.getAccessibility1().getMaxLength()) continue;
				for (size_t i2 = 0; i2< hybridZ_left.size2(); i2++) {
					// check complementarity of boundary
					if ( Z_equal(hybridZ_left(i1,i2), 0.0) ) continue;
					// ensure max interaction length in seq 2
					for (size_t j2 = 0; j2 < hybridZ_right.size2() ; j2++) {
						// check interaction length
						if (sj2+j2-si2+i2 >= energy.getAccessibility2().getMaxLength()) continue;
						// check complementarity of boundary
						if (Z_equal(hybridZ_right(j1,j2),0.0)) continue;
						// compute partition function given the current seed
						updateZ(si1-i1, sj1+j1, si2-i2, sj2+j2, hybridZ_left(i1,i2) * seedZ * hybridZ_right(j1,j2), true);
					} // dj2
				} // di2
			} // dj1
		} // di1

	} // si1 / si2

	// update ensemble mfe
	for (std::unordered_map<size_t, ZPartition >::const_iterator it = Z_partitions.begin(); it != Z_partitions.end(); ++it)
	{
		PredictorMfe::updateOptima( it->second.i1, it->second.j1, it->second.i2, it->second.j2, energy.getE(it->second.partZ), true );
	}

	LOG(DEBUG) << "Overall Z: " << getHybridZ();

	// report mfe interaction
	reportOptima( outConstraint );

}

//////////////////////////////////////////////////////////////////////////

E_type
PredictorMfeEns2dSeedExtension::
getNonOverlappingEnergy( const size_t si1, const size_t si2, const size_t si1p, const size_t si2p ) {

#if INTARNA_IN_DEBUG_MODE
	// check indices
	if( !seedHandler.isSeedBound(si1,si2) )
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::getNonOverlappingEnergy( si "+toString(si1)+","+toString(si2)+",..) is no seed bound");
	if( !seedHandler.isSeedBound(si1p,si2p) )
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::getNonOverlappingEnergy( sip "+toString(si1p)+","+toString(si2p)+",..) is no seed bound");
#endif

	// sanity check
	if( ! seedHandler.areLoopOverlapping(si1,si2,si1p,si2p) || si1 == si1p ) {
		return E_type(0);
	}

	// trace S
	Interaction interaction = Interaction(energy.getAccessibility1().getSequence(), energy.getAccessibility2().getAccessibilityOrigin().getSequence());
	interaction.basePairs.push_back( energy.getBasePair(si1, si2) );
	seedHandler.traceBackSeed( interaction, si1, si2 );

	E_type fullE = 0;
	size_t k1old = si1, k2old = si2;
	for (size_t i = 1; i < interaction.basePairs.size(); i++) {
		// get index of current base pair
		size_t k1 = energy.getIndex1(interaction.basePairs[i]);
		// check if overlap done
		if (k1 > si1p) break;
		size_t k2 = energy.getIndex2(interaction.basePairs[i]);
		// add hybridization energy
		fullE += energy.getE_interLeft(k1old,k1,k2old,k2);
		// store
		k1old = k1;
		k2old = k2;
	}
	return fullE;
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

				E_type seedE = seedHandler.getSeedE(i1, i2);
				if (E_isNotINF(seedE)) {
					// iterate seeds in S region
					// todo maybe rename in si1prime etc
					size_t sj1 = RnaSequence::lastPos, sj2 = RnaSequence::lastPos;
					// todo replace with left-most loop-overlapping seed
					size_t si1overlap = i1, si2overlap = i2;
					std::map<size_t, std::tuple<size_t, size_t>> overlappingSeeds;
					// find left-most loop-overlapping seed
					while( seedHandler.updateToNextSeed(sj1,sj2
							, i1, i1+seedHandler.getSeedLength1(i1, i2)-2
							, i2, i2+seedHandler.getSeedLength2(i1, i2)-2) )
					{
						// check for overlapping seeds
						if (seedHandler.areLoopOverlapping(i1, i2, sj1, sj2) && i1 < sj1) {
							// update left-most loop-overlapping seed
							overlappingSeeds.insert(std::pair<size_t, std::tuple<size_t, size_t>>(std::min(sj1, sj2), std::make_tuple(sj1, sj2)));
						}
					}

					//if (i1 < si1overlap) {
					if (overlappingSeeds.size() > 0) {
					  si1overlap = std::get<0>(overlappingSeeds.begin()->second);
						si2overlap = std::get<1>(overlappingSeeds.begin()->second);
						// compute Energy of loop S \ S'
						E_type nonOverlapE = getNonOverlappingEnergy(i1, i2, si1overlap, si2overlap);
						LOG(DEBUG) << "overlap | Energy = " << nonOverlapE;
						// substract energy.getBoltzmannWeight( nonOverlapE ) * hybridZ_left( si1overlap, si2overlap bis anchor seed [==1 falls gleich])
//						subtractList.push_back(std::make_tuple(j1-i1, j2-i2, nonOverlapE));
						// TODO: case where overlap == anchor
						hybridZ_left(j1-i1, j2-i2) -= energy.getBoltzmannWeight( nonOverlapE ) * hybridZ_left( j1-si1overlap, j2-si2overlap);
					} else {
						// if no S'
						LOG(DEBUG) << "non overlap | Energy = " << energy.getBoltzmannWeight(seedE) * hybridZ_left( j1-(i1+seedHandler.getSeedLength1(i1, i2)-1), j2-(i2+seedHandler.getSeedLength2(i1, i2)-1));
						// substract seedZ * hybridZ_left(right end seed bis anchor seed)
//						subtractList.push_back(std::make_tuple(j1-i1, j2-i2, hybridZ_left(j1-i1, j2-i2)));
            hybridZ_left(j1-i1, j2-i2) -= energy.getBoltzmannWeight(seedE) * hybridZ_left( j1-(i1+seedHandler.getSeedLength1(i1, i2)-1), j2-(i2+seedHandler.getSeedLength2(i1, i2)-1));
					}

				}

			} // complementary

		} // i2
	} // i1

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

	size_t si1 = RnaSequence::lastPos, si2 = RnaSequence::lastPos;
	while( seedHandler.updateToNextSeed(si1,si2
			, i1,j1+1-seedHandler.getConstraint().getBasePairs()
			, i2,j2+1-seedHandler.getConstraint().getBasePairs() ) )
	{
		E_type seedE = seedHandler.getSeedE(si1, si2);

		const size_t sl1 = seedHandler.getSeedLength1(si1, si2);
		const size_t sl2 = seedHandler.getSeedLength2(si1, si2);
		const size_t sj1 = si1+sl1-1;
		const size_t sj2 = si2+sl2-1;
		const size_t maxMatrixLen1 = energy.getAccessibility1().getMaxLength()-sl1+1;
		const size_t maxMatrixLen2 = energy.getAccessibility2().getMaxLength()-sl2+1;

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
