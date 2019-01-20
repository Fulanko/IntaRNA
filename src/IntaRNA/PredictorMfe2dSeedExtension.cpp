
#include "IntaRNA/PredictorMfe2dSeedExtension.h"

namespace IntaRNA {

//////////////////////////////////////////////////////////////////////////

PredictorMfe2dSeedExtension::
PredictorMfe2dSeedExtension(
		const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker
		, SeedHandler * seedHandlerInstance )
 :
	PredictorMfe2d(energy,output,predTracker)
	, seedHandler(seedHandlerInstance)
	, hybridE_right( 0,0 )
{
	assert( seedHandler.getConstraint().getBasePairs() > 1 );
}

//////////////////////////////////////////////////////////////////////////

PredictorMfe2dSeedExtension::
~PredictorMfe2dSeedExtension()
{
}

//////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dSeedExtension::
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
		throw std::runtime_error("PredictorMfe2dSeedExtension : the enumeration of non-overlapping suboptimal interactions is not supported in this prediction mode");
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
	LOG(DEBUG) << "interaction size: " << interaction_size1 << " / " << interaction_size2 ;

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
		LOG(DEBUG) << "FOUND SEED: " << seedHandler.getSeedE(si1, si2) << " at: " <<si1 <<"~" << energy.getBasePair(si1,si2).first <<", " <<si2 <<"~" << energy.getBasePair(si1,si2).second ;
		const size_t sl1 = seedHandler.getSeedLength1(si1, si2)-1;
		const size_t sl2 = seedHandler.getSeedLength2(si1, si2)-1;
		const size_t sj1 = si1+sl1;
		const size_t sj2 = si2+sl2;
		LOG(DEBUG) << "SEED: " << sj1 <<","<<sj2<<" length "<<sl1<<","<<sl2;
		// check if seed fits into interaction range
		if (sj1 > interaction_size1 || sj2 > interaction_size2)
			continue;

		// EL
		LOG(DEBUG) <<"compute EL";
		hybridE_pq.resize( si1+1, si2+1 );
		fillHybridE(si1, si2, outConstraint, 0, 0);

		// ER
		LOG(DEBUG) <<"compute ER";
		hybridE_right.resize( interaction_size1-sj1, interaction_size2-sj2);
		fillHybridE_right(sj1, sj2, outConstraint, interaction_size1-1, interaction_size2-1);

		// update Optimum for all boundary combinations
		for (int i1 = 0; i1<=si1; i1++) {
			// ensure max interaction length in seq 1
			for (int j1 = 0; j1 < hybridE_right.size1() ; j1++) {
				if(si1-i1+sl1+j1 > energy.getAccessibility1().getMaxLength()) continue;
				for (int i2 = 0; i2<=si2; i2++) {
					if (E_isINF(hybridE_pq(i1,i2))) continue;
					// ensure max interaction length in seq 2
					for (int j2 = 0; j2 < hybridE_right.size2() ; j2++) {
						if(si2-i2+sl2+j2 > energy.getAccessibility2().getMaxLength()) continue;
						if (E_isINF(hybridE_right(j1,j2))) continue;

						E_type fullE = seedE;
						fullE += hybridE_pq(i1,i2); // contains E_init
						fullE += hybridE_right(j1,j2);
						// LOG(DEBUG) << i1 << ":" << sj1+j1 << ":" << i2 << ":" << sj2+j2 << "__" << hybridE_pq(i1,i2) << "___" <<hybridE_right(j1,j2) << " --> " << energy.getE(i1, sj1+j1, i2, sj2+j2, fullE);
						PredictorMfe::updateOptima( i1, sj1+j1, i2, sj2+j2, fullE, true );
					} // j2
				} // i2
			} // j1
		} // i1

	} // si2
	} // si1

	// report mfe interaction
	reportOptima( outConstraint );

}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dSeedExtension::
fillHybridE_right( const size_t i1, const size_t i2
			, const OutputConstraint & outConstraint
			, const size_t j1max, const size_t j2max )
{

	// global vars to avoid reallocation
	size_t j1,j2,k1,k2;
	//////////  FIRST ROUND : COMPUTE HYBRIDIZATION ENERGIES ONLY  ////////////

	// current minimal value
	E_type curMinE = E_INF;
	// iterate over all window starts j1 (seq1) and j2 (seq2)
	for (j1=i1; j1 <= j1max; j1++ ) {
		// w1 width check obsolete due to hybridErangeRight setup
		// screen for left boundaries in seq2
		for (j2=i2; j2 <= j2max; j2++ ) {

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
				// updateOptima( i1,j1,i2,j2, hybridE_pq_right(j1,j2), true );
				continue;
			}
		}
	}

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dSeedExtension::
traceBack( Interaction & interaction, const OutputConstraint & outConstraint  )
{
	// check if something to trace
	if (interaction.basePairs.size() < 2) {
		return;
	}

#if INTARNA_IN_DEBUG_MODE
	// sanity checks
	if ( interaction.basePairs.size() != 2 ) {
		throw std::runtime_error("PredictorMfe2dSeedExtension::traceBack() : given interaction does not contain boundaries only");
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
		throw std::runtime_error("PredictorMfe2dSeedExtension::traceBack() : given interaction not valid");
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
		throw std::runtime_error("PredictorMfe2dSeedExtension::traceBack() : given boundaries "+toString(interaction)+" can not hold a seed of "+toString(seedHandler.getConstraint().getBasePairs())+" base pairs");
	}
#endif

	LOG(DEBUG) << "traceback" ;

	LOG(DEBUG) << i1 << "/" << i2 << "/" << j1 << "/" << j2;

	E_type fullE = interaction.energy;

	LOG(DEBUG) << "interactionEnergy" << fullE;

	for (size_t si1 = i1; si1 <= j1+1-seedHandler.getConstraint().getBasePairs(); si1++) {
		for (size_t si2 = i2; si2 <= j2+1-seedHandler.getConstraint().getBasePairs(); si2++) {
			E_type seedE = seedHandler.getSeedE(si1, si2);
			if ( E_isNotINF( seedE ) ) {

				const size_t sl1 = seedHandler.getSeedLength1(si1, si2)-1;
				const size_t sl2 = seedHandler.getSeedLength2(si1, si2)-1;
				const size_t sj1 = si1+sl1;
				const size_t sj2 = si2+sl2;

				LOG(DEBUG) << "seed at " << si1 << "  " << si2;
				LOG(DEBUG) << "start at " << i1 << "  " << i2;

				hybridE_pq.resize( si1+1, si2+1 );
				fillHybridE( si1, si2, outConstraint, 0, 0 );
				hybridE_right.resize( j1-sj1+1, j2-sj2+1 );
				fillHybridE_right( sj1, sj2, outConstraint, j1, j2 );

				LOG(DEBUG) << seedE << "::" << hybridE_pq( i1, i2 ) << "::" << hybridE_right( j1-sj1, j2-sj2 );
				LOG(DEBUG) << energy.getE(i1, j1, i2, j2, seedE + hybridE_pq( i1, i2 ) + hybridE_right( j1-sj1, j2-sj2 ));

				if ( E_equal( fullE,
						(energy.getE(i1, j1, i2, j2, seedE + hybridE_pq( i1, i2 ) + hybridE_right( j1-sj1, j2-sj2 )))))
				{
					// found seed -> traceback
					LOG(DEBUG) << "found seed";

					// the currently traced value for i1-si1, i2-si2
					E_type curE = hybridE_pq(i1,i2);

					// trace back left
					while( i1 != si1 ) {

						// check if just internal loop
						if ( E_equal( curE, (energy.getE_interLeft(i1,si1,i2,si2) + hybridE_pq(si1,si2)) ) )
						{
							break;
						}
						// check all interval splits
						if ( (si1-i1) > 1 && (si2-i2) > 1) {
							// temp variables
							size_t k1,k2;
							bool traceNotFound = true;
							// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
							for (k1=std::min(si1-1,i1+energy.getMaxInternalLoopSize1()+1); traceNotFound && k1>i1; k1--) {
							for (k2=std::min(si2-1,i2+energy.getMaxInternalLoopSize2()+1); traceNotFound && k2>i2; k2--) {
								// check if (k1,k2) are valid left boundary
								if ( E_isNotINF( hybridE_pq(k1,k2) ) ) {
									// LOG(DEBUG) << (energy.getE_interLeft(i1,k1,i2,k2) + hybridE_pq(k1,k2));
									if ( E_equal( curE,
											(energy.getE_interLeft(i1,k1,i2,k2) + hybridE_pq(k1,k2)) ) )
									{
										// stop searching
										traceNotFound = false;
										// store splitting base pair
										LOG(DEBUG) << "update left at: " << k1 << ":" << k2;
										interaction.basePairs.push_back( energy.getBasePair(k1,k2) );
										// trace right part of split
										i1=k1;
										i2=k2;
										curE = hybridE_pq(i1,i2);
									}
								}
							}
							}
						}

				  } // traceback left

					// trace seed
					LOG(DEBUG) << "seed at:" << si1 << ":" << si2;
					if (si1 > i1 && si2 > i2) {
						interaction.basePairs.push_back( energy.getBasePair(si1,si2) );
					}
					seedHandler.traceBackSeed( interaction, si1, si2 );
					if (sj1 < j1 && sj2 < j2) {
						interaction.basePairs.push_back( energy.getBasePair(sj1,sj2) );
					}

					// the currently traced value for sj1-j1, sj2-j2
					curE = hybridE_right(j1-sj1,j2-sj2);

					// trace back right
					while( j1 != sj1 ) {

						// check if just internal loop
						if ( E_equal( curE, (energy.getE_interLeft(sj1,j1,sj2,j2) + hybridE_right(0,0)) ) )
						{
							break;
						}
						// check all interval splits
						if ( (j1-sj1) > 1 && (j2-sj2) > 1) {
							// temp variables
							size_t k1,k2;
							bool traceNotFound = true;
							// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
							for (k1=std::min(j1-1,sj1+energy.getMaxInternalLoopSize1()+1); traceNotFound && k1>sj1; k1--) {
							for (k2=std::min(j2-1,sj2+energy.getMaxInternalLoopSize2()+1); traceNotFound && k2>sj2; k2--) {
								// check if (k1,k2) are valid left boundary
								if ( E_isNotINF( hybridE_right(k1-sj1,k2-sj2) ) ) {
									if ( E_equal( curE,
											(energy.getE_interLeft(k1,j1,k2,j2) + hybridE_right(k1-sj1,k2-sj2)) ) )
									{
										// stop searching
										traceNotFound = false;
										// store splitting base pair
										LOG(DEBUG) << "update right at: " << k1 << ":" << k2;
										interaction.basePairs.push_back( energy.getBasePair(k1,k2) );
										// trace right part of split
										j1=k1;
										j2=k2;
										curE = hybridE_right(j1-sj1,j2-sj2);
									}
								}
							}
							}
						}
					}  // traceback right

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
PredictorMfe2dSeedExtension::
getNextBest( Interaction & curBest )
{
	INTARNA_NOT_IMPLEMENTED("PredictorMfe2dSeedExtension::getNextBest() not implemented yet");
}

//////////////////////////////////////////////////////////////////////////


} // namespace
