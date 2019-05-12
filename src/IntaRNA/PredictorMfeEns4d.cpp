
#include "IntaRNA/PredictorMfeEns4d.h"

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////////

PredictorMfeEns4d::
PredictorMfeEns4d(
		const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker )
 : PredictorMfeEns(energy,output,predTracker)
	, hybridZ( 0,0 )
	, Z(0.0)
	, maxProbInteraction(energy.getAccessibility1().getSequence()
			,energy.getAccessibility2().getAccessibilityOrigin().getSequence())
{
}


////////////////////////////////////////////////////////////////////////////

PredictorMfeEns4d::
~PredictorMfeEns4d()
{
	// clean up
	this->clear();
}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns4d::
predict( const IndexRange & r1
		, const IndexRange & r2
		, const OutputConstraint & outConstraint
		)
{
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"predicting maximally probable interactions in O(n^4) space..."; }
	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

	if (outConstraint.reportMax > 1) {
		INTARNA_NOT_IMPLEMENTED("PredictorMfeEns4d::predict(reportMax > 1) : not implemented");
	}

#if INTARNA_IN_DEBUG_MODE
	// check indices (both regions ascending due to reversing of seq2)
	if (!(r1.isAscending() && r2.isAscending()) )
		throw std::runtime_error("PredictorMfeEns4d::predict("+toString(r1)+","+toString(r2)+") is not sane");
#endif

	// clear data
	clear();

	// setup index offset
	energy.setOffset1( r1.from );
	energy.setOffset2( r1.from );

	// resize matrix
	hybridZ.resize( std::min( energy.size1()
						, (r1.to==RnaSequence::lastPos?energy.size1()-1:r1.to)-r1.from+1 )
				, std::min( energy.size2()
						, (r2.to==RnaSequence::lastPos?energy.size2()-1:r2.to)-r2.from+1 ) );

	size_t debug_count_cells_null=0, debug_count_cells_nonNull = 0, debug_cellNumber=0;

	size_t maxWidthFori1i2 = 0;

	bool i1blocked, i1or2blocked;
	// initialize 3rd and 4th dimension of the matrix
	for (size_t i1=0; i1<hybridZ.size1(); i1++) {
		// check if i1 is blocked for interaction
		i1blocked = !energy.isAccessible1(i1);
	for (size_t i2=0; i2<hybridZ.size2(); i2++) {
		// check whether i1 or i2 is blocked for interaction
		i1or2blocked = i1blocked || !energy.isAccessible2(i2);

		if (hybridZ.size1()-i1 < hybridZ.size2()-i2) {
			maxWidthFori1i2 = getMaxInteractionWidth( hybridZ.size1()-i1, energy.getMaxInternalLoopSize1() );
		} else {
			maxWidthFori1i2 = getMaxInteractionWidth( hybridZ.size2()-i2, energy.getMaxInternalLoopSize2() );
		}

		debug_cellNumber =
				/*w1 = */ std::min(energy.getAccessibility1().getMaxLength(), std::min( hybridZ.size1()-i1, maxWidthFori1i2) )
			*	/*w2 = */ std::min(energy.getAccessibility2().getMaxLength(), std::min( hybridZ.size2()-i2, maxWidthFori1i2) );

		// check if i1 and i2 are not blocked and can form a base pair
		if ( ! i1or2blocked
			&& energy.areComplementary( i1, i2 ))
		{
			// create new 2d matrix for different interaction site widths
			hybridZ(i1,i2) = new Z2dMatrix(
				/*w1 = */ std::min(energy.getAccessibility1().getMaxLength(), std::min( hybridZ.size1()-i1, maxWidthFori1i2) ),
				/*w2 = */ std::min(energy.getAccessibility2().getMaxLength(), std::min( hybridZ.size2()-i2, maxWidthFori1i2) ));

			debug_count_cells_nonNull += debug_cellNumber;
		} else {
			// reduce memory consumption and avoid computation for this start index combination
			hybridZ(i1,i2) = NULL;
			debug_count_cells_null += debug_cellNumber;
		}
	}
	}
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ LOG(DEBUG) <<"init 4d matrix : "<<debug_count_cells_nonNull <<" to be filled ("
				<<((double)debug_count_cells_nonNull/(double)(debug_count_cells_nonNull+debug_count_cells_null))
				<<"%) and "<<debug_count_cells_null <<" not allocated"; }

	// initialize max prob interaction for updates
	initOptima( outConstraint );
	// initialize overall partition function for updates
	initZ( outConstraint );

	// fill matrix
	fillHybridZ( );

	// update ensemble mfe
	for (std::unordered_map<size_t, ZPartition >::const_iterator it = Z_partitions.begin(); it != Z_partitions.end(); ++it)
	{
		// if partition function is > 0
		if (Z_isNotINF(it->second.partZ) && it->second.partZ > 0) {
			//LOG(DEBUG) << "partZ: " << it->second.partZ;
			updateOptima( it->second.i1, it->second.j1, it->second.i2, it->second.j2, it->second.partZ, true );
		}
	}

	LOG(DEBUG) <<"Overall Z = "<< getOverallZ();
	LOG(DEBUG) <<"Overall E = "<<E_2_Ekcal(energy.getE(getOverallZ()));

	// report interaction site with maximal probability
	reportOptima( outConstraint );
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns4d::
clear()
{
	// delete 3rd and 4th dimension of the matrix
	for (Z4dMatrix::iterator1 iRows = hybridZ.begin1(); iRows != hybridZ.end1(); iRows++) {
		for (Z4dMatrix::iterator2 ijEntry = iRows.begin(); ijEntry != iRows.end(); ijEntry++) {
			// delete 2d matrix for current ij
			 INTARNA_CLEANUP(*ijEntry);
		}
	}
	// clear matrix, free data
	hybridZ.clear();

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns4d::
fillHybridZ()
{

	// global vars to avoid reallocation
	size_t i1,i2,j1,j2,w1,w2,k1,k2;

	//////////  FIRST ROUND : COMPUTE HYBRIDIZATION ENERGIES ONLY  ////////////

	// reset overall partition function
	Z = 0.0;

	// current Z value
	Z_type curZ = 0.0;
	// iterate increasingly over all window sizes w1 (seq1) and w2 (seq2)
	for (w1=0; w1<energy.getAccessibility1().getMaxLength(); w1++) {
	for (w2=0; w2<energy.getAccessibility2().getMaxLength(); w2++) {
		// iterate over all window starts i1 (seq1) and i2 (seq2)
		// TODO PARALLELIZE THIS DOUBLE LOOP ?!
		for (i1=0; i1+w1<hybridZ.size1(); i1++) {
		for (i2=0; i2+w2<hybridZ.size2(); i2++) {
			// check if left boundary is complementary
			// and widths are possible
			if (hybridZ(i1,i2) == NULL || hybridZ(i1,i2)->size1()<=w1 || hybridZ(i1,i2)->size2()<=w2) {
				// interaction not possible: nothing to do
				continue;
			}
			// check if interaction exceeds possible width due to max-loop-length
			if ( getMaxInteractionWidth( 1+w1, energy.getMaxInternalLoopSize1() ) < w2
				|| getMaxInteractionWidth( 1+w2, energy.getMaxInternalLoopSize2() ) < w1)
			{
				// ignore this entry
				(*hybridZ(i1,i2))(w1,w2) = 0;
				continue;
			}

			// get window ends j (seq1) and l (seq2)
			j1=i1+w1;
			j2=i2+w2;

			// check if right boundary is complementary
			if (hybridZ(j1,j2) == NULL) {
				// not complementary -> ignore this entry
				(*hybridZ(i1,i2))(w1,w2) = 0;
				continue;
			}
			// compute entry
			curZ = 0;

			// either interaction initiation
			if ( w1==0 && w2==0 )  {
				curZ += energy.getBoltzmannWeight(energy.getE_init());
			} else {
			// or only internal loop energy (nothing between i and j)
				if ( (w1+1) <= energy.getMaxInternalLoopSize1() && (w2+1) <= energy.getMaxInternalLoopSize2()) {
					curZ += energy.getBoltzmannWeight(energy.getE_interLeft(i1,j1,i2,j2))
						* (*hybridZ(j1,j2))(0,0);
				}
			}

			if (w1 > 1 && w2 > 1) {
				// sum all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				for (k1=std::min(j1-1,i1+energy.getMaxInternalLoopSize1()+1); k1>i1; k1--) {
				for (k2=std::min(j2-1,i2+energy.getMaxInternalLoopSize2()+1); k2>i2; k2--) {
					// check if (k1,k2) are complementary
					if (hybridZ(k1,k2) != NULL && hybridZ(k1,k2)->size1()>(j1-k1) && hybridZ(k1,k2)->size2()>(j2-k2)) {
						curZ += energy.getBoltzmannWeight(energy.getE_interLeft(i1,k1,i2,k2))
								* ((*hybridZ(k1,k2))(j1-k1,j2-k2));
					}
				}
				}
			}
			// store value
			(*hybridZ(i1,i2))(w1,w2) = curZ;
			// update max prob interaction
			updateZ( i1,j1,i2,j2, (*hybridZ(i1,i2))(w1,w2), true );
		}
		}
	}
	}

}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns4d::
initOptima( const OutputConstraint & outConstraint )
{
	// initialize max prob interaction (partition function value)
	maxProbInteraction.energy = 0.0;
	// reset boundary base pairs
	maxProbInteraction.r1.from = RnaSequence::lastPos;
	maxProbInteraction.r1.to = RnaSequence::lastPos;
	maxProbInteraction.r2.from = RnaSequence::lastPos;
	maxProbInteraction.r2.to = RnaSequence::lastPos;

	if (outConstraint.reportMax > 1) {
		INTARNA_NOT_IMPLEMENTED("PredictorMfeEns4d::initOptima(reportMax > 1)");
	}
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns4d::
updateOptima( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const E_type E
		, const bool isHybridE )
{
	INTARNA_NOT_IMPLEMENTED("this function should not be called ..");
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns4d::
updateOptima( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const Z_type interZ
		, const bool isHybridZ )
{
//		{ LOG(DEBUG) <<"Z( "<<i1<<"-"<<j1<<", "<<i2<<"-"<<j2<<" ) = "
//						<<curZ
//						<<" = " <<(eH + eE + eD); }

	// add Boltzmann weights of all penalties
	Z_type curZ = isHybridZ ? interZ * energy.getBoltzmannWeight( energy.getE(i1,j1,i2,j2,0.0) ) : interZ;
	E_type curE = energy.getE(curZ);

	// report call if needed
	if (predTracker != NULL) {
		// inform about prediction
		predTracker->updateOptimumCalled( i1 + energy.getOffset1()
										, j1 + energy.getOffset1()
										, i2 + energy.getOffset2()
										, j2 + energy.getOffset2()
										, curE );
	}

	// update overall partition function
	Z += (double)curZ;

	if (curE < maxProbInteraction.energy) {
		// store new global min
		maxProbInteraction.energy = curE;
		// store interaction boundaries
		// left
		maxProbInteraction.r1.from = i1+energy.getOffset1();
		maxProbInteraction.r2.from = energy.getAccessibility2().getReversedIndex(i2+energy.getOffset2());
		// right
		maxProbInteraction.r1.to = j1+energy.getOffset1();
		maxProbInteraction.r2.to = energy.getAccessibility2().getReversedIndex(j2+energy.getOffset2());
	}
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns4d::
traceBack( Interaction & interaction, const OutputConstraint & outConstraint  )
{
	// check if something to trace
	if (interaction.basePairs.size() < 2) {
		return;
	}

#if INTARNA_IN_DEBUG_MODE
	// sanity checks
	if ( interaction.basePairs.size() != 2 ) {
		throw std::runtime_error("PredictorMfeEns2d::traceBack() : given interaction does not contain boundaries only");
	}
#endif

	// ensure sorting
	interaction.sort();

	// check for single base pair interaction
	if (interaction.basePairs.at(0).first == interaction.basePairs.at(1).first) {
		// delete second boundary (identical to first)
		interaction.basePairs.resize(1);
		// update done
		return;
	}

#if INTARNA_IN_DEBUG_MODE
	// sanity checks
	if ( ! interaction.isValid() ) {
		throw std::runtime_error("PredictorMfeEns2d::traceBack() : given interaction is not valid");
	}
#endif

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns4d::
getNextBest( Interaction & curBest )
{
	curBest.energy = E_INF;
	curBest.basePairs.clear();
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns4d::
reportOptima( const OutputConstraint & outConstraint )
{

	if (outConstraint.reportMax == 0) {
		return;
	}

	if (outConstraint.reportMax > 1) {
		INTARNA_NOT_IMPLEMENTED("PredictorMfeEns4d::reportOptima(reportMax > 1)");
	}

	// push to output handler
	output.add( maxProbInteraction, outConstraint );

}

////////////////////////////////////////////////////////////////////////////


} // namespace
