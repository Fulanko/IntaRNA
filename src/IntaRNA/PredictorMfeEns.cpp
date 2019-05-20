
#include "IntaRNA/PredictorMfeEns.h"

#include <iostream>
#include <algorithm>

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////////

PredictorMfeEns::PredictorMfeEns(
		const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker
		)
	: PredictorMfe(energy,output,predTracker)
	, overallZ(0)
{

}

////////////////////////////////////////////////////////////////////////////

PredictorMfeEns::~PredictorMfeEns()
{
}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns::
initZ( const OutputConstraint & outConstraint )
{
	// reinit overall partition function
	overallZ = Z_type(0);

	// reinit mfe information
	initOptima( outConstraint );
}

////////////////////////////////////////////////////////////////////////////

Z_type
PredictorMfeEns::
getOverallZ() const
{
	return overallZ;
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns::
checkKeyBoundaries( const size_t maxLength )
{
	// check if getMaxLength > sqrt3(size_t) -> error
	if (maxLength > cbrt(std::numeric_limits<size_t>::max())) {
		throw std::runtime_error("PredictorMfeEns maxLength too big for key generation (out of bounds)");
	}
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns::
updateZ( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const Z_type partZ
		, const bool isHybridZ )
{
	// check if something to be done
	if (Z_equal(partZ,0) || Z_isINF(overallZ))
		return;
	// update overall hybridization partition function
	if (isHybridZ) {
		// add ED terms etc. before update
		overallZ += partZ * energy.getBoltzmannWeight( energy.getE(i1,j1,i2,j2,E_type(0)) );
		if (Z_isINF(overallZ))
			LOG(WARNING) <<"!!! PredictorMfeEns::overallZ = infinity (partition function overflow)" ;
	} else {
		// just update
		overallZ += partZ;
	}

	// store partial Z
	size_t maxLength = std::max(energy.getAccessibility1().getMaxLength(), energy.getAccessibility2().getMaxLength());
	size_t key = 0;
	key += i1;
	key += j1 * pow(maxLength, 1);
	key += i2 * pow(maxLength, 2);
	key += j2 * pow(maxLength, 3);
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

////////////////////////////////////////////////////////////////////////////


} // namespace
