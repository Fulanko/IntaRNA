
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
}

////////////////////////////////////////////////////////////////////////////


} // namespace
