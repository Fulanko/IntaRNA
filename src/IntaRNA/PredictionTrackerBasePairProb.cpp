
#include "IntaRNA/PredictionTrackerBasePairProb.h"

extern "C" {
	#include <ViennaRNA/vrna_config.h>
	#include <ViennaRNA/plotting/probabilities.h>
}

#include <boost/foreach.hpp>
#include <boost/regex.hpp>

namespace IntaRNA {

//////////////////////////////////////////////////////////////////////

PredictionTrackerBasePairProb::
PredictionTrackerBasePairProb(
		const InteractionEnergy & energy
		, const std::string & fileName
	)
 :	PredictionTracker()
	, energy(energy)
	, fileName(fileName)
	, probabilityThreshold(0.0001)
{
}

//////////////////////////////////////////////////////////////////////

PredictionTrackerBasePairProb::
~PredictionTrackerBasePairProb()
{
}

//////////////////////////////////////////////////////////////////////

void
PredictionTrackerBasePairProb::
updateOptimumCalled( const size_t i1, const size_t j1
					, const size_t i2, const size_t j2
					, const E_type curE
					)
{
}

//////////////////////////////////////////////////////////////////////

void
PredictionTrackerBasePairProb::
updateZ( PredictorMfeEns *predictor, SeedHandler *seedHandler )
{

	// sequence strings
	const std::string & rna1 = energy.getAccessibility1().getSequence().asString();
	const std::string & rna2 = energy.getAccessibility2().getSequence().asString();

	size_t s1 = energy.size1();
	size_t s2 = energy.size2();
	size_t n1 = energy.getAccessibility1().getMaxLength();
	size_t n2 = energy.getAccessibility2().getMaxLength();

	// calculate initial probabilities (for max window length)
	for (size_t i1 = 0; i1 < s1-n1+1; i1++) {
		for (size_t i2 = 0; i2 < s2-n2+1; i2++) {
			LOG(DEBUG) << i1 << ":" << i1 + n1 - 1 << ":" << i2 << ":" << i2 + n2 - 1 << " = " << predictor->getHybridZ(i1, i1 + n1 - 1, i2, i2 + n2 - 1);
			//if (!Z_equal(predictor->getHybridZ(i1, i1 + n1 - 1, i2, i2 + n2 - 1), 0)) {
				Interaction::Boundary key(i1, i1 + n1 -1, i2, i2 + n2 - 1);
				LOG(DEBUG) << ( predictor->getHybridZ(i1, i1 + n1 - 1, i2, i2 + n2 - 1) * energy.getBoltzmannWeight(energy.getED1(i1, i1 + n1 - 1) + energy.getED2(i2, i2 + n2 - 1)) ) / predictor->getZall();
				structureProbs[key] = ( predictor->getHybridZ(i1, i1 + n1 - 1, i2, i2 + n2 - 1) * energy.getBoltzmannWeight(energy.getED1(i1, i1 + n1 - 1) + energy.getED2(i2, i2 + n2 - 1)) ) / predictor->getZall();
			//}
		}
	}

	// loop over window lengths
	for (size_t w1 = n1; w1 > 0; w1--) {
		for (size_t w2 = n2; w2 > 0; w2--) {
			// skip initial probabilities
			if (w1 == n1 && w2 == n2) continue;
			// shift window over sequence length
			for (size_t i1 = 0; i1 < s1-w1+1; i1++) {
				for (size_t i2 = 0; i2 < s2-w2+1; i2++) {
					size_t j1 = i1 + w1 - 1;
					size_t j2 = i2 + w2 - 1;
					Z_type prob = 0.0;

					// if seed-based prediction, compute missing Z values
					/*if (seedHandler != NULL && i1 != j1 && i2 != j2 && Z_equal(predictor->getHybridZ(i1, j1, i2, j2), 0)) {
						LOG(DEBUG) << "missing Z at " << i1 << ":" << j1 << ":" << i2 << ":" << j2;
						computeMissingZ(i1, j1, i2, j2, predictor, seedHandler);
					} else {
						LOG(DEBUG) << "found Z at " << i1 << ":" << j1 << ":" << i2 << ":" << j2 << " with Z = " << getHybridZ(i1, j1, i2, j2, predictor);
					}*/

					LOG(DEBUG) << "window " << i1 << ":"  << j1 << ":"  << i2 << ":"  << j2;

					// TODO: check for case distinction
					// |||
					// | |
					// are not treated separately, resulting in
					// (1,1) having higher probability than (0,0) and (2,2)
					// in case of RNA sequences GGG and CCC

					// loop over combinations of (l..i), (j..r)
					/*
					for (size_t l1 = i1+1; l1-- > 0 && i1 - l1 < 2; ) {
						for (size_t l2 = i2+1; l2-- > 0 && i2 - l2 < 2; ) {
							for (size_t r1 = j1; r1 <= std::min(s1-1, j1+1); r1++) {
								for (size_t r2 = j2; r2 <= std::min(s2-1, j2+1); r2++) {
					*/
					for (size_t l1 = i1+1; l1-- > 0; ) {
						for (size_t l2 = i2+1; l2-- > 0; ) {
							for (size_t r1 = j1; r1 < s1; r1++) {
								// ensure maximal loop length
								if (r1-l1 > energy.getMaxInternalLoopSize1()+1) break;
								for (size_t r2 = j2; r2 < s2; r2++) {
									// ensure maximal loop length
									if (r2-l2 > energy.getMaxInternalLoopSize2()+1) break;

									// LOG(DEBUG) << "...check " << l1 << ":"  << r1 << ":"  << l2 << ":"  << r2;

									if (l1 == i1 && l2 == i2 && j1 == r1 && j2 == r2) {
										prob += ( predictor->getHybridZ(i1, j1, i2, j2) * energy.getBoltzmannWeight(energy.getED1(i1, j1) + energy.getED2(i2, j2)) ) / predictor->getZall();
									} else {
										// get outer probability
										Interaction::Boundary key(l1, r1, l2, r2);
										if ( structureProbs.find(key) != structureProbs.end()) {
											if (!Z_equal(predictor->getHybridZ(l1, r1, l2, r2), 0)) {
												prob += (( predictor->getHybridZ(l1, r1, l2, r2) * energy.getBoltzmannWeight(energy.getED1(l1, r1) + energy.getED2(l2, r2)) ) / predictor->getZall()
																 * (l1 == i1 && l2 == i2 ? 1 : energy.getBoltzmannWeight(energy.getE_interLeft(l1,i1,l2,i2)))
																 * getHybridZ(i1, j1, i2, j2, predictor)
																 * (j1 == r1 && j2 == r2 ? 1 : energy.getBoltzmannWeight(energy.getE_interLeft(j1,r1,j2,r2)))
															 ) / getHybridZ(l1, r1, l2, r2, predictor);
											}
										} else {
											throw std::runtime_error("Missing outer structure probability at: " + toString(l1) + ":" + toString(r1) + ":" + toString(l2) + ":" + toString(r2));
										}
									}
								} // r2
							} // r1
						} // l2
					} // l1

					// store structure probability
					Interaction::Boundary key(i1, j1, i2, j2);
 					structureProbs[key] = prob;

				} // i2
			} // i1
		} // w2
	} // w1

	// create plist
	struct vrna_elem_prob_s plist1[s1*s2+1];
	struct vrna_elem_prob_s plist2[s1*s2+1];
	size_t i = 0;

	for (auto it = structureProbs.begin(); it != structureProbs.end(); ++it)
	{
		if (it->first.i1 == it->first.j1 && it->first.i2 == it->first.j2 && it->second > probabilityThreshold) {
			LOG(DEBUG) << "prob at " << it->first.i1 << ":" << it->first.j1 << ":" << it->first.i2 << ":" << it->first.j2 << " = " << it->second;
			plist1[i].i = it->first.i1 + 1;
			plist1[i].j = s1 + it->first.i2 + 1;
			plist1[i].p = it->second;
			plist1[i].type = 0; // base-pair prob
			i++;
		}
	}

	// create dot plot
	std::string reverseRna2(rna2);
	std::reverse(reverseRna2.begin(), reverseRna2.end());
	char *rna = strdup((rna1 + reverseRna2).c_str());
	char *name = strdup(fileName.c_str());
	std::stringstream comment;
	comment << "Intermolecural Base-pair probabilities generated by "
	        << INTARNA_PACKAGE_STRING
          << " using Vienna RNA package "
					<< VRNA_VERSION;
	PS_dot_plot_list(rna, name, plist1, plist2, &(comment.str())[0u]);

}

////////////////////////////////////////////////////////////////////////////

void
PredictionTrackerBasePairProb::
computeMissingZ( const size_t i1, const size_t j1
							, const size_t i2, const size_t j2
							, PredictorMfeEns *predictor
							, SeedHandler* seedHandler )
{

	// search for known partition
	bool foundOuter = false;
	size_t l1, r1, l2, r2;
	for (l1 = i1+1; !foundOuter && l1-- > 0;) {
		for (l2 = i2+1; !foundOuter && l2-- > 0;) {
			for (r1 = j1-1; !foundOuter && r1++ < energy.size1()-1;) {
				for (r2 = j2-1; !foundOuter && r2++ < energy.size2()-1;) {
					if (!Z_equal(getHybridZ(l1, r1, l2, r2, predictor), 0)) {
						foundOuter = true;
					}
				}
			}
		}
	}

	if (!foundOuter) {
		throw std::runtime_error("Could not compute missing Z: no outer region found");
		return;
	}

	Z_type partZ = 0.0;

	// Check if full seed in subregion
	if (isFullSeedinRegion(i1, j1, i2, j2, seedHandler)) {

		// Case 1
		LOG(DEBUG) << "case1";

		// compute Z(l-i), Z(j-r)

		if (Z_equal(getHybridZ(l1, i1, l2, i2, predictor), 0)) {
			partZ = (
					Z_equal(getHybridZ(l1, r1, l2, r2, predictor), 0)
				/ Z_equal(getHybridZ(i1, r1, i2, r2, predictor), 0)
			) * energy.getBoltzmannWeight(energy.getED1(l1, i1) + energy.getED2(l2, i2));
		}
		updateHybridZ(l1, i1, l2, i2, partZ);

		if (Z_equal(getHybridZ(j1, r1, j2, r2, predictor), 0)) {
			partZ = (
					Z_equal(getHybridZ(l1, r1, l2, r2, predictor), 0)
				/ Z_equal(getHybridZ(l1, j1, l2, j2, predictor), 0)
			) * energy.getBoltzmannWeight(energy.getED1(j1, r1) + energy.getED2(j2, r2));
		}
		updateHybridZ(j1, r1, j2, r2, partZ);

	} else {

		// check if full seed outside of subregion (left or right)
		if (isFullSeedinRegion(l1, i1, l2, i2, seedHandler)) {
			// full seed left of subregion

			// Case 2.1 (left)
			LOG(DEBUG) << "case2.1 left";

			if (!Z_equal(getHybridZ(i1, j1, i2, j2, predictor), 0)) return;
			partZ = (
				  Z_equal(getHybridZ(l1, r1, l2, r2, predictor), 0)
				/ Z_equal(getHybridZ(l1, i1, l2, i2, predictor), 0)
			) * energy.getBoltzmannWeight(energy.getED1(i1, j1) + energy.getED2(i2, j2));
			updateHybridZ(i1, j1, i2, j2, partZ);

		} else if (isFullSeedinRegion(j1, r1, j2, r2, seedHandler)) {
			// full seed right of subregion

			// Case 2.1 (right)
			LOG(DEBUG) << "case2.1 right";

			if (!Z_equal(getHybridZ(i1, j1, i2, j2, predictor), 0)) return;
			partZ = (
				  Z_equal(getHybridZ(l1, r1, l2, r2, predictor), 0)
				/ Z_equal(getHybridZ(j1, r1, j2, r2, predictor), 0)
			) * energy.getBoltzmannWeight(energy.getED1(i1, j1) + energy.getED2(i2, j2));
			updateHybridZ(i1, j1, i2, j2, partZ);

		} else {

			// Case 2.2
			LOG(DEBUG) << "case2.2";

			// check side of k
			size_t k1, k2;
			bool kOnRight;
			if (i1 != l1 || i2 != l2) {
				// k at the left of region
				kOnRight = false;
				k1 = i1;
				k2 = i2;
			} else {
				kOnRight = true;
				k1 = j1;
				k2 = j2;
			}

			// find leftmost overlapping seeds
			std::vector< std::pair <size_t, size_t> > overlappingSeeds = getLeftMostSeedsAtK(k1, k2, seedHandler);

			// loop over overlapping seeds
			for(std::vector< std::pair <size_t, size_t> >::iterator it = overlappingSeeds.begin(); it != overlappingSeeds.end(); ++it) {
				// calculate missing partition function
				const size_t sl1 = seedHandler->getSeedLength1(it->first, it->second);
				const size_t sl2 = seedHandler->getSeedLength2(it->first, it->second);
				if (kOnRight) {
					partZ += (
						getHybridZ(i1, it->first+sl1-1, i2, it->second+sl2-1, predictor)
						* energy.getBoltzmannWeight(energy.getED1(i1, k1) + energy.getED2(i2, k2))
					) / energy.getBoltzmannWeight(getPartialSeedEnergy(it->first,it->second,k1,r1,k2,r2,seedHandler));
					//partZ *= energy.getBoltzmannWeight(energy.getE_init());
				} else {
					partZ += (
						getHybridZ(it->first, j1, it->second, j2, predictor)
						* energy.getBoltzmannWeight(energy.getED1(k1, j1) + energy.getED2(k2, j2))
					) / energy.getBoltzmannWeight(getPartialSeedEnergy(it->first,it->second,l1,k1,l2,k2,seedHandler));
				}
	    }

			if (kOnRight) {
				updateHybridZ(i1, k1, i2, k2, partZ);
			} else {
				updateHybridZ(k1, j1, k2, j2, partZ);
			}

			// TODO: if no overlapping seeds
			// what now? :|

		}

	}

	LOG(DEBUG) << "Z: " << partZ;

}

////////////////////////////////////////////////////////////////////////////

bool
PredictionTrackerBasePairProb::
isFullSeedinRegion( const size_t i1, const size_t j1
				        	, const size_t i2, const size_t j2
				        	, SeedHandler* seedHandler )
{
	size_t si1 = RnaSequence::lastPos, si2 = RnaSequence::lastPos;
	while( seedHandler->updateToNextSeed(si1,si2
			, i1, j1+1-seedHandler->getConstraint().getBasePairs()
			, i2, j2+1-seedHandler->getConstraint().getBasePairs()) )
	{
		const size_t sl1 = seedHandler->getSeedLength1(si1, si2);
		const size_t sl2 = seedHandler->getSeedLength2(si1, si2);
		const size_t sj1 = si1+sl1-1;
		const size_t sj2 = si2+sl2-1;
		if (sj1 <= j1 && sj2 <= j2) {
			return true;
		}
	}
	return false;
}

////////////////////////////////////////////////////////////////////////////

std::vector< std::pair <size_t, size_t> >
PredictionTrackerBasePairProb::
getLeftMostSeedsAtK( const size_t k1, const size_t k2
									 , SeedHandler* seedHandler )
{
  std::vector< std::pair <size_t, size_t> > seeds;
	size_t s1 = energy.size1();
	size_t s2 = energy.size2();
	size_t maxSeedLength = seedHandler->getConstraint().getBasePairs();

	size_t i1 = (k1 < maxSeedLength) ? 0 : k1 - maxSeedLength;
	size_t i2 = (k2 < maxSeedLength) ? 0 : k2 - maxSeedLength;
	size_t j1 = (k1 + maxSeedLength >= s1) ? s1-1 : k1 + maxSeedLength;
	size_t j2 = (k2 + maxSeedLength >= s2) ? s2-1 : k2 + maxSeedLength;

	size_t si1 = RnaSequence::lastPos, si2 = RnaSequence::lastPos;
	while( seedHandler->updateToNextSeed(si1, si2, i1, j1, i2, j2) ) {
		// check if a loop-overlapping seed already exists
		bool found = false;
		for(std::vector< std::pair <size_t, size_t> >::iterator it = seeds.begin(); it != seeds.end(); ++it) {
			if (seedHandler->areLoopOverlapping(si1, si2, it->first, it->second)) {
				found = true;
				break;
			}
    }
		if (!found) {
			seeds.push_back(std::make_pair(si1, si2));
		}
	}
	return seeds;
}

////////////////////////////////////////////////////////////////////////////

E_type
PredictionTrackerBasePairProb::
getPartialSeedEnergy( const size_t si1, const size_t si2
										, const size_t i1, const size_t j1
										, const size_t i2, const size_t j2
										, SeedHandler* seedHandler )
{
	// trace S
	Interaction interaction = Interaction(energy.getAccessibility1().getSequence(), energy.getAccessibility2().getAccessibilityOrigin().getSequence());
	const size_t sl1 = seedHandler->getSeedLength1(si1, si2);
	const size_t sl2 = seedHandler->getSeedLength2(si1, si2);
	interaction.basePairs.push_back( energy.getBasePair(si1, si2) );
	seedHandler->traceBackSeed( interaction, si1, si2 );
	interaction.basePairs.push_back( energy.getBasePair(si1+sl1-1, si2+sl2-1) );

	E_type partE = 0;
	size_t i1old = i1;
	size_t i2old = i2;

	for (size_t i = 0; i < interaction.basePairs.size(); i++) {
		// get index of current base pair
		size_t s1 = energy.getIndex1(interaction.basePairs[i]);
		size_t s2 = energy.getIndex2(interaction.basePairs[i]);

		if (s1 < i1 || s2 < i2) continue;
		if (s1 > j1 && s2 < j2) break;

		if (!E_isINF(energy.getE_interLeft(i1old,s1,i2old,s2))) {
			// add hybridization energy
			partE += energy.getE_interLeft(i1old,s1,i2old,s2);
			// store
			i1old = s1;
			i2old = s2;
		}
	}
	return partE;
}

////////////////////////////////////////////////////////////////////////////

Z_type
PredictionTrackerBasePairProb::
getHybridZ( const size_t i1, const size_t j1
					, const size_t i2, const size_t j2
					, PredictorMfeEns *predictor)
{
	Z_type partZ = predictor->getHybridZ(i1, j1, i2, j2);
	if (Z_equal(partZ, 0)) {
		Interaction::Boundary key(i1, j1, i2, j2);
		if ( Z_partition.find(key) == Z_partition.end() ) {
			partZ = Z_type(0);
		} else {
			partZ = Z_partition[key];
		}
	}
	return partZ;
}

////////////////////////////////////////////////////////////////////////////

void
PredictionTrackerBasePairProb::
updateHybridZ( const size_t i1, const size_t j1
						 , const size_t i2, const size_t j2
						 , const Z_type partZ )
{
	Interaction::Boundary key(i1,j1,i2,j2);
	auto keyEntry = Z_partition.find(key);
	keyEntry->second = partZ;
}

////////////////////////////////////////////////////////////////////////////

} // namespace
