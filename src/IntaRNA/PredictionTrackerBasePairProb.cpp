
#include "PredictionTrackerBasePairProb.h"

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
	, maxDotPlotSize(640)
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
	PredictorMfeEns2dSeedExtension* seedPredictor = dynamic_cast<PredictorMfeEns2dSeedExtension*>(predictor);
	isSeedPredictor = (seedPredictor != nullptr);

	// sequence strings
	const std::string & rna1 = energy.getAccessibility1().getSequence().asString();
	const std::string & reverseRna2 = energy.getAccessibility2().getAccessibilityOrigin().getSequence().asString();

	size_t s1 = energy.size1();
	size_t s2 = energy.size2();
	size_t n1 = energy.getAccessibility1().getMaxLength();
	size_t n2 = energy.getAccessibility2().getMaxLength();

	Z_type maxZ = 0.0;
	Interaction::Boundary interactionBoundary;

	// initialize Z_partition
	const PredictorMfeEns::Site2Z_hash & Z_partition = predictor->getZPartition();

	// create index of left/right boundaries
	for (auto z = Z_partition.begin(); z != Z_partition.end(); ++z) {
		// identify best interaction boundary
		Z_type Zstruct = z->second * energy.getBoltzmannWeight(energy.getE(z->first.i1, z->first.j1, z->first.i2, z->first.j2, E_type(0)));
		if (Zstruct > maxZ) {
			maxZ = Zstruct;
			interactionBoundary = z->first;
		}

		Interaction::BasePair iBP(z->first.i1, z->first.i2);
		Interaction::BasePair jBP(z->first.j1, z->first.j2);

		// create left and jBP index
		if (z->first.i1 != z->first.j1 && z->first.i2 != z->first.j2) {
			// encode iBP/jBP boundary
			// create left index
			rightExt[iBP].insert(jBP);
			// create right index
			if (seedHandler != NULL) {
				leftExt[jBP].insert(iBP);
			}
		}

	} // it (Z_partition)

	// Compute base-pair probabilities via combinations
	if (!isSeedPredictor) {
		computeBasePairProbsNoSeed(predictor);
	} else {
		computeBasePairProbs(seedPredictor, seedHandler);
	}

	// build plist
	struct vrna_elem_prob_s plist[structureProbs.size()+1];
	size_t i = 0;
	const Z_type Zall = predictor->getZall();
	for (auto sp = structureProbs.begin(); sp != structureProbs.end(); ++sp) {
		if ( (sp->second /Zall)  > probabilityThreshold) {
			plist[i].i = sp->first.first + 1;
			plist[i].j = sp->first.second + 1;
			plist[i].p = (sp->second /Zall);
			plist[i].type = 0; // base-pair prob
			i++;
		}
	}
	plist[i].i = 0; // list end

	// create dot plot
	char *name = strdup(fileName.c_str());
	std::string comment =
	  "Intermolecular base-pair probabilities generated by "
	  INTARNA_PACKAGE_STRING
    " using Vienna RNA package "
    VRNA_VERSION;

	generateDotPlotSvg(strdup(rna1.c_str()), strdup(reverseRna2.c_str()), name, plist, comment.c_str(), interactionBoundary, energy);

}

////////////////////////////////////////////////////////////////////////////

void
PredictionTrackerBasePairProb::
computeBasePairProbsNoSeed( const PredictorMfeEns *predictor )
{
	for (auto z = predictor->getZPartition().begin(); z != predictor->getZPartition().end(); ++z) {
		assert( !Z_equal(z->second, 0) );

		Z_type bpProb = z->second * energy.getBoltzmannWeight(energy.getE(z->first.i1, z->first.j1, z->first.i2, z->first.j2, E_type(0)));

		Interaction::BasePair iBP(z->first.i1, z->first.i2);
		Interaction::BasePair jBP(z->first.j1, z->first.j2);

		// left end and single bp
		updateProb(iBP, bpProb);

		if (iBP < jBP) {
			// right end
			updateProb(jBP, bpProb);

			// inner (rightExt > jBP)
			for (auto right = rightExt[jBP].begin(); right != rightExt[jBP].end(); ++right) {
				Z_type innerProb = z->second * getHybridZ(z->first.j1, right->first, z->first.j2, right->second, predictor)
									 * energy.getBoltzmannWeight(energy.getE(z->first.i1, right->first, z->first.i2, right->second, -energy.getE_init()));
				updateProb(jBP, innerProb);
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////

void
PredictionTrackerBasePairProb::
computeBasePairProbs( const PredictorMfeEns2dSeedExtension *predictor, const SeedHandler* seedHandler )
{
	auto ZL_partition = predictor->getZLPartition();
	size_t i1, j1, i2, j2;
	for (auto z = predictor->getZPartition().begin(); z != predictor->getZPartition().end(); ++z) {
		i1 = z->first.i1;
		j1 = z->first.j1;
		i2 = z->first.i2;
		j2 = z->first.j2;

		Z_type bpZ;
		// i external left
		bpZ = getHybridZ(i1, j1, i2, j2, predictor)
				* energy.getBoltzmannWeight(energy.getE(i1, j1, i2, j2, E_type(0)));
		updateProb(Interaction::BasePair(i1,i2), bpZ);

		// j external right (and not single bp)
		if (i1!=j1) {
			updateProb(Interaction::BasePair(j1,j2), bpZ);
		}

		Z_type tempZ;

		// loop internal k (inbetween i and j)
		for (size_t k1 = i1+1; k1 < j1; k1++) {
			for (size_t k2 = i2+1; k2 < j2; k2++) {

				Interaction::BasePair kBP(k1, k2);
				bpZ = 0;
				tempZ = 0;

				// k internal
				Z_type ZSleft = getHybridZ(i1, k1, i2, k2, predictor);
				Z_type ZSright = getHybridZ(k1, j1, k2, j2, predictor);

				// ... ZS:ZS
				tempZ = ZSleft * ZSright / energy.getBoltzmannWeight(energy.getE_init());
				bpZ += tempZ;

				// ... ZS:ZN
				tempZ = ZSleft * getZRPartition(predictor, seedHandler, k1, j1, k2, j2);
				bpZ += tempZ;

				// ... ZN:ZS
				tempZ = getZPartitionValue(&ZL_partition, Interaction::Boundary(i1,k1,i2,k2), false) * ZSright / energy.getBoltzmannWeight(energy.getE_init());
				bpZ += tempZ;

				// seeds overlapping k
				size_t si1 = RnaSequence::lastPos, si2 = RnaSequence::lastPos;
				if (seedHandler->getConstraint().getBasePairs() > 2) {
					while( seedHandler->updateToNextSeedWithK(si1, si2, k1, k2, false))
					{
						size_t sj1 = si1+seedHandler->getSeedLength1(si1,si2)-1;
						size_t sj2 = si2+seedHandler->getSeedLength2(si1,si2)-1;
						// seed region mismatch
						if ((sj1 == j1 && sj2 != j2) || (sj1 != j1 && sj2 == j2)) {
							continue;
						}
						// check if still in region
						if (si1 < i1 || si2 < i2 || sj1 > j1 || sj2 > j2) {
							continue;
						}
						// ZN:seed:ZS
						tempZ = getZPartitionValue(&ZL_partition, Interaction::Boundary(i1,si1,i2,si2), false) // contains E_init
									* energy.getBoltzmannWeight(seedHandler->getSeedE(si1, si2))
									* getHybridZ(sj1, j1, sj2, j2, predictor)
									/ energy.getBoltzmannWeight(energy.getE_init());
						bpZ += tempZ;

						// ZS:seed:ZN
						tempZ = getHybridZ(i1, si1, i2, si2, predictor)
									* energy.getBoltzmannWeight(seedHandler->getSeedE(si1, si2))
									* getZRPartition(predictor, seedHandler, sj1, j1, sj2, j2);
						bpZ += tempZ;

						// ZN:seed:ZP
						tempZ = getZPartitionValue(&ZL_partition, Interaction::Boundary(i1,si1,i2,si2), false) // contains E_init
									* energy.getBoltzmannWeight(seedHandler->getSeedE(si1, si2))
									* (getZHPartition(predictor, seedHandler, sj1, j1, sj2, j2) - getHybridZ(sj1, j1, sj2, j2, predictor))
									/ energy.getBoltzmannWeight(energy.getE_init());
						bpZ += tempZ;

						// ZNL:seed':ZNR
						size_t spi1 = RnaSequence::lastPos, spi2 = RnaSequence::lastPos;
						while( seedHandler->updateToNextSeedWithK(spi1,spi2,k1,k2))
						{
							if (seedHandler->areLoopOverlapping(spi1, spi2, si1, si2)) {
								size_t spj1 = spi1+seedHandler->getSeedLength1(spi1,spi2)-1;
								size_t spj2 = spi2+seedHandler->getSeedLength2(spi1,spi2)-1;
								if (k1 < spj1 && k2 < spj2) {
									continue;
								}
								// check if still in region
								if (spi1 < i1 || spi2 < i2 ||spj1 > j1 || spj2 > j2) {
									continue;
								}
								tempZ = getZPartitionValue(&ZL_partition, Interaction::Boundary(i1,spi1,i2,spi2), false) // contains E_init
											* energy.getBoltzmannWeight(PredictorMfeEns2dSeedExtension::getNonOverlappingEnergy(spi1, spi2, si1, si2, energy, *seedHandler, false))
											* energy.getBoltzmannWeight(seedHandler->getSeedE(si1, si2))
											* getZRPartition(predictor, seedHandler, sj1, j1, sj2, j2);
								bpZ += tempZ;
							}
						}
						
					}
				}

				// add ED values, dangling ends, etc.
				bpZ *= energy.getBoltzmannWeight(energy.getE(i1, j1, i2, j2, E_type(0)));
				updateProb(kBP, bpZ);

			} // k2
		} // k1
	} // Z_partitions
}

////////////////////////////////////////////////////////////////////////////

Z_type
PredictionTrackerBasePairProb::
getZPartitionValue( const Site2Z_hash *Zpartition, const Interaction::Boundary & boundary, const bool addZInit )
{
	auto keyEntry = Zpartition->find(boundary);
	if ( Zpartition->find(boundary) == Zpartition->end() ) {
		return 0;
	} else {
		if (addZInit) {
      return keyEntry->second * energy.getBoltzmannWeight(energy.getE_init());
		} else {
			return keyEntry->second;
		}
	}
}

////////////////////////////////////////////////////////////////////////////

Z_type
PredictionTrackerBasePairProb::
getZHPartition( const PredictorMfeEns2dSeedExtension *predictor, const SeedHandler* seedHandler
              , const size_t i1, const size_t j1
	            , const size_t i2, const size_t j2 )
{
	// sanity check
	if (!energy.areComplementary(i1,i2) || !energy.areComplementary(j1,j2)) {
		return Z_type(0);
	}
	// single bp boundary with ZH == 1
	if (i1==j1 && i2==j2) {
		return Z_type(energy.getBoltzmannWeight(energy.getE_init()));
	}
	Interaction::Boundary boundary(i1,j1,i2,j2);

	// check if ZH available from predictor
	Z_type ZH = getZPartitionValue(&predictor->getZHPartition(), boundary, true);
	if (Z_equal(ZH, 0.0)) {
		// check if ZH available in missing ZH partitions
		auto ZHentry = ZH_partition_missing.find(boundary);
		if ( ZHentry != ZH_partition_missing.end() ) {
			ZH = ZHentry->second;
		} else {
			ZH = fillHybridZ(i1, j1, i2, j2, seedHandler);
		}
	}

	return ZH;
}

////////////////////////////////////////////////////////////////////////////

Z_type
PredictionTrackerBasePairProb::
getZRPartition( const PredictorMfeEns2dSeedExtension *predictor, const SeedHandler* seedHandler
              , const size_t i1, const size_t j1
	            , const size_t i2, const size_t j2 )
{
	// memoization
	Interaction::Boundary boundary(i1,j1,i2,j2);
	auto keyEntry = ZR_partition.find(boundary);
	if ( ZR_partition.find(boundary) != ZR_partition.end() ) {
		return ZR_partition[boundary];
	}

	// single bp boundary with ZNR == 1
	if (i1==j1 && i2==j2 && energy.areComplementary(i1,i2)) {
		return Z_type(1);
	}

	Z_type partZ = getZHPartition(predictor, seedHandler, i1, j1, i2, j2);
	Z_type ZS = getHybridZ(boundary, predictor);
	assert( ZS <= partZ );
	partZ -= ZS;

	// remove Einit
	partZ /= energy.getBoltzmannWeight(energy.getE_init());

	// iterate all seeds that overlap anchor seed sa on the right side
  size_t si1 = RnaSequence::lastPos, si2 = RnaSequence::lastPos;
	while( seedHandler->updateToNextSeedWithK(si1, si2, i1, i2, false))
	{
		size_t sj1 = si1 + seedHandler->getSeedLength1(si1, si2) - 1;
		size_t sj2 = si2 + seedHandler->getSeedLength2(si1, si2) - 1;
		// check if still in region
		if ( j1 < sj1 || j2 < sj2 ) continue;
		E_type Eoverlap = seedHandler->getSeedE(si1, si2)
						- PredictorMfeEns2dSeedExtension::getNonOverlappingEnergy(si1, si2, i1, i2, energy, *seedHandler, false);
		Z_type corrZterm = energy.getBoltzmannWeight(Eoverlap) * getZRPartition(predictor, seedHandler, sj1, j1, sj2, j2);
		assert(corrZterm <= partZ);
		partZ -= corrZterm;
	}

	// store ZR_partition
	ZR_partition[boundary] = partZ;

	return partZ;
}

////////////////////////////////////////////////////////////////////////////

Z_type
PredictionTrackerBasePairProb::
getHybridZ( const size_t i1, const size_t j1
					, const size_t i2, const size_t j2
					, const PredictorMfeEns *predictor)
{
	Interaction::Boundary boundary(i1, j1, i2, j2);
	return getHybridZ(boundary, predictor);
}

////////////////////////////////////////////////////////////////////////////

Z_type
PredictionTrackerBasePairProb::
getHybridZ( const Interaction::Boundary & boundary
					, const PredictorMfeEns *predictor)
{
	// check in original data
	if ( predictor->getZPartition().find(boundary) != predictor->getZPartition().end() ) {
		return predictor->getZPartition().find(boundary)->second;
	} else {
		// fall back
		return Z_type(0);
	}
}

////////////////////////////////////////////////////////////////////////////

Z_type
PredictionTrackerBasePairProb::
getBasePairProb( const size_t i1, const size_t i2
					     , const PredictorMfeEns *predictor)
{
	Interaction::BasePair bp(i1, i2);
	if ( structureProbs.find(bp) == structureProbs.end() ) {
		return Z_type(0);
	} else {
		return structureProbs[bp] / predictor->getZall();
	}
}

////////////////////////////////////////////////////////////////////////////

Z_type
PredictionTrackerBasePairProb::
fillHybridZ( const size_t ll1, const size_t si1, const size_t ll2, const size_t si2, const SeedHandler* seedHandler )
{
#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!energy.areComplementary(si1,si2) )
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::fillHybridZ("+toString(si1)+","+toString(si2)+",..) are not complementary");
#endif

  hybridZ.resize( si1-ll1+1, si2-ll2+1 );

	// global vars to avoid reallocation
	size_t i1,i2,k1,k2;

	// determine whether or not lonely base pairs are allowed or if we have to
	// ensure a stacking to the right of the left boundary (i1,i2)
	const size_t noLpShift = seedHandler->getConstraint().isLpAllowed() ? 0 : 1;
	Z_type iStackZ = Z_type(1);

	Z_type overallZ = Z_type(0);

	// iterate over all window starts i1 (seq1) and i2 (seq2)
	for (size_t l1=0; l1 < hybridZ.size1(); l1++) {
		for (size_t l2=0; l2 < hybridZ.size2(); l2++) {
			i1 = si1-l1;
			i2 = si2-l2;

			// referencing cell access
			Z_type & curZ = hybridZ(si1-i1,si2-i2);

			// init current cell (0 if not just right-most (j1,j2) base pair)
			curZ = (i1==si1 && i2==si2) ? energy.getBoltzmannWeight(energy.getE_init()) : 0.0;

			// check if complementary (use global sequence indexing)
			if( i1<si1 && i2<si2 && energy.areComplementary(i1,i2) ) {

				// right-stacking of i if no-LP
				if (!seedHandler->getConstraint().isLpAllowed()) {
					// skip if no stacking possible
					if (!energy.areComplementary(i1+noLpShift,i2+noLpShift)) {
						continue;
					}
					// get stacking energy to avoid recomputation in recursion below
					iStackZ = energy.getBoltzmannWeight(energy.getE_interLeft(i1,i1+noLpShift,i2,i2+noLpShift));
					// check just stacked
					curZ += iStackZ + hybridZ(l1-noLpShift,l2-noLpShift);
				}

				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				for (k1=i1+noLpShift; k1++ < si1; ) {
					// ensure maximal loop length
					if (k1-i1-noLpShift > energy.getMaxInternalLoopSize1()+1) break;
					for (k2=i2+noLpShift; k2++ < si2; ) {
						// ensure maximal loop length
						if (k2-i2-noLpShift > energy.getMaxInternalLoopSize2()+1) break;
						// check if (k1,k2) are valid left boundary
						if ( ! Z_equal(hybridZ(si1-k1,si2-k2), 0.0) ) {
							curZ += (iStackZ
									* energy.getBoltzmannWeight(energy.getE_interLeft(i1+noLpShift,k1,i2+noLpShift,k2))
									* hybridZ(si1-k1,si2-k2));
						}
					} // k2
				} // k1
			} // complementary

			// store partial Z
			Interaction::Boundary key(i1,si1,i2,si2);
			auto keyEntry = ZH_partition_missing.find(key);
			if ( ZH_partition_missing.find(key) == ZH_partition_missing.end() ) {
				ZH_partition_missing[key] = curZ;
			}
			// store Z of full region for return value
			if (i1 == ll1 && i2 == ll2) {
				overallZ = curZ;
			}

		} // i2
	} // i1

	return overallZ;
}

////////////////////////////////////////////////////////////////////////////

bool
PredictionTrackerBasePairProb::
generateDotPlot( const char *seq1, const char *seq2, const char *fileName
							 , const plist *pl, const char *comment
						   , const Interaction::Boundary interactionBoundary )
{
	FILE *file;
	file = fopen(fileName,"w");
	if (file == NULL) return false; /* failure */

	int bbox[4];
	bbox[0] = 0;
	bbox[1] = 0;
	bbox[2] = 72 * (strlen(seq1) + 3);
	bbox[3] = 72 * (strlen(seq2) + 3);

	size_t maxSize = std::max(bbox[2], bbox[3]);
	float scale = 1;
	if (maxSize > maxDotPlotSize) {
		scale = maxDotPlotSize / (float)maxSize;
		bbox[2] *= scale;
		bbox[3] *= scale;
	}

	fprintf(file,
          "%%!PS-Adobe-3.0 EPSF-3.0\n"
          "%%%%Creator: IntaRNA\n"
          "%%%%Title: RNA Dot Plot\n"
          "%%%%BoundingBox: %d %d %d %d\n"
          "%%%%DocumentFonts: Helvetica\n"
          "%%%%Pages: 1\n"
          "%%%%EndComments\n",
          bbox[0], bbox[1], bbox[2], bbox[3]);

	// scaling
	fprintf(file,
		      "%%%%BeginProcSet: epsffit 1 0\n"
		      "gsave\n"
		      "%f 0 translate\n"
	        "%f %f scale\n"
		      "%%%%EndProcSet\n\n",
				  scale, scale, scale);

  // comment

	if (comment) {
    fprintf(file, "%%%% %s\n", comment);
  }

	fprintf(file, "/DPdict 100 dict def\n");
	fprintf(file, "DPdict begin\n");

	// ps template

	fprintf(file, dotplotTemplate);
  fprintf(file, "end\n");
	fprintf(file, "DPdict begin\n");

	// sequences

	unsigned int i, length;
  length = strlen(seq1);
  fprintf(file, "/sequence1 { (\\\n");
  i = 0;
  while (i < length) {
    fprintf(file, "%.255s\\\n", seq1 + i);  /* no lines longer than 255 */
    i += 255;
  }
  fprintf(file, ") } def\n");
	fprintf(file, "/len { sequence1 length } bind def\n\n");
	length = strlen(seq2);
  fprintf(file, "/sequence2 { (\\\n");
  i = 0;
  while (i < length) {
    fprintf(file, "%.255s\\\n", seq2 + i);  /* no lines longer than 255 */
    i += 255;
  }
  fprintf(file, ") } def\n");
	fprintf(file, "/len2 { sequence2 length } bind def\n\n");

	fprintf(file, "72 72 translate\n"
	               "72 72 scale\n");

  fprintf(file, "/Helvetica findfont 0.95 scalefont setfont\n\n");

	// basepair data

	fprintf(file,"drawseq1\n");
	fprintf(file,"drawseq2\n");

	fprintf(file,"%%data starts here\n");

	fprintf(file,"%%start of base pair probability data\n");

  fprintf(file, "/coor [\n");

	for (const plist *pl1 = pl; pl1->i > 0; pl1++) {
    if (pl1->type == 0) {
      fprintf(file, "%1.9f %d %d boxgray\n", sqrt(pl1->p), pl1->i, pl1->j);
    }
  }

  fprintf(file, "] def\n");

	fprintf(file, "0.25 0.25 0.25 setrgbcolor\n");

	fprintf(file, "\n%%draw the grid\ndrawgrid\n\n");

	// print frame
	fprintf(file,
		      "0.03 setlinewidth\n\
           %1.1f %1.1f %zu %zu rectangle\n\
					 0 0 0 setrgbcolor\n\
           stroke\n", 0.5, 0.5, strlen(seq1), strlen(seq2));

	// print best interaction outline
	fprintf(file,
		      "0.03 setlinewidth\n\
           %1.1f %1.1f %zu %zu rectangle\n\
					 1 0 0 setrgbcolor\n\
           stroke\n", (float)interactionBoundary.i1 + 0.5, (float)interactionBoundary.i2 + 0.5, interactionBoundary.j1 - interactionBoundary.i1 + 1, interactionBoundary.j2 - interactionBoundary.i2 + 1);

	fprintf(file, "showpage\n"
	            "end\n"
              "%%%%EOF\n");

	fclose(file);
	return true; /* success */
}

////////////////////////////////////////////////////////////////////////////

bool
PredictionTrackerBasePairProb::
generateDotPlotSvg( const char *seq1, const char *seq2, const char *fileName
							 , const plist *pl, const char *comment
						   , const Interaction::Boundary interactionBoundary
							 , const InteractionEnergy & energy )
{
	FILE *file;
	file = fopen(fileName, "w");
	if (file == NULL) return false; /* failure */

	// file information
	fprintf(file, "<!--\n\
%s \n\n\
Sequence 1: %s\n\
Sequence 2: %s\n"
	, comment, seq1, seq2);

	fprintf(file, "-->\n\n");

  const size_t unitSize = 1;
	const size_t boxSize = 2 * unitSize;
	const size_t maxWidth = (2.5+strlen(seq1)) * boxSize;
	const size_t maxHeight = (2.5+strlen(seq2)) * boxSize;

	fprintf(file,
          "<svg viewBox='0 0 %d %d' width='%d' height='%d' version='1.1' xmlns='http://www.w3.org/2000/svg'>\n",
					maxWidth+boxSize, maxHeight+boxSize, maxDotPlotSize, maxDotPlotSize);

	// draw dots
	fprintf(file, "\n<!-- Base pair probabilities -->\n\n");
	for (const plist *pl1 = pl; pl1->i > 0; pl1++) {
    if (pl1->type == 0) {
			std::ostringstream message;
			message << "(" << (pl1->i) << "," << (strlen(seq2)-pl1->j+1) << ") = " << pl1->p;
      fprintf(file, drawSvgSquare(pl1->i, strlen(seq2)-pl1->j+1, boxSize, pl1->p, "bp", message.str().c_str()));
    }
  }

	fprintf(file, "\n<!-- Unpaired probabilities -->\n\n");
	// unpaired probs seq1
	for (size_t i = 0; seq1[i] != '\0'; i++) {
		float acc = 1-energy.getBoltzmannWeight(energy.getED1(i, i));
		if (acc > 0) {
      std::ostringstream message;
		  message << "(" << (i+1) << ") = " << acc;
		  fprintf(file, drawSvgSquare(i+1, -0.5, boxSize, acc, "unpaired", message.str().c_str()));
		}
  }

	// unpaired probs seq2
	for (size_t i = 0; seq2[i] != '\0'; i++) {
		float acc = 1-energy.getBoltzmannWeight(energy.getED2(i, i));
		if (acc > 0) {
      std::ostringstream message;
		  message << "(" << (strlen(seq2)-i) << ") = " << acc;
		  fprintf(file, drawSvgSquare(-0.5, strlen(seq2)-i, boxSize, acc, "unpaired", message.str().c_str()));
		}
  }

	// draw grid
	const float strokeWidth = 0.02 * boxSize;
	fprintf(file, "\n<!-- Grid -->\n\n");
	for (size_t i = 0; i <= strlen(seq1); i++) {
		if (i % 5 == 0) {
      fprintf(file,
		          "<line x1='%d' y1='%d' x2='%d' y2='%d' class='sep5' style='stroke:#555;stroke-width:%f;opacity:0.2;'/>\n"
					  	, size_t((i+2.5)*2*unitSize), 5*unitSize, size_t((i+2.5)*2*unitSize), maxHeight, strokeWidth);
		}
		if (i % 10 == 0) {
      fprintf(file,
		          "<line x1='%d' y1='%d' x2='%d' y2='%d' class='sep10' style='stroke:#555;stroke-width:%f;opacity:0.5;'/>\n"
					  	, size_t((i+2.5)*2*unitSize), 5*unitSize, size_t((i+2.5)*2*unitSize), maxHeight, strokeWidth);
		}
		if (i % 50 == 0) {
      fprintf(file,
		          "<line x1='%d' y1='%d' x2='%d' y2='%d' class='sep50' style='stroke:#555;stroke-width:%f;opacity:1;'/>\n"
					  	, size_t((i+2.5)*2*unitSize), 5*unitSize, size_t((i+2.5)*2*unitSize), maxHeight, strokeWidth);
		}
	}
	for (size_t i = 0; i <= strlen(seq2); i++) {
		if (i % 5 == 0) {
      fprintf(file,
		          "<line x1='%d' y1='%d' x2='%d' y2='%d' class='sep5' style='stroke:#555;stroke-width:%f;opacity:0.2;'/>\n"
						  , 5*unitSize, size_t((i+2.5)*2*unitSize), maxWidth, size_t((i+2.5)*2*unitSize), strokeWidth);
		}
		if (i % 10 == 0) {
      fprintf(file,
		          "<line x1='%d' y1='%d' x2='%d' y2='%d' class='sep10' style='stroke:#555;stroke-width:%f;opacity:0.5;'/>\n"
						  , 5*unitSize, size_t((i+2.5)*2*unitSize), maxWidth, size_t((i+2.5)*2*unitSize), strokeWidth);
		}
		if (i % 50 == 0) {
      fprintf(file,
		          "<line x1='%d' y1='%d' x2='%d' y2='%d' class='sep50' style='stroke:#555;stroke-width:%f;opacity:1;'/>\n"
						  , 5*unitSize, size_t((i+2.5)*2*unitSize), maxWidth, size_t((i+2.5)*2*unitSize), strokeWidth);
		}
	}

	// draw frames
	fprintf(file,
		      "<rect x='%d' y='%d' width='%d' height='%d' class='frame' style='pointer-events:none;stroke:black;stroke-width:%f;fill-opacity:0;'/>\n"
					, 5*unitSize, 5*unitSize, maxWidth-5*unitSize, (maxHeight-5*unitSize), 2*strokeWidth);
	fprintf(file,
		      "<rect x='%d' y='%d' width='%d' height='%d' class='frame' style='pointer-events:none;stroke:black;stroke-width:%f;fill-opacity:0;'/>\n"
					, 5*unitSize, boxSize, maxWidth-5*unitSize, boxSize, 2*strokeWidth);
	fprintf(file,
		      "<rect x='%d' y='%d' width='%d' height='%d' class='frame' style='pointer-events:none;stroke:black;stroke-width:%f;fill-opacity:0;'/>\n"
					, boxSize, 5*unitSize, boxSize, maxHeight-5*unitSize, 2*strokeWidth);

	fprintf(file, "\n<!-- Sequences -->\n\n");
	// draw sequence 1
	for (size_t i = 0; seq1[i] != '\0'; i++) {
    fprintf(file,
		        "<text x='%d' y='%d' class='nt' style='text-anchor:middle;font-size:%dpx;font-family:Arial;'>%c<title>(%d)</title></text>\n"
						, (i+3) * boxSize, boxSize, boxSize, seq1[i], i+1);
		fprintf(file,
		        "<text x='%d' y='%d' class='nt' style='text-anchor:middle;font-size:%dpx;font-family:Arial;'>%c<title>(%d)</title></text>\n"
						, (i+3) * boxSize, maxHeight + boxSize, boxSize, seq1[i], i+1);
  }

	// draw sequence 2
	for (size_t i = 0; seq2[i] != '\0'; i++) {
    fprintf(file,
		        "<text x='%d' y='%d' class='nt' style='text-anchor:middle;font-size:%dpx;font-family:Arial;'>%c<title>(%d)</title></text>\n"
						, unitSize, size_t((i+3.5)*2*unitSize), boxSize, seq2[i], i+1);
		fprintf(file,
		        "<text x='%d' y='%d' class='nt' style='text-anchor:middle;font-size:%dpx;font-family:Arial;'>%c<title>(%d)</title></text>\n"
						, maxWidth + unitSize, size_t((i+3.5)*2*unitSize), boxSize, seq2[i], i+1);
  }

	// draw best interaction outline
	fprintf(file, "\n<!-- Mfe -->\n\n");
	fprintf(file,
		      "<rect x='%d' y='%d' width='%d' height='%d' class='mfe' style='pointer-events:none;stroke:red;stroke-width:%f;fill-opacity:0;'/>\n"
					, size_t((interactionBoundary.i1+2.5)*2*unitSize)
					, size_t((strlen(seq2)-interactionBoundary.j2+1.5)*2*unitSize)
					, (interactionBoundary.j1-interactionBoundary.i1+1) * boxSize
					, (interactionBoundary.j2-interactionBoundary.i2+1) * boxSize
					, 2*strokeWidth);

	fprintf(file, "\n<!-- Styling -->\n\n");
  fprintf(file, "<style type='text/css'>\n\
    <![CDATA[\n\
      .nt {\n\
				text-anchor:middle !important;\n\
				font-size:%dpx !important;\n\
				font-family:Arial !important;\n\
			}\n\
      .bp {\n\
				fill:blue !important;\n\
				width:%dpx !important;\n\
				height:%dpx !important;\n\
			}\n\
			.unpaired {\n\
				fill:red !important;\n\
				width:%dpx !important;\n\
				height:%dpx !important;\n\
			}\n\
      .sep5 {\n\
				stroke:#555 !important;\n\
				stroke-width:%f !important;\n\
				opacity:0.2 !important;\n\
			}\n\
      .sep10 {\n\
				stroke:#555 !important;\n\
				stroke-width:%f !important;\n\
				opacity:0.5 !important;\n\
			}\n\
      .sep50 {\n\
				stroke:#555 !important;\n\
				stroke-width:%f !important;\n\
				opacity:1 !important;\n\
			}\n\
      .frame {\n\
				pointer-events:none !important;\n\
				stroke:black !important;\n\
				stroke-width:%f !important;\n\
				fill-opacity:0 !important;\n\
			}\n\
			.mfe {\n\
				pointer-events:none !important;\n\
				stroke:blue !important;\n\
				stroke-width:%f !important;\n\
				fill-opacity:0 !important;\n\
			}\n\
    ]]>\n\
  </style>\n", boxSize, boxSize, boxSize, boxSize, boxSize, strokeWidth, strokeWidth, strokeWidth, 2*strokeWidth, 4*strokeWidth);
	fprintf(file, "</svg>\n");

	fclose(file);
	return true; /* success */
}

////////////////////////////////////////////////////////////////////////////

const char*
PredictionTrackerBasePairProb::
drawSvgSquare(const float x, const float y, const size_t size, const float probability, const char* className, const char* tooltip)
{
	std::ostringstream svg;
  svg << "<rect x='" << size_t((x+1.5) * size) << "' y='" << size_t((y+1.5) * size) << "' width='" << size << "' height='" << size << "' class='" << className << "' style='fill:black;fill-opacity:" << sqrt(probability) << ";'>";
	svg << "<title>" << tooltip << "</title>";
	svg << "</rect>\n";
	return svg.str().c_str();
}

////////////////////////////////////////////////////////////////////////////

} // namespace
