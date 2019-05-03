
#include "IntaRNA/OutputHandlerCsv.h"

#if INTARNA_MULITHREADING
	#include <omp.h>
#endif

#include <boost/algorithm/string.hpp>

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////

std::map<OutputHandlerCsv::ColType,std::string> OutputHandlerCsv::colType2string;

////////////////////////////////////////////////////////////////////////

OutputHandlerCsv::OutputHandlerCsv(
		  std::ostream & out
		, const InteractionEnergy & energy
		, const ColTypeList colOrder
		, const std::string& colSep
		, const bool printHeader
		, const std::string& listSep
		)
 :	out(out)
	, energy(energy)
	, colOrder(colOrder)
	, colSep(colSep)
	, listSep(listSep)
{
	// init mapping of coltypes to string
	initColType2string();

	// print CSV header of column names
	if (printHeader) {
		// ensure outputs do not intervene
#if INTARNA_MULITHREADING
		#pragma omp critical(intarna_omp_outputStreamUpdate)
#endif
		{
			out <<getHeader(colOrder,colSep);
		}
	}
}

////////////////////////////////////////////////////////////////////////

OutputHandlerCsv::~OutputHandlerCsv()
{
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_outputStreamUpdate)
#endif
	{
		// force output
		out.flush();
	}
}

////////////////////////////////////////////////////////////////////////

void
OutputHandlerCsv::
add( const Interaction & i, const OutputConstraint & outConstraint )
{
#if INTARNA_IN_DEBUG_MODE
	// debug checks
	if ( i.basePairs.size() > 0 && ! i.isValid() ) {
		throw std::runtime_error("OutputHandlerCsv::add() : given interaction is not valid : "+toString(i));
	}
#endif

	// special handling if no base pairs present
	if (i.basePairs.size() == 0) {
		return;
	}

	// get interaction start/end per sequence
	const size_t i1 = i.basePairs.begin()->first;
	const size_t j1 = i.basePairs.rbegin()->first;
	const size_t i2 = i.basePairs.begin()->second;
	const size_t j2 = i.basePairs.rbegin()->second;

	// get individual energy contributions
	InteractionEnergy::EnergyContributions contr = energy.getE_contributions(i);

	// ensure outputs do not intervene
	{
		std::stringstream outTmp;

		for (auto col = colOrder.begin(); col != colOrder.end(); col++) {
			// print separator if needed
			if (col != colOrder.begin()) {
				outTmp <<colSep;
			}
			// print this column information
			switch ( *col ) {

			case id1:
				// ensure no colSeps are contained
				outTmp <<boost::replace_all_copy(energy.getAccessibility1().getSequence().getId(), colSep, "_");
				break;

			case id2:
				// ensure no colSeps are contained
				outTmp <<boost::replace_all_copy(energy.getAccessibility2().getSequence().getId(), colSep, "_");
				break;

			case seq1:
				outTmp <<energy.getAccessibility1().getSequence().asString();
				break;

			case seq2:
				outTmp <<energy.getAccessibility2().getAccessibilityOrigin().getSequence().asString();
				break;

			case subseq1:
				outTmp <<energy.getAccessibility1().getSequence().asString().substr(i1, j1-i1+1);
				break;

			case subseq2:
				outTmp <<energy.getAccessibility2().getAccessibilityOrigin().getSequence().asString().substr(j2, i2-j2+1);
				break;

			case subseqDP:
				outTmp <<energy.getAccessibility1().getSequence().asString().substr(i1, j1-i1+1)
					<<'&'
					<<energy.getAccessibility2().getAccessibilityOrigin().getSequence().asString().substr(j2, i2-j2+1);
				break;

			case subseqDB:
				outTmp <<(i1+1)
					<<energy.getAccessibility1().getSequence().asString().substr(i1, j1-i1+1)
					<<'&'
					<<(j2+1)
					<<energy.getAccessibility2().getAccessibilityOrigin().getSequence().asString().substr(j2, i2-j2+1);
				break;

			case start1:
				outTmp <<(i1+1);
				break;

			case end1:
				outTmp <<(j1+1);
				break;

			case start2:
				outTmp <<(j2+1);
				break;

			case end2:
				outTmp <<(i2+1);
				break;

			case hybridDP:
				outTmp <<Interaction::dotBracket( i );
				break;

			case hybridDB:
				outTmp <<Interaction::dotBar( i );
				break;

			case hybridDPfull:
				outTmp <<Interaction::dotBracket( i, '(', ')', true );
				break;

			case hybridDBfull:
				outTmp <<Interaction::dotBar( i, true );
				break;

			case E:
				outTmp <<E_2_Ekcal(i.energy);
				break;

			case ED1:
				outTmp <<E_2_Ekcal(contr.ED1);
				break;

			case ED2:
				outTmp <<E_2_Ekcal(contr.ED2);
				break;

			case Pu1:
				outTmp <<E_2_Ekcal(Z_2_E(Z_exp( - E_2_Z(contr.ED1) / energy.getRT() )));
				break;

			case Pu2:
				outTmp <<E_2_Ekcal(Z_2_E(Z_exp( - E_2_Z(contr.ED2 / energy.getRT() ))));
				break;

			case E_init:
				outTmp <<E_2_Ekcal(contr.init);
				break;

			case E_loops:
				outTmp <<E_2_Ekcal(contr.loops);
				break;

			case E_dangleL:
				outTmp <<E_2_Ekcal(contr.dangleLeft);
				break;

			case E_dangleR:
				outTmp <<E_2_Ekcal(contr.dangleRight);
				break;

			case E_endL:
				outTmp <<E_2_Ekcal(contr.endLeft);
				break;

			case E_endR:
				outTmp <<E_2_Ekcal(contr.endRight);
				break;

			case E_hybrid:
				outTmp <<E_2_Ekcal(i.energy - contr.ED1 - contr.ED2);
				break;

			case E_norm:
				outTmp <<E_2_Ekcal(i.energy) / std::log( energy.size1() * energy.size2() );
				break;

			case E_hybridNorm:
				outTmp <<E_2_Ekcal(i.energy - contr.ED1 - contr.ED2) / std::log( energy.size1() * energy.size2() );
				break;

			case E_add:
				outTmp <<E_2_Ekcal(contr.energyAdd);
				break;

			case seedStart1:
				if (i.seed == NULL) {
					outTmp <<std::numeric_limits<E_kcal_type>::signaling_NaN();
				} else {
					if (!outConstraint.bestSeedOnly) {
						// generate list
						outTmp << std::accumulate( std::next( i.seed->begin() )
										, i.seed->end()
										, toString( i.seed->begin()->bp_i.first +1 ) // start with first element
										, [&](std::string a, Interaction::Seed s) {
											 return std::move(a) + listSep + toString( s.bp_i.first +1 ); // extend list
										 });
					} else {
						outTmp <<i.seed->begin()->bp_i.first +1;
					}
				}
				break;

			case seedEnd1:
				if (i.seed == NULL) {
					outTmp <<std::numeric_limits<E_kcal_type>::signaling_NaN();
				} else {
					if (!outConstraint.bestSeedOnly) {
						// generate list
						outTmp << std::accumulate( std::next( i.seed->begin() )
										, i.seed->end()
										, toString( i.seed->begin()->bp_j.first +1 ) // start with first element
										, [&](std::string a, Interaction::Seed s) {
											 return std::move(a) + listSep + toString( s.bp_j.first +1 ); // extend list
										 });
					} else {
						outTmp <<i.seed->begin()->bp_j.first +1;
					}
				}
				break;

			case seedStart2:
				if (i.seed == NULL) {
					outTmp <<std::numeric_limits<E_kcal_type>::signaling_NaN();
				} else {
					if (!outConstraint.bestSeedOnly) {
						// generate list
						outTmp << std::accumulate( std::next( i.seed->begin() )
										, i.seed->end()
										, toString( i.seed->begin()->bp_j.second +1 ) // start with first element
										, [&](std::string a, Interaction::Seed s) {
											 return std::move(a) + listSep + toString( s.bp_j.second +1 ); // extend list
										 });
					} else {
						outTmp <<i.seed->begin()->bp_j.second +1;
					}
				}
				break;

			case seedEnd2:
				if (i.seed == NULL) {
					outTmp <<std::numeric_limits<E_kcal_type>::signaling_NaN();
				} else {
					if (!outConstraint.bestSeedOnly) {
						// generate list
						outTmp << std::accumulate( std::next( i.seed->begin() )
										, i.seed->end()
										, toString( i.seed->begin()->bp_i.second +1 ) // start with first element
										, [&](std::string a, Interaction::Seed s) {
											 return std::move(a) + listSep + toString( s.bp_i.second +1 ); // extend list
										 });
					} else {
						outTmp <<i.seed->begin()->bp_i.second +1;
					}
				}
				break;

			case seedE:
				if (i.seed == NULL) {
					outTmp <<std::numeric_limits<E_kcal_type>::signaling_NaN();
				} else {
					outTmp << E_2_Ekcal(i.seed->begin()->energy);
					if (!outConstraint.bestSeedOnly) {
						// generate list
						// print via std::for_each instead of std::accumulate due to rounding issues of boost::lexical_cast or std::to_string
						// since sometimes (float(int)/100.0) gives strings with 10E-5 deviations of expected value
						std::for_each( ++(i.seed->begin()), i.seed->end(), [&]( const Interaction::Seed & s) {
									 outTmp << listSep << E_2_Ekcal(s.energy);
									});
					}
				}
				break;

			case seedED1:
				if (i.seed == NULL) {
					outTmp <<std::numeric_limits<E_kcal_type>::signaling_NaN();
				} else {
					outTmp << E_2_Ekcal(energy.getED1( i.seed->begin()->bp_i.first, i.seed->begin()->bp_j.first ));
					if (!outConstraint.bestSeedOnly) {
						// generate list
						std::for_each( ++(i.seed->begin()), i.seed->end(), [&]( const Interaction::Seed & s) {
									 outTmp << listSep << E_2_Ekcal(energy.getED1( s.bp_i.first, s.bp_j.first ));
									});
					}
				}
				break;

			case seedED2:
				if (i.seed == NULL) {
					outTmp <<std::numeric_limits<E_kcal_type>::signaling_NaN();
				} else {
					outTmp << E_2_Ekcal(energy.getAccessibility2().getAccessibilityOrigin().getED( i.seed->begin()->bp_j.second, i.seed->begin()->bp_i.second ));
					if (!outConstraint.bestSeedOnly) {
						// generate list
						std::for_each( ++(i.seed->begin()), i.seed->end(), [&]( const Interaction::Seed & s) {
							outTmp << listSep << E_2_Ekcal(energy.getAccessibility2().getAccessibilityOrigin().getED( s.bp_j.second, s.bp_i.second ));
						});
					}
				}
				break;

			case seedPu1:
				if (i.seed == NULL) {
					outTmp <<std::numeric_limits<E_kcal_type>::signaling_NaN();
				} else {
					if (!outConstraint.bestSeedOnly) {
						// generate list
						outTmp << std::accumulate( std::next( i.seed->begin() )
										, i.seed->end()
										, toString( energy.getBoltzmannWeight( energy.getED1( i.seed->begin()->bp_i.first, i.seed->begin()->bp_j.first )) ) // start with first element
										, [&](std::string a, Interaction::Seed s) {
											 return std::move(a) + listSep + toString( energy.getBoltzmannWeight( energy.getED1( s.bp_i.first, s.bp_j.first ) ) ); // extend list
										 });
					} else {
						outTmp <<energy.getBoltzmannWeight( energy.getED1( i.seed->begin()->bp_i.first, i.seed->begin()->bp_j.first ));
					}
				}
				break;

			case seedPu2:
				if (i.seed == NULL) {
					outTmp <<std::numeric_limits<E_kcal_type>::signaling_NaN();
				} else {
					if (!outConstraint.bestSeedOnly) {
						// generate list
						outTmp << std::accumulate( std::next( i.seed->begin() )
										, i.seed->end()
										, toString( energy.getBoltzmannWeight( energy.getAccessibility2().getAccessibilityOrigin().getED( i.seed->begin()->bp_j.second, i.seed->begin()->bp_i.second )) ) // start with first element
										, [&](std::string a, Interaction::Seed s) {
											 return std::move(a) + listSep + toString( energy.getBoltzmannWeight( energy.getAccessibility2().getAccessibilityOrigin().getED( s.bp_j.second, s.bp_i.second ) ) ); // extend list
										 });
					} else {
						outTmp <<energy.getBoltzmannWeight( energy.getAccessibility2().getAccessibilityOrigin().getED( i.seed->begin()->bp_j.second, i.seed->begin()->bp_i.second ));
					}
				}
				break;

			default : throw std::runtime_error("OutputHandlerCsv::add() : unhandled ColType '"+colType2string[*col]+"'");
			}
		}
		outTmp <<'\n';
	#if INTARNA_MULITHREADING
		#pragma omp critical(intarna_omp_outputStreamUpdate)
	#endif
		{
			out << outTmp.str();
		} // omp critical(intarna_omp_outputStreamUpdate)
	}

}

////////////////////////////////////////////////////////////////////////

} // namespace
