#ifndef SIPMCALIB_INVSQCALC_INVSQFUNC_HPP
#define SIPMCALIB_INVSQCALC_INVSQFUNC_HPP

#include <cmath>

/**
 * @brief The inverse square law formula with just the z coordinate as the
 * variable.
 *
 * Pure C++ implementation of the inverse square law function:
 *
 * L(z; L0, x0, z0, P ) = L0 (z-z0)/(x0^2 + (z-z0)^2)^3/2 + P
 *
 * As the scale and uncertainty of L0 is typically much larger than the other
 * parameters, which causes issues with plotting and error estimations, and
 * additional parameter of N0 has been provided such that the L0 in the
 * function actually is L0 = L*N0. As this parameter has perfect correlation
 * with N, this parameter should only ever be set as a constant.
 */
inline double
InvSq_Z( const double*vz, const double*param )
{
  const double z  = vz[0];
  const double z0 = param[0];
  const double o  = param[1];
  const double L  = param[2];
  const double P  = param[3];
  const double N0 = param[4];

  const double D2 = ( o * o )+(( z-z0 ) * ( z-z0 ));

  // return N/D2 + P;
  return ( N0 * L * ( z-z0 ) ) / ( D2 * sqrt( D2 ) )+P;
}


#endif
