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
 * Where the variable is in the z direction.
 */
inline double
InvSq_Z( const double*vz, const double*param )
{
  const double z  = vz[0];
  const double z0 = param[0];
  const double o  = param[1];
  const double L  = param[2];
  const double P  = param[3];

  const double D2 = ( o * o )+(( z-z0 ) * ( z-z0 ));

  // return N/D2 + P;
  return ( L * ( z-z0 ) ) / ( D2 * sqrt( D2 ) )+P;
}


#endif
