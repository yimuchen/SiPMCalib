#include "TMath.h"
#include "UserUtils/MathUtils/interface/Measurement/Measurement.hpp"

double
LinearModel( const double* xx, const double* par )
{
  const double x      = xx[0];
  const double gain   = par[0];
  const double offset = par[1];
  const double N      = par[2];
  return gain * x * N+offset;
}


/**
 * @brief Leading order model for nonlinearity
 * @ingroup SiPMCalc
 *
 * - x is the estimated number of effective incoming photons (PDE multiplied).
 *
 * - par[0] is the gain.
 * - par[1] is the total number of pixels in the device.
 * - par[2] is the readout pedestal.
 *
 * The return is the expected readout as estimated by the leading order
 * occupancy model.
 */
double
LOModel( const double* xx, const double* par )
{
  const double x      = xx[0];
  const double gain   = par[0];
  const double N      = par[1];
  const double offset = par[2];
  return gain * N * ( 1-TMath::Exp( -x / N ) )+offset;
}


/**
 * @brief Next-to-Leading order model for nonlinearity
 * @ingroup SiPMCalc
 *
 * - x is the estimated number of effective incoming photons (PDE multiplied).
 *
 * - par[0] is the gain.
 * - par[1] is the total number of pixels.
 * - par[2] is the readout pedestal.
 * - par[3] Probability of secondary pulsing
 * - par[4] Recovery time factor
 *
 * The return is the expected readout as estimated by the leading order
 * occupancy model.
 */
double
NLOModel( const double* xx, const double* par )
{
  const double x     = xx[0];
  const double gain  = par[0];
  const double N     = par[1];
  const double ped   = par[2];
  const double alpha = par[3];
  const double beta  = par[4];

  // Calculating number of fired pixel by LO approximation
  const double LOpar[3]  = {1.0, N, 0.0};
  const double Nfired_LO = LOModel( xx, LOpar );

  // Simple linear correction
  const double NLOLinear  = (( 1-alpha ) * Nfired_LO )+( alpha * x );
  const double NonLinCorr = ( beta+1 ) / ( beta+( x / Nfired_LO ));

  return ( gain * NLOLinear * NonLinCorr )+ped;
}


/**
 * @brief Model of the luminosity power as a function of the bias voltage
 *
 * The model is done using the a exponential model with a vertical offset
 * (pedestal):
 *
 * F(x) = exp( a*x ) + ped
 *
 * -x is the input bias voltage.
 * - par[0] is the natural log of the exponent.
 * - par[1] is the vertical offset.
 */
double
BiasModel( const double* xx, const double* par  )
{
  const double x   = xx[0];
  const double ex  = par[0];
  const double ped = par[1];

  return TMath::Exp( ex * x )+ped;
}


// From collected data.
usr::Measurement
A( const unsigned idx )
{
  static const double Araw_lo[7] = {
    0.7081, 1668.4, 2508.0, 2977.1, 4195.3, 5273.5, 6578.1};

  static const double Araw_hi[7] = {
    0.7015, 1670.6, 2511.8, 2983.0, 4200.5, 5286.7, 6604.4};
  if( idx < 1 || idx > 6 ){ return usr::Measurement( 0, 0, 0 ); }
  if( idx == 6 ){ return usr::Measurement( 1, 0, 0 ); }
  const double lo   = ( Araw_lo[idx]-Araw_hi[0] ) / ( Araw_hi[6]-Araw_lo[0] );
  const double hi   = ( Araw_hi[idx]-Araw_lo[0] ) / ( Araw_lo[6]-Araw_hi[0] );
  const double diff = fabs( hi-lo );

  return usr::Measurement( ( lo+hi ) / 2, diff, diff );
}


usr::Measurement
B( const unsigned idx )
{
  static const double Braw_lo[7] = {
    0.696, 31069, 8633.1, 2480.4, 701.10, 97.002, 20.379};

  static const double Braw_hi[7] = {
    0.700, 31328, 8691.0, 2497.4, 704.95, 97.330, 20.470};
  if( idx < 1 || idx > 6 ){ return usr::Measurement( 0, 0, 0 ); }
  if( idx == 1 ){ return usr::Measurement( 1, 0, 0 ); }
  const double lo   = ( Braw_lo[idx]-Braw_hi[0] ) / ( Braw_hi[1]-Braw_lo[0] );
  const double hi   = ( Braw_hi[idx]-Braw_lo[0] ) / ( Braw_lo[1]-Braw_hi[0] );
  const double diff = fabs( hi-lo ) / 2;

  return usr::Measurement( ( lo+hi ) / 2, diff, diff );
}
