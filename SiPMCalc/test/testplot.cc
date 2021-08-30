#include "SiPMCalib/SiPMCalc/interface/SiPMPdf.hpp"
#include "UserUtils/Common/interface/Maths.hpp"
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooRealVar.h"

#include "TRandom3.h"
#include <algorithm>
#include <iostream>

double
est_enf( const double mean,
         const double lambda,
         const double alpha,
         const double ins0  )
{
  static const unsigned npoints  = 1e6;
  static const unsigned n_events = 1e6;
  static const double gain       = 100;
  static const double s1         = 1;
  static const double beta       = gain * 0.5 ;
  const double s0                = ins0 * gain;

  std::vector<double> gen_poisson_int;

  // Generating the integral
  gen_poisson_int.push_back( 0 );

  for( unsigned i = 1; i < npoints; ++i ){
    gen_poisson_int.push_back( gen_poisson_int.back()+  SiPMPdf::GeneralPoissonProb( i-1, mean, lambda ) );
  }

  gen_poisson_int.push_back( 1 );

  // Generating the random events
  std::vector<double> values;
  TRandom3 rand( 65536 );

  for( unsigned i = 0; i < n_events; ++i ){
    double erand = rand.Uniform( 1 );

    // Making a random number of pes based on generalized poisson
    unsigned npe = std::upper_bound( gen_poisson_int.begin()
                                   , gen_poisson_int.end()
                                   , erand  ) - gen_poisson_int.begin() -1;

    double sk = TMath::Sqrt( s0*s0 + npe * s1 * s1 );
    double v  = npe * gain + rand.Gaus( 0, sk );

    unsigned nap = rand.Binomial( npe, alpha );

    for( unsigned j = 0; j < nap; ++j ){
      v += rand.Exp( beta ) + rand.Gaus( 0, s0 );
    }

    values.push_back( v );
  }

  const double sdev = usr::StdDev( values );
  const double m    = usr::Mean( values );

  return ( sdev*sdev )/( m*m ) * mean;
}


int
main( int argc, char* argv[] )
{

  for( double mean = 0.5; mean <= 20.0; mean += 0.5 ){
    for( double lambda = 0.0; lambda <= 0.125; lambda += 0.025 ){
      for( double alpha = 0.0; alpha <= 0.04; alpha += 0.01 ){
        for( double s0 : {0.01, 0.2, 0.5, 1.0} ){
          std::cout <<
            usr::fstr( "%8.2lf %8.5lf %8.5lf %8.5lf %8.5lf"
                     ,  mean
                     , lambda
                     , alpha
                     , s0
                     , est_enf( mean, lambda, alpha, s0 ) ) << std::endl;
        }
      }
    }
  }

  return 0;
}
