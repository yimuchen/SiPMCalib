#include "SiPMCalib/SiPMCalc/interface/SiPMDarkFunc.hpp"
#include "UserUtils/Common/interface/Maths.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>

#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_fft_real.h>

#include "TMath.h"

static const unsigned reserve_size = 4 * 1024 * 1024;

MDistro::MDistro() :
  loEdge( 0 ),
  hiEdge( 0 ),
  epsilon( 0 ),
  width( 0 ),
  spline( ROOT::Math::Interpolation::kCSPLINE ),
  spline_acc( ROOT::Math::Interpolation::kCSPLINE )
{
  xArray.reserve( reserve_size );
  convArray.reserve( reserve_size );
  convTempArray.reserve( reserve_size );
  mfuncArray.reserve( reserve_size );
  gaussArray.reserve( reserve_size );
}

MDistro::MDistro(
  const double lo,
  const double hi,
  const double ep,
  const double w ) :
  loEdge( 0 ),
  hiEdge( 0 ),
  epsilon( 0 ),
  width( 0 ),
  spline( ROOT::Math::Interpolation::kCSPLINE ),
  spline_acc( ROOT::Math::Interpolation::kCSPLINE )
{
  xArray.reserve( reserve_size );
  convArray.reserve( reserve_size );
  convTempArray.reserve( reserve_size );
  mfuncArray.reserve( reserve_size );
  gaussArray.reserve( reserve_size );
  SetParam( lo, hi, ep, w );
}

MDistro::~MDistro()
{
}

void
MDistro::SetParam(
  const double lo,
  const double hi,
  const double ep,
  const double w )
{
  loEdge  = std::min( lo, hi );
  hiEdge  = std::max( lo, hi );
  epsilon = fabs( ep );
  width   = fabs( w );

  ParamHash();
}


void
MDistro::ParamHash()
{
  const uint64_t hashval = usr::HashValue( loEdge )
                           ^ usr::HashValue( hiEdge )
                           ^ usr::HashValue( epsilon )
                           ^ usr::HashValue( width );

  if( hashval != paramHash ){
    // Saving Hash, recalculating FFT arrays.
    paramHash = hashval;
    MakeFFTArray();
  }
}

void
MDistro::MakeFFTArray()
{
  const unsigned nbins = usr::RoundUpToP2( std::max( {
    std::ceil( ( xMax() - xMin() )/epsilon ) +1,
    2048.
  } ) );
  const double opepsilon = ( xMax() - xMin() )/( (double)( nbins-1 ) );
  const double xcen      = ( xMax() + xMin() ) / 2;

  convTempArray.resize( nbins );
  mfuncArray.resize( nbins );
  gaussArray.resize( nbins );
  xArray.resize( nbins );
  convArray.resize( nbins );
  accArray.resize( nbins );

  for( unsigned i = 0; i < nbins; ++i ){
    const double x = xMin() + i * opepsilon;
    xArray[i]     = x;
    mfuncArray[i] = MFuncEval( x );
    gaussArray[i] = TMath::Gaus( x, xcen, width, kTRUE );
  }


  gsl_fft_real_radix2_transform( mfuncArray.data(), 1, nbins );
  gsl_fft_real_radix2_transform( gaussArray.data(), 1, nbins );

  for( unsigned i = 0; i < nbins; ++i ){
    convTempArray[i] = mfuncArray[i] * gaussArray[i];
  }

  gsl_fft_halfcomplex_radix2_inverse( convTempArray.data(), 1, nbins );


  // Shifting by half period (circular convolution theorem)
  // opepsilon for proper normalization.
  for( unsigned i = 0; i < nbins; ++i ){
    convArray[i] = convTempArray[( i+nbins/2 )%nbins]  * opepsilon;

    if( i > 0 ){
      accArray[i] = accArray[i-1] + convArray[i] * opepsilon ;
    } else {
      accArray[0] = convArray[0] * opepsilon ;
    }
  }

  spline.SetData( xArray, convArray );
  spline_acc.SetData( xArray, accArray );
}

double
MDistro::EdgeDist() const
{
  return fabs( hiEdge - loEdge );
}

static const double _edge_mult  = 5;
static const double _width_mult = 10;


double
MDistro::xMin() const
{
  return loEdge - std::max( _edge_mult*EdgeDist(), _width_mult*width );
}

double
MDistro::xMax() const
{
  return hiEdge + std::max( _edge_mult*EdgeDist(), _width_mult*width );
}

double
MDistro::Evaluate( const double x ) const
{
  if( x < xArray.front() || xArray.back() < x ){
    return 0;
  } else {
    return std::max( spline.Eval( x ), 0. );
  }
}

double
MDistro::EvaluateAccum( const double x ) const
{
  if( x < xArray.front() ){
    return 0;
  } else if( x > xArray.back() ){
    return 1;
  } else {
    return std::max( spline_acc.Eval( x ), 0.0 );
  }
}


double
MDistro::MFuncEval( const double x ) const
{
  if( loEdge + epsilon < x && x < hiEdge - epsilon ){
    return ( 1 / ( x-loEdge ) + 1 / ( hiEdge - x ) )/
           ( 2*std::log( ( EdgeDist() - epsilon ) / epsilon ) );
  } else {
    return 0;
  }
}
