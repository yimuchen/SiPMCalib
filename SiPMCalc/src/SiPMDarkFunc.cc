#include "SiPMCalib/SiPMCalc/interface/SiPMDarkFunc.hpp"
#include "UserUtils/Common/interface/Maths.hpp"

#include <cmath>
#include <complex>
#include <algorithm>

#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

#include "TMath.h"

MDistro::MDistro() :
  loEdge(0), hiEdge(0), epsilon(0), width(0),
  spline(nullptr)
{
}

MDistro::MDistro(
    const double lo,
    const double hi,
    const double ep,
    const double w) :
    loEdge(0), hiEdge(0), epsilon(0), width(0),
    spline(nullptr)
{
  SetParam(lo, hi, ep, w);
}

MDistro::~MDistro()
{
}


void MDistro::SetParam(
    const double lo,
    const double hi,
    const double ep,
    const double w)
{
  loEdge = std::min(lo, hi);
  hiEdge = std::max(lo, hi);
  epsilon = fabs(ep);
  width = fabs(w);


  ParamHash();
}


void MDistro::ParamHash()
{
  const uint64_t hashval = usr::HashValue(loEdge) ^ usr::HashValue(hiEdge)
                         ^ usr::HashValue(epsilon) ^ usr::HashValue(width);

  if (hashval != paramHash) {
    // Saving Hash, recalculating FFT arrays.
    paramHash = hashval;
    MakeFFTArray();
  }
}

void MDistro::MakeFFTArray()
{
  const unsigned nbins = usr::RoundUpToP2( std::max({
    std::ceil( (xMax() - xMin())/epsilon ) +1,
    1024.
  }));
  const double opepsilon = (xMax() - xMin())/((double)(nbins-1));
  const double xcen      = (xMax() + xMin() ) / 2 ;
  std::vector<double> convtemp( nbins );
  std::vector<double> mfunc   ( nbins );
  std::vector<double> gauss   ( nbins );
  xArray   .resize( nbins );
  convArray.resize( nbins );

  for( unsigned i = 0 ; i < nbins; ++i ){
    const double x = xMin() + i * opepsilon;
    xArray[i] = x;
    mfunc[i] = MFuncEval(x);
    gauss[i] = TMath::Gaus( x, xcen, width, kTRUE );
  }

  gsl_fft_real_radix2_transform(mfunc.data(),1,nbins);
  gsl_fft_real_radix2_transform(gauss.data(),1,nbins);

  for( unsigned i = 0 ; i < nbins; ++i ){
    convtemp[i] = mfunc[i] * gauss[i];
  }

  gsl_fft_halfcomplex_radix2_inverse(convtemp.data(),1,nbins);

  // Shifting by half period (circular convolution theorem)
  for( unsigned i = 0 ; i < nbins; ++i ){
    convArray[i] = convtemp[(i+nbins/2)%nbins];
  }

  spline = std::make_unique<ROOT::Math::Interpolator>( xArray, convArray );
}

double MDistro::EdgeDist() const
{
  return fabs(hiEdge - loEdge);
}

static const double _edge_mult = 20 ;
static const double _width_mult = 20;


double MDistro::xMin() const
{
  return loEdge - std::max( _edge_mult*EdgeDist(), _width_mult*width) ;
}

double MDistro::xMax() const
{
  return hiEdge + std::min( _edge_mult*EdgeDist(), _width_mult*width);
}

#include <iostream>
double MDistro::Evaluate( const double x ) const
{
  //double ans ;
  //int    spline_err ;
  if( xArray.front() < x && x < xArray.back() ){
    return  std::max( spline->Eval( x ), 0. );
  } else {
    return 0;
  }
}

double MDistro::MFuncEval( const double x ) const
{
  if( loEdge + epsilon < x && x < hiEdge - epsilon ){
    return  (1 / (x-loEdge) + 1 / (hiEdge - x))/
            (2*std::log( (EdgeDist() - epsilon) / epsilon) );
  } else {
    return 0;
  }
}


/******************************************************************************/

SiPMDarkFunc::SiPMDarkFunc(
    const double _ped,
    const double _gain,
    const double _s0,
    const double _s1,
    const double _acfrac1,
    const double _acfrac2,
    const double _acshift
)
{
  SetParam(_ped,_gain,_s0,_s1,_acfrac1,_acfrac2,_acshift);
}

SiPMDarkFunc::~SiPMDarkFunc(){}

void SiPMDarkFunc::SetParam(
    const double _ped,
    const double _gain,
    const double _s0,
    const double _s1,
    const double _acfrac1,
    const double _acfrac2,
    const double _acshift
)
{
  ped     = _ped;
  gain    = _gain;
  s0      = _s0;
  s1      = _s1;
  acfrac1 = _acfrac1;
  acfrac2 = _acfrac2;
  acshift = _acshift;
  _m_primary.SetParam(  low_edge(),high_edge(),epsilon(),w1());
  _m_secondary.SetParam(low_edge(),      ped  ,epsilon(),w1());
}

double SiPMDarkFunc::Evaluate( const double x ) const
{
  if( acshift == 0 ){
    return (1-acfrac1) * TMath::Gaus(x,ped,s0, true )
           + acfrac1 * _m_primary.Evaluate(x);
  } else {
    return (1-acfrac1-acfrac2) * TMath::Gaus( x, ped, s0, true )
           + acfrac1 * _m_primary.Evaluate(x)
           + acfrac2 * _m_secondary.Evaluate(x);
  }
}

double SiPMDarkFunc::EvalM1( const double x ) const
{
  return _m_primary.Evaluate(x);
}

double SiPMDarkFunc::EvalM2( const double x ) const
{
  return _m_secondary.Evaluate(x);
}

double SiPMDarkFunc::low_edge() const { return ped - acshift; }

double SiPMDarkFunc::high_edge() const { return ped + gain; }

double SiPMDarkFunc::epsilon() const {
  return std::max( acshift/1024., 0.1) ;
}

double SiPMDarkFunc::w1() const { return sqrt(s0*s0+s1*s1);}