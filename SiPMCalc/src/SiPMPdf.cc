#include "SiPMCalib/SiPMCalc/interface/SiPMPdf.hpp"

#include "RooConstVar.h"
#include "RooRealVar.h"
#include "Rtypes.h"

#include "TH1D.h"
#include "TMath.h"

#include <algorithm>
#include <iostream>

// ClassImp( SiPMPdf );

// Full model construction
SiPMPdf::SiPMPdf( const char* name, const char* title,
                  RooAbsReal& _x,
                  RooAbsReal& _ped,
                  RooAbsReal& _gain,
                  RooAbsReal& _s0,
                  RooAbsReal& _s1,
                  RooAbsReal& _mean,
                  RooAbsReal& _lambda,
                  RooAbsReal& _alpha,
                  RooAbsReal& _beta,
                  RooAbsReal& _dcfrac,
                  RooAbsReal& _epsilon
                  ) :
  RooAbsPdf( name, title ),
  x(          "x",      "obs",          this, _x ),
  ped(        "ped",    "pedestal",     this, _ped ),
  gain(       "gain",   "gain",         this, _gain ),
  s0(         "s0",     "comnoise",     this, _s0 ),
  s1(         "s1",     "pixnoise",     this, _s1 ),
  mean(       "mean",   "mean",         this, _mean ),
  lambda(     "lambda", "crosstalk",    this, _lambda ),
  alpha(      "alpha",  "alpha",        this, _alpha ),
  beta(       "beta",   "beta",         this, _beta ),
  dcfraction( "dcfrac", "darkfraction", this, _dcfrac ),
  epsilon(    "eps",    "epsilon",      this, _epsilon ),
  mdistro( ped, ped+gain, epsilon, sqrt( s0*s0 + s1*s1 ) )
{
}

// No Dark current model
SiPMPdf::SiPMPdf( const char* name, const char* title,
                  RooAbsReal& _x,
                  RooAbsReal& _ped,
                  RooAbsReal& _gain,
                  RooAbsReal& _s0,
                  RooAbsReal& _s1,
                  RooAbsReal& _mean,
                  RooAbsReal& _lambda,
                  RooAbsReal& _alpha,
                  RooAbsReal& _beta
                  ) :
  RooAbsPdf( name, title ),
  x(          "x",      "obs",          this, _x ),
  ped(        "ped",    "pedestal",     this, _ped ),
  gain(       "gain",   "gain",         this, _gain ),
  s0(         "s0",     "comnoise",     this, _s0 ),
  s1(         "s1",     "pixnoise",     this, _s1 ),
  mean(       "mean",   "mean",         this, _mean ),
  lambda(     "lambda", "crosstalk",    this, _lambda ),
  alpha(      "alpha",  "alpha",        this, _alpha ),
  beta(       "beta",   "beta",         this, _beta ),
  dcfraction( "dcfrac", "darkfraction", this, RooFit::RooConst( 0 ) ),
  epsilon(    "eps",    "epsilon",      this, RooFit::RooConst( 0.01 ) ),
  mdistro( ped, ped+gain, epsilon, sqrt( s0*s0 + s1*s1 ) )
{}

SiPMPdf::SiPMPdf( const char* name, const char* title,
                  RooAbsReal& _x,
                  RooAbsReal& _ped,
                  RooAbsReal& _gain,
                  RooAbsReal& _s0,
                  RooAbsReal& _s1,
                  RooAbsReal& _mean,
                  RooAbsReal& _lambda
                  ) :
  RooAbsPdf( name, title ),
  x(          "x",      "obs",        this, _x ),
  ped(        "ped",    "pedestal",   this, _ped ),
  gain(       "gain",   "gain",       this, _gain ),
  s0(         "s0",     "comnoise",   this, _s0 ),
  s1(         "s1",     "pixnoise",   this, _s1 ),
  mean(       "mean",   "mean",       this, _mean ),
  lambda(     "lambda", "crosstalk",  this, _lambda ),
  alpha(      "alpha",  "alpha",      this, RooFit::RooConst( 0 ) ),
  beta(       "beta",   "beta",       this, RooFit::RooConst( 1000 ) ),
  dcfraction( "dcfrac", "dcfraction", this, RooFit::RooConst( 0 ) ),
  epsilon(    "eps",    "epsilon",    this, RooFit::RooConst( 0.01 ) ),
  mdistro( ped, ped+gain, epsilon, sqrt( s0*s0 + s1*s1 ) )
{}


SiPMPdf::SiPMPdf( const SiPMPdf& other, const char* name ) :
  RooAbsPdf( other, name ),
  x( "x", this, other.x ),
  ped( "ped", this,  other.ped ),
  gain( "gain", this, other.gain ),
  s0( "s0", this, other.s0 ),
  s1( "s1", this, other.s1 ),
  mean( "mean", this, other.mean ),
  lambda( "lambda", this, other.lambda ),
  alpha( "alpha", this, other.alpha ),
  beta( "beta",  this, other.beta ),
  dcfraction( "dcfrac", this, other.dcfraction ),
  epsilon( "eps", this, other.epsilon ),
  mdistro( ped, ped+gain, epsilon, sqrt( s0*s0 + s1*s1 ) )
{}

SiPMPdf::~SiPMPdf(){}

TObject*
SiPMPdf::clone( const char* name ) const
{
  return new SiPMPdf( *this, name );
}

double
SiPMPdf::evaluate() const
{
  double prob = gen_poisson( 0 ) * gauss_k( 0 );

  for( int k = 1; k < mean + 2*TMath::Sqrt( mean ) + 15; ++k ){
    double probk = binomial_prob( k, 0 ) * gauss_k( k );

    for( int i = 1; i <= k; ++i ){
      probk += binomial_prob( k, i ) * ap_eff( k, i );
    }

    prob += gen_poisson( k ) * probk;
  }

  if( prob <= 0 ){
    prob = std::numeric_limits<double>::min();
    // Forcing non-zero to avoid fit crashing.
  }

  return prob;
}

double
SiPMPdf::gen_poisson( const int k ) const
{
  if( lambda == 0 ){
    return TMath::Poisson( k, mean );
  }

  const double y = ( mean + k * lambda );
  double prod    = 1;

  for( int i = 1; i <= k; ++i ){
    prod *= y;
    prod /= (double)( i );
  }

  prod *= mean / y;
  return prod * TMath::Exp( -y );
}

double
SiPMPdf::ap_eff( const int k, const int i ) const
{
  const double pk = ped + gain * k;
  const double sk = TMath::Sqrt( s0*s0 + k*s1*s1 );
  const double y  = x - pk;

  if( i > 1 ){
    if( y < 0 ){ return 0; }

    double ans = TMath::Exp( -y / beta ) / beta;

    for( int j = 1; j <= i-1; ++j ){
      ans *= ( y / beta );
      ans /= double(j);
    }

    return ans;
  } else {
    const double ap    = TMath::Exp( -y / beta ) / beta;
    const double smear = 0.5 * ( 1+TMath::Erf( y/TMath::Sqrt( 2*sk*sk ) ) );
    return ap * smear;
  }
}

double
SiPMPdf::gauss_k( const int k ) const
{
  const double pk = ped + gain * k;
  const double sk = TMath::Sqrt( s0*s0 + k*s1*s1 );

  // Adding dark current to all Geiger discharge peaks
  if( dcfraction == 0. || k > 0  ){
    return TMath::Gaus( x, pk, sk, true );
  } else {
    mdistro.SetParam( ped, ped+gain, epsilon, TMath::Sqrt( s0*s0+s1*s1 ) );
    return ( 1-dcfraction ) * TMath::Gaus( x, pk, sk, true )
           + dcfraction * mdistro.Evaluate( x-pk );
  }
}

double
SiPMPdf::binomial_prob( const int k, const int i ) const
{
  if( alpha > 0 ){
    return TMath::BinomialI( alpha, k, i )
           - TMath::BinomialI( alpha, k, i+1 );
  } else {
    return i == 0 ? 1 : 0;
  }
}
