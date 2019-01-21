#include "SiPMCalib/SiPMCalc/interface/SiPMPdf.hpp"

#include "RooConstVar.h"
#include "RooRealVar.h"
#include "Rtypes.h"

#include "TH1D.h"
#include "TMath.h"

#include <algorithm>
#include <iostream>

// ClassImp( SiPMPdf );

SiPMPdf::SiPMPdf( const char* name, const char* title,
                  RooRealVar& _x,
                  RooRealVar& _ped,
                  RooRealVar& _gain,
                  RooRealVar& _s0,
                  RooRealVar& _s1,
                  RooRealVar& _mean,
                  RooRealVar& _lambda,
                  RooRealVar& _acfrac,
                  RooRealVar& _acshift,
                  RooRealVar& _acwidth,
                  RooRealVar& _alpha,
                  RooRealVar& _beta
                  ) :
  RooAbsPdf( name, title ),
  x( "x", "obs", this, _x ),
  ped( "ped", "pedestal", this, _ped ),
  gain( "gain", "gain", this, _gain ),
  s0( "s0", "comnoise", this, _s0 ),
  s1( "s1", "pixnoise", this, _s1 ),
  mean( "mean", "mean", this, _mean ),
  lambda( "lambda", "crosstalk", this, _lambda ),
  acfrac( "acfrac", "ac_fraction", this, _acfrac ),
  acshift( "acshift", "ac_shift", this, _acshift ),
  acwidth( "acwidth", "ac_width", this, _acwidth ),
  alpha( "alpha", "alpha", this, _alpha ),
  beta( "beta", "beta", this, _beta )
{
}

SiPMPdf::SiPMPdf( const char* name, const char* title,
                  RooRealVar& _x,
                  RooRealVar& _ped,
                  RooRealVar& _gain,
                  RooRealVar& _s0,
                  RooRealVar& _s1,
                  RooRealVar& _mean,
                  RooRealVar& _lambda,
                  RooRealVar& _alpha,
                  RooRealVar& _beta
                  ) :
  RooAbsPdf( name, title ),
  x( "x", "obs", this, _x ),
  ped( "ped", "pedestal", this, _ped ),
  gain( "gain", "gain", this, _gain ),
  s0( "s0", "comnoise", this, _s0 ),
  s1( "s1", "pixnoise", this, _s1 ),
  mean( "mean", "mean", this, _mean ),
  lambda( "lambda", "crosstalk", this, _lambda ),
  acfrac( "acfrac", "ac_fraction", this, RooFit::RooConst(0) ),
  acshift( "acshift", "ac_shift", this, RooFit::RooConst(0) ),
  acwidth( "acwidth", "ac_width", this, RooFit::RooConst(1) ),
  alpha( "alpha", "alpha", this, _alpha ),
  beta( "beta", "beta", this, _beta )
{}

SiPMPdf::SiPMPdf( const char* name, const char* title,
                  RooRealVar& _x,
                  RooRealVar& _ped,
                  RooRealVar& _gain,
                  RooRealVar& _s0,
                  RooRealVar& _s1,
                  RooRealVar& _mean,
                  RooRealVar& _lambda,
                  RooRealVar& _acfrac,
                  RooRealVar& _acshift,
                  RooRealVar& _acwidth
                  ) :
  RooAbsPdf( name, title ),
  x( "x", "obs", this, _x ),
  ped( "ped", "pedestal", this, _ped ),
  gain( "gain", "gain", this, _gain ),
  s0( "s0", "comnoise", this, _s0 ),
  s1( "s1", "pixnoise", this, _s1 ),
  mean( "mean", "mean", this, _mean ),
  lambda( "lambda", "crosstalk", this, _lambda ),
  acfrac( "acfrac", "ac_fraction", this, _acfrac ),
  acshift( "acshift", "ac_shift", this, _acshift ),
  acwidth( "acwidth", "ac_width", this, _acwidth ),
  alpha( "alpha", "alpha", this, RooFit::RooConst( 0 ) ),
  beta( "beta", "beta", this, RooFit::RooConst( _x.getMax()*1e3 ) )
{}

SiPMPdf::SiPMPdf( const char* name, const char* title,
                  RooRealVar& _x,
                  RooRealVar& _ped,
                  RooRealVar& _gain,
                  RooRealVar& _s0,
                  RooRealVar& _s1,
                  RooRealVar& _mean,
                  RooRealVar& _lambda
                  ) :
  RooAbsPdf( name, title ),
  x( "x", "obs", this, _x ),
  ped( "ped", "pedestal", this, _ped ),
  gain( "gain", "gain", this, _gain ),
  s0( "s0", "comnoise", this, _s0 ),
  s1( "s1", "pixnoise", this, _s1 ),
  mean( "mean", "mean", this, _mean ),
  lambda( "lambda", "crosstalk", this, _lambda ),
  acfrac( "acfrac", "ac_fraction", this, RooFit::RooConst( 0 ) ),
  acshift( "acshift", "ac_shift", this, RooFit::RooConst( 0 ) ),
  acwidth( "acwidth", "ac_width", this, RooFit::RooConst( 1 ) ),
  alpha( "alpha", "alpha", this, RooFit::RooConst( 0 ) ),
  beta( "beta", "beta", this, RooFit::RooConst( _x.getMax()*1e3 ) )
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
  acfrac( "acfrac", this, other.acfrac ),
  acshift( "acshift", this, other.acshift ),
  acwidth( "acwith", this, other.acwidth ),
  alpha( "alpha", this, other.alpha ),
  beta( "beta",  this, other.beta )
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

  for( int k = 1; k < mean + TMath::Sqrt( mean ) + 8; ++k ){
    double probk = binomial_prob( k, 0 ) * gauss_k( k );

    for( int i = 1; i <= k; ++i ){
      probk += binomial_prob( k, i ) * ap_eff( k, i );
    }

    prob += gen_poisson( k ) * probk;
  }

  return prob;
}


double
SiPMPdf::gen_poisson( const int k ) const
{
  const double y = ( mean + k * lambda );
  double prod    = 1/y;

  for( int i = 1; i <= k; ++i ){
    prod *= y;
    prod /= (double)( i );
  }

  return mean * prod * TMath::Exp( -y );
}

double
SiPMPdf::ap_eff( const int k, const int i ) const
{
  const double pk = ped + gain * k;
  const double sk = TMath::Sqrt( s0*s0 + k*s1*s1 );
  const double y  = x - pk;

  if( y < 0 ){
    return 0;
  } else if( i > 1 ){
    double prod = 1 / beta;

    for( int j = 1; j <= i-1; ++j ){
      prod *= y / ( j* beta );
    }

    return prod * TMath::Exp( -y / beta );

  } else {
    const double num = TMath::Exp( -y / beta );
    const double den = TMath::Sqrt( 2 * TMath::Pi() ) * sk * beta;
    const double err = ( 1 + TMath::Erf( y / ( TMath::Sqrt( 2 ) * sk ) ) )/ 2.;
    return err * num /  den;
  }
}

double
SiPMPdf::gauss_k( const int k ) const
{
  const double pk = ped + gain * k;
  const double sk = TMath::Sqrt( s0*s0 + k*s1*s1 );

  if( k == 0 ){
    const double acpk = ped + gain * k - acshift;
    return ( 1- acfrac ) * TMath::Gaus( x, pk, sk )
           + acfrac * TMath::Gaus( x, acpk, acwidth );
  } else {
    return TMath::Gaus( x, pk, sk );
  }
}

double
SiPMPdf::binomial_prob( const int k, const int i ) const
{
  return TMath::BinomialI( alpha, k, i )
         - TMath::BinomialI( alpha, k, i+1 );
}
