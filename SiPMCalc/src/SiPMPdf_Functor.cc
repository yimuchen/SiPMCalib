#include "SiPMCalib/SiPMCalc/interface/SiPMPdf_Functor.hpp"
#include "TMath.h"

#include <iostream>

SiPMPdf_Functor::SiPMPdf_Functor(){}

SiPMPdf_Functor::~SiPMPdf_Functor(){}

ROOT::Math::IBaseFunctionMultiDim*
SiPMPdf_Functor::Clone() const { return 0; }

double
SiPMPdf_Functor::g_poisson( const int k, const double* _x ) const
{
  const double y = ( mean(_x) + k * lambda(_x) );
  double prod    = 1/y;

  for( int i = 1; i <= k; ++i ){
    prod *= y;
    prod /= (double)( i );
  }

  return mean(_x) * prod * TMath::Exp( -y );
}

double
SiPMPdf_Functor::after_pulse_eff(
  const int k,
  const int i,
  const double* _x
  ) const
{
  const double pk = ped(_x) + gain(_x) * k;
  const double sk = TMath::Sqrt( s0(_x) * s0(_x) + k * s1(_x) * s1(_x) );
  const double y  = x(_x) - pk;
  if( y < 0 ){
    return 0;
  } else if( i > 1 ){
    double prod = 1 / beta(_x);

    for( int j = 1; j <= i-1; ++j ){
      prod *= y / ( j* beta(_x) );
    }

    return prod * TMath::Exp( -y / beta(_x) );

  } else {
    const double num = TMath::Exp( -y / beta(_x) );
    const double den = TMath::Sqrt( 2 * TMath::Pi() ) * sk * beta(_x);
    const double err = ( 1 + TMath::Erf( y / ( TMath::Sqrt( 2 ) * sk ) ) )/ 2.;
    return err * num /  den;
  }
}

double
SiPMPdf_Functor::gauss_k( const int k, const double* _x ) const
{
  const double pk = ped(_x) + gain(_x) * k;
  const double sk = TMath::Sqrt( s0(_x)* s0(_x) + k * s1(_x)*s1(_x) );
  return TMath::Gaus( x(_x), pk, sk );
}

double
SiPMPdf_Functor::DoEval( const double* _x ) const
{
  double prob = g_poisson( 0, _x ) * gauss_k( 0, _x );

  for( int k = 1; k < mean(_x) * 8; ++k ){
    double probk = TMath::BinomialI( alpha(_x), k, 0 )
                    * gauss_k( k, _x );

    for( int i = 1; i <= k; ++i ){
      probk += TMath::BinomialI( alpha(_x), k, i )
                * after_pulse_eff( k, i, _x );
    }

    prob += g_poisson( k, _x ) * probk;
  }

  return prob;
}
