#include "SiPMCalib/SiPMCalc/interface/SiPMSumGaus_Functor.hpp"

#include "TMath.h"

SiPMSumGaus_Functor::SiPMSumGaus_Functor( const unsigned n ) :
  _ndim( 5+n ){}

SiPMSumGaus_Functor::~SiPMSumGaus_Functor(){};

ROOT::Math::IBaseFunctionMultiDim*
SiPMSumGaus_Functor::Clone() const
{ return new SiPMSumGaus_Functor( _ndim ); }


double
SiPMSumGaus_Functor::gauss_k( const int k, const double* _x ) const
{
  const double pk = ped( _x ) + gain( _x ) * k;
  const double sk = TMath::Sqrt( s0( _x )* s0( _x ) + k * s1( _x )*s1( _x ) );
  return TMath::Gaus( x( _x ), pk, sk );
}

double
SiPMSumGaus_Functor::DoEval( const double* x ) const
{
  double prob = 0;

  for( unsigned i = 0; i < NPeak(); ++i ){
    prob += amp( x, i ) * gauss_k( i, x );
  }

  return prob;
}
