#include "SiPMCalib/SiPMCalc/interface/SiPMDarkPdf.hpp"
#include "TMath.h"

SiPMDarkPdf::SiPMDarkPdf( const char* name, const char* title,
                          RooAbsReal& _x,
                          RooAbsReal& _ped,
                          RooAbsReal& _gain,
                          RooAbsReal& _s0,
                          RooAbsReal& _s1,
                          RooAbsReal& _dcfrac,
                          RooAbsReal& _epsilon
                          ) :
  RooAbsPdf( name, title ),
  x( "x", "obs", this, _x ),
  ped( "ped", "pedestal", this, _ped ),
  gain( "gain", "gain", this, _gain ),
  s0( "s0", "comnoise", this, _s0 ),
  s1( "s1", "pixnoise", this, _s1 ),
  dcfrac( "dcfrac1", "dcfrac1", this, _dcfrac ),
  epsilon( "epsilon", "epsilon", this, _epsilon ),
  mdistro( ped, ped+gain, epsilon, sqrt( s0*s0 + s1*s1 ) )
{
}

SiPMDarkPdf::SiPMDarkPdf( const SiPMDarkPdf& other, const char* name ) :
  RooAbsPdf( other, name ),
  x(       "x",        this, other.x    ),
  ped(     "ped",      this, other.ped  ),
  gain(    "gain",     this, other.gain ),
  s0(      "s0",       this, other.s0   ),
  s1(      "s1",       this, other.s1   ),
  dcfrac(  "dcfrac1",  this, other.dcfrac ),
  epsilon( "epsilon",  this, other.epsilon ),
  mdistro( ped, gain, epsilon, sqrt( s0*s0 + s1*s1 ) )
{
}

TObject*
SiPMDarkPdf::clone( const char* name ) const
{
  return new SiPMDarkPdf( *this, name );
}

SiPMDarkPdf::~SiPMDarkPdf(){}

double
SiPMDarkPdf::evaluate() const
{
  mdistro.SetParam( ped, ped+gain, epsilon, sqrt( s0*s0+s1*s1 ) );
  return ( 1-dcfrac ) * TMath::Gaus( x, ped, s0, kTRUE )
         + dcfrac * mdistro.Evaluate( x );
}
