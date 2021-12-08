#include "SiPMCalib/SiPMCalc/interface/CrossTalkPdf.hpp"

#include "RooRealVar.h"
#include "TMath.h"

CrossTalkPdf::CrossTalkPdf( const char* name,
                            const char* title,
                            RooRealVar& _x,
                            RooRealVar& _x0,
                            RooRealVar& _gain,
                            RooRealVar& _s0,
                            RooRealVar& _s1,
                            RooRealVar& _prob ) :
  RooAbsPdf( name, title ),
  x        (    "x", "obs", this, _x ),
  x0       (   "x0", "x0", this, _x0 ),
  gain     ( "gain", "gain", this, _gain ),
  s0       (   "s0", "comnoise", this, _s0 ),
  s1       (   "s1", "pixnoise", this, _s1 ),
  prob     ( "prob", "prob", this, _prob )
{}


CrossTalkPdf::CrossTalkPdf( const CrossTalkPdf& other, const char*name ) :
  RooAbsPdf( other, name ),
  x        (       "x",       this, other.x    ),
  x0       (      "x0",      this, other.x0  ),
  gain     (    "gain",    this, other.gain ),
  s0       (      "s0",      this, other.s0   ),
  s1       (      "s1",      this, other.s1   ),
  prob     (    "prob",    this, other.prob )
{}


CrossTalkPdf::~CrossTalkPdf(){}

TObject*
CrossTalkPdf::clone( const char*name ) const
{
  return new CrossTalkPdf( *this, name );
}


double
CrossTalkPdf::evaluate() const
{
  double ans = 0;

  for( int i = 0; i <= 4; ++i ){
    ans += cross_prob( i ) * gauss_k( i );
  }

  return ans;
}


double
CrossTalkPdf::cross_prob( const int k ) const
{
  if( k == 0 ){
    return 1 / ( 1-prob );
  } else {
    return pow( prob, k );
  }
}


double
CrossTalkPdf::gauss_k( const int k ) const
{
  return TMath::Gaus( x,
                      x0+k * gain,
                      sqrt( s0 * s0+( k+1 ) * s1 * s1 ),
                      kTRUE );
}
