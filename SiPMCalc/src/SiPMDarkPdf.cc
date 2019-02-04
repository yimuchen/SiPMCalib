#include "SiPMCalib/SiPMCalc/interface/SiPMDarkPdf.hpp"

#include "RooRealVar.h"

SiPMDarkPdf::SiPMDarkPdf( const char* name, const char* title,
                          RooRealVar& _x,
                          RooRealVar& _ped,
                          RooRealVar& _gain,
                          RooRealVar& _s0,
                          RooRealVar& _s1,
                          RooRealVar& _dcfrac1,
                          RooRealVar& _dcfrac2,
                          RooRealVar& _acshift
                          ) :
  RooAbsPdf( name, title ),
  x( "x", "obs", this, _x ),
  ped( "ped", "pedestal", this, _ped ),
  gain( "gain", "gain", this, _gain ),
  s0( "s0", "comnoise", this, _s0 ),
  s1( "s1", "pixnoise", this, _s1 ),
  dcfrac1( "dcfrac1", "dcfrac1", this, _dcfrac1 ),
  dcfrac2( "dcfrac2", "dcfrac2", this, _dcfrac2 ),
  acshift( "acshift", "acshift", this, _acshift ),
  func( ped, gain, s0, s1, dcfrac1, dcfrac2, acshift )
{
}

SiPMDarkPdf::SiPMDarkPdf( const SiPMDarkPdf& other, const char* name ) :
  RooAbsPdf( other, name ),
  x(       "x",        this, other.x    ),
  ped(     "ped",      this, other.ped  ),
  gain(    "gain",     this, other.gain ),
  s0(      "s0",       this, other.s0   ),
  s1(      "s1",       this, other.s1   ),
  dcfrac1( "dcfrac1",   this, other.dcfrac1 ),
  dcfrac2( "dcfrac2",  this, other.dcfrac2 ),
  acshift( "acshift",  this, other.acshift ),
  func( ped, gain, s0, s1, dcfrac1, dcfrac2, acshift )
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
  func.SetParam( ped, gain, s0, s1, dcfrac1, dcfrac2, acshift );
  return func.Evaluate( x );
}
