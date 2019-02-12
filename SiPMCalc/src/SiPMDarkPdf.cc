#include "SiPMCalib/SiPMCalc/interface/SiPMDarkPdf.hpp"

#include "RooRealVar.h"

SiPMDarkPdf::SiPMDarkPdf( const char* name, const char* title,
                          RooRealVar& _x,
                          RooRealVar& _ped,
                          RooRealVar& _gain,
                          RooRealVar& _s0,
                          RooRealVar& _s1,
                          RooRealVar& _dcfrac
                          ) :
  RooAbsPdf( name, title ),
  x( "x", "obs", this, _x ),
  ped( "ped", "pedestal", this, _ped ),
  gain( "gain", "gain", this, _gain ),
  s0( "s0", "comnoise", this, _s0 ),
  s1( "s1", "pixnoise", this, _s1 ),
  dcfrac( "dcfrac1", "dcfrac1", this, _dcfrac ),
  func( ped, gain, s0, s1, dcfrac )
{
}

SiPMDarkPdf::SiPMDarkPdf( const SiPMDarkPdf& other, const char* name ) :
  RooAbsPdf( other, name ),
  x(       "x",        this, other.x    ),
  ped(     "ped",      this, other.ped  ),
  gain(    "gain",     this, other.gain ),
  s0(      "s0",       this, other.s0   ),
  s1(      "s1",       this, other.s1   ),
  dcfrac( "dcfrac1",   this, other.dcfrac ),
  func( ped, gain, s0, s1, dcfrac )
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
  func.SetParam( ped, gain, s0, s1, dcfrac );
  return func.Evaluate( x );
}
