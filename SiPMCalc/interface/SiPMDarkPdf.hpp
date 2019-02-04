#ifndef SIPMCALIB_SIPMCALC_SIPMDARKPDF_HPP
#define SIPMCALIB_SIPMCALC_SIPMDARKPDF_HPP

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooTrace.h"

#include "SiPMCalib/SiPMCalc/interface/SiPMDarkFunc.hpp"

class SiPMDarkPdf : public RooAbsPdf
{
public:
  SiPMDarkPdf( const char*, const char*,
               RooRealVar& x,
               RooRealVar& ped,
               RooRealVar& gain,
               RooRealVar& s0,
               RooRealVar& s1,
               RooRealVar& dcfrac
               );

  ~SiPMDarkPdf();

  SiPMDarkPdf( const SiPMDarkPdf&, const char* name = 0 );
  virtual TObject* clone( const char* name ) const;

protected:
  RooRealProxy x;
  RooRealProxy ped;
  RooRealProxy gain;
  RooRealProxy s0;
  RooRealProxy s1;
  RooRealProxy dcfrac;

  mutable SiPMDarkFunc func;// Since evaluation is a const function

  double evaluate() const;
};


#endif