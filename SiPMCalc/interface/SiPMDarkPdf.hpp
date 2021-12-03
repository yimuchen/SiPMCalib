#ifndef SIPMCALIB_SIPMCALC_SIPMDARKPDF_HPP
#define SIPMCALIB_SIPMCALC_SIPMDARKPDF_HPP

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooTrace.h"

#include "SiPMCalib/SiPMCalc/interface/SiPMDarkFunc.hpp"

class SiPMDarkPdf : public RooAbsPdf
{
public:
  SiPMDarkPdf( const char*,
               const char*,
               RooAbsReal& x,
               RooAbsReal& ped,
               RooAbsReal& gain,
               RooAbsReal& s0,
               RooAbsReal& s1,
               RooAbsReal& dcfrac,
               RooAbsReal& epsilon );
  ~SiPMDarkPdf();
  SiPMDarkPdf( const SiPMDarkPdf&, const char* name = 0 );
  virtual TObject* clone( const char* name ) const;

  void RunEstimate( const RooAbsData&, const std::string& plot = "" );

protected:
  RooRealProxy x;
  RooRealProxy ped;
  RooRealProxy gain;
  RooRealProxy s0;
  RooRealProxy s1;
  RooRealProxy dcfrac;
  RooRealProxy epsilon;

  mutable MDistro mdistro;// Since evaluation is a const function

  double evaluate() const;
};


#endif
