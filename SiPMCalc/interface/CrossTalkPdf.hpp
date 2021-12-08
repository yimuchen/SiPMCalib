#ifndef SIPMCALIB_SIPMCALC_CROSSTALKPDF_HPP
#define SIPMCALIB_SIPMCALC_CROSSTALKPDF_HPP

#include "SiPMCalib/SiPMCalc/interface/SiPMDarkFunc.hpp"

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooTrace.h"

class CrossTalkPdf : public RooAbsPdf
{
public:
  // Full model constructor
  CrossTalkPdf( const char*,
                const char*,
                RooRealVar& x,
                RooRealVar& x0,
                RooRealVar& gain,
                RooRealVar& s0,
                RooRealVar& s1,
                RooRealVar& prob );
  CrossTalkPdf( const CrossTalkPdf&, const char*name = 0 );
  virtual ~CrossTalkPdf();

  virtual TObject* clone( const char*name ) const;

  double gauss_k( const int k  ) const;
  double cross_prob( const int k ) const;
  inline double
  Eval() const { return evaluate(); }

protected:
  RooRealProxy x;
  RooRealProxy x0;
  RooRealProxy gain;
  RooRealProxy s0;
  RooRealProxy s1;
  RooRealProxy prob;

  double evaluate() const;

private:
//  ClassDef(CrossTalkPdf,1);
};

#endif
