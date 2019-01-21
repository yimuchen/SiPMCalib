#ifndef SIPMCALIB_SIPMCALC_SIPMPDF_HPP
#define SIPMCALIB_SIPMCALC_SIPMPDF_HPP

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooTrace.h"
#include <vector>

class SiPMPdf : public RooAbsPdf
{
public:
  SiPMPdf( const char*, const char*,
           RooRealVar& x,
           RooRealVar& ped,
           RooRealVar& gain,
           RooRealVar& s0,
           RooRealVar& s1,
           RooRealVar& mean,
           RooRealVar& lambda,
           RooRealVar& acfrac,
           RooRealVar& acshift,
           RooRealVar& acwidth,
           RooRealVar& alpha,
           RooRealVar& beta
           );
  SiPMPdf( const char*, const char*,
           RooRealVar& x,
           RooRealVar& ped,
           RooRealVar& gain,
           RooRealVar& s0,
           RooRealVar& s1,
           RooRealVar& mean,
           RooRealVar& lambda,
           RooRealVar& alpha,
           RooRealVar& beta
           );
  SiPMPdf( const char*, const char*,
           RooRealVar& x,
           RooRealVar& ped,
           RooRealVar& gain,
           RooRealVar& s0,
           RooRealVar& s1,
           RooRealVar& mean,
           RooRealVar& lambda,
           RooRealVar& acfrac,
           RooRealVar& acshift,
           RooRealVar& acawidth
           );
  SiPMPdf( const char*, const char*,
           RooRealVar& x,
           RooRealVar& ped,
           RooRealVar& gain,
           RooRealVar& s0,
           RooRealVar& s1,
           RooRealVar& mean,
           RooRealVar& lambda
           );

  SiPMPdf( const SiPMPdf&, const char* name = 0 );

  virtual
  ~SiPMPdf();

  virtual TObject* clone( const char* name ) const;

protected:
  RooRealProxy x;
  RooRealProxy ped;
  RooRealProxy gain;
  RooRealProxy s0;
  RooRealProxy s1;
  RooRealProxy mean;
  RooRealProxy lambda;
  RooRealProxy acfrac;
  RooRealProxy acshift;
  RooRealProxy acwidth;
  RooRealProxy alpha;
  RooRealProxy beta;

  double evaluate() const;

  double gen_poisson( const int k ) const;
  double ap_eff( const int k, const int i ) const;
  double gauss_k( const int k  ) const;
  double binomial_prob( const int k, const int i ) const;

private:
//  ClassDef(SiPMPdf,1);
};

#endif
