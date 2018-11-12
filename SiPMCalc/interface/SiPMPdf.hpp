#ifndef SIPMCALIB_SIPMCALC_SIPMPDF_HPP
#define SIPMCALIB_SIPMCALC_SIPMPDF_HPP

#include "RooFunctorBinding.h"

class SiPMPdf : public RooFunctorPdfBinding
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
           RooRealVar& lambda
           );
  SiPMPdf( const char*, const char*,
           RooRealVar& x,
           RooRealVar& ped,
           RooRealVar& gain,
           RooRealVar& s0,
           RooRealVar& s1,
           RooRealVar& mean
           );
  ~SiPMPdf();

private:
  static ROOT::Math::IBaseFunctionMultiDim* _sipm_pdf;
};

#endif