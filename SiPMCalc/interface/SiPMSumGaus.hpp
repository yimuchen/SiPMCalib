#ifndef SIPMCALIB_SIPMCALC_SIPMSUMGAUS_HPP
#define SIPMCALIB_SIPMCALC_SIPMSUMGAUS_HPP

#include "RooFunctorBinding.h"

class SiPMSumGaus : public RooFunctorPdfBinding
{
public:
  SiPMSumGaus( const char*, const char*,
           RooRealVar& x,
           RooRealVar& ped,
           RooRealVar& gain,
           RooRealVar& s0,
           RooRealVar& s1,
           RooRealVar& a0,
           RooRealVar& a1,
           RooRealVar& a2,
           RooRealVar& a3
           );
  SiPMSumGaus( const char*, const char*,
           RooRealVar& x,
           RooRealVar& ped,
           RooRealVar& gain,
           RooRealVar& s0,
           RooRealVar& s1,
           RooRealVar& a0,
           RooRealVar& a1,
           RooRealVar& a2
           );
  SiPMSumGaus( const char*, const char*,
           RooRealVar& x,
           RooRealVar& ped,
           RooRealVar& gain,
           RooRealVar& s0,
           RooRealVar& s1,
           RooRealVar& a0,
           RooRealVar& a1
           );
  ~SiPMSumGaus();
};

#endif