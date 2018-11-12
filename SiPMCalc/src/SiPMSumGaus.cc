#include "SiPMCalib/SiPMCalc/interface/SiPMSumGaus.hpp"
#include "SiPMCalib/SiPMCalc/interface/SiPMSumGaus_Functor.hpp"

#include "RooRealVar.h"
SiPMSumGaus::SiPMSumGaus( const char* name, const char* title,
                          RooRealVar& x,
                          RooRealVar& ped,
                          RooRealVar& gain,
                          RooRealVar& s0,
                          RooRealVar& s1,
                          RooRealVar& a0,
                          RooRealVar& a1,
                          RooRealVar& a2,
                          RooRealVar& a3
                          ) :
  RooFunctorPdfBinding( name, title, SiPMSumGaus_Functor( 4 ),
                        RooArgList( x, ped, gain, s0, s1, a0, a1, a2, a3 )
                        ){}
SiPMSumGaus::SiPMSumGaus( const char* name, const char* title,
                          RooRealVar& x,
                          RooRealVar& ped,
                          RooRealVar& gain,
                          RooRealVar& s0,
                          RooRealVar& s1,
                          RooRealVar& a0,
                          RooRealVar& a1,
                          RooRealVar& a2
                          ) :
  RooFunctorPdfBinding( name, title, SiPMSumGaus_Functor( 3 ),
                        RooArgList( x, ped, gain, s0, s1, a0, a1, a2 )
                        ){}
SiPMSumGaus::SiPMSumGaus( const char* name, const char* title,
                          RooRealVar& x,
                          RooRealVar& ped,
                          RooRealVar& gain,
                          RooRealVar& s0,
                          RooRealVar& s1,
                          RooRealVar& a0,
                          RooRealVar& a1
                          ) :
  RooFunctorPdfBinding( name, title, SiPMSumGaus_Functor( 2 ),
                        RooArgList( x, ped, gain, s0, s1, a0, a1 )
                        ){}

SiPMSumGaus::~SiPMSumGaus(){}
