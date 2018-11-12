#include "SiPMCalib/SiPMCalc/interface/SiPMPdf.hpp"
#include "SiPMCalib/SiPMCalc/interface/SiPMPdf_Functor.hpp"

#include "RooRealVar.h"
#include "RooConstVar.h"

ROOT::Math::IBaseFunctionMultiDim* SiPMPdf::_sipm_pdf = new SiPMPdf_Functor();

SiPMPdf::SiPMPdf( const char* name, const char* title,
                  RooRealVar& x,
                  RooRealVar& ped,
                  RooRealVar& gain,
                  RooRealVar& s0,
                  RooRealVar& s1,
                  RooRealVar& mean,
                  RooRealVar& lambda,
                  RooRealVar& alpha,
                  RooRealVar& beta
                  ) :
  RooFunctorPdfBinding(
    name, title, *_sipm_pdf,
    RooArgList( x, ped, gain, s0, s1, mean, lambda, alpha, beta ) ){}

SiPMPdf::SiPMPdf( const char* name, const char* title,
                  RooRealVar& x,
                  RooRealVar& ped,
                  RooRealVar& gain,
                  RooRealVar& s0,
                  RooRealVar& s1,
                  RooRealVar& mean,
                  RooRealVar& lambda
                  ) :
  RooFunctorPdfBinding(
    name, title, *_sipm_pdf,
    RooArgList(
      x, ped, gain, s0, s1, mean, lambda,
      RooFit::RooConst( 0 ),
      RooFit::RooConst( x.getMax() * 1e3 ) ) ){}

SiPMPdf::SiPMPdf( const char* name, const char* title,
                  RooRealVar& x,
                  RooRealVar& ped,
                  RooRealVar& gain,
                  RooRealVar& s0,
                  RooRealVar& s1,
                  RooRealVar& mean
                  ):
  RooFunctorPdfBinding(
    name, title, *_sipm_pdf,
    RooArgList(
      x, ped, gain, s0, s1, mean,
      RooFit::RooConst( 0 ),
      RooFit::RooConst( 0 ),
      RooFit::RooConst( x.getMax() * 1e3 ) ) ){}
SiPMPdf::~SiPMPdf(){}
