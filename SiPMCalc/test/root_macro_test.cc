#include "RooRealVar.h"
#include "TSystem.h"
#include "RooDataSet.h"

#include "SiPMCalib/SiPMCalc/interface/SiPMPdf.hpp"
#include <cstdlib>

void
root_macro_test()
{
  // This might need to be loaded before calling the macro file.
  // If you are running a script without the compiling flag
  const std::string libfile
    = std::getenv("CMSSW_BASE") + std::string("/lib/")
      + std::getenv("SCRAM_ARCH") + std::string("/libSiPMCalibSiPMCalc.so");
  gSystem->Load( libfile.c_str() );

  // Declaring the Variables.
  RooRealVar x( "x", "x", -100, 1600 );
  RooRealVar ped( "ped", "ped", -300, 300 );
  RooRealVar gain( "gain", "gain", 0, 1000 );
  RooRealVar s0( "s0", "s0", 0.001, 1000 );
  RooRealVar s1( "s1", "s1", 0, 100 );
  RooRealVar mean( "mean", "mean", 0.001, 50 );
  RooRealVar lambda( "lambda", "lambda", 0, 0.50 );
  RooRealVar alpha( "alpha", "alpha", 0, 0.5 );
  RooRealVar beta( "beta", "beta", 10, 20000 );
  RooRealVar dcfrac( "dcfrac", "dcfrac", 0, 0.4 );
  RooRealVar eps( "eps", "eps", 1e-5, 1e-1 );

  SiPMPdf p( "p", "p", x, ped, gain, s0, s1,
             mean, lambda,
             alpha, beta,
             dcfrac, eps  );

  RooDataSet data("data","data", RooArgSet(x) );

  p.fitTo( data) ;
}
