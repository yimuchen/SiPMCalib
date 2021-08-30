#include "SiPMCalib/SiPMCalc/interface/SiPMPdf.hpp"
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"

#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooRealVar.h"

#include <fstream>
#include <iostream>

int
main( int argc, char* argv[] )
{
  // Parsing the input
  const std::string output = argv[1];
  const double nevents     = std::stof( argv[2] );
  const double mean_in     = std::stof( argv[3] );
  const double lamb_in     = std::stof( argv[4] );
  const double alph_in     = std::stof( argv[5] );
  const double s0_in       = std::stof( argv[6] );

  // The main variable
  RooRealVar x(      "x",       "x", -500, 10000 );

  // Operation parameters would typically be fixed at constant values
  RooConstVar ped(       "ped",    "ped",   0.00  );
  RooConstVar gain(     "gain",   "gain", 100.00  );
  RooConstVar lambda( "lambda", "lambda", lamb_in );
  RooConstVar alpha(   "alpha",  "alpha", alph_in );
  RooConstVar beta(     "beta",   "beta", 150.00  );
  RooConstVar dcfrac( "dcfrac", "dcfrac",   0.00001  );
  RooConstVar eps(       "eps",    "eps",   0.004 );

  RooConstVar mean(   "mean",   "mean",  mean_in       );
  RooConstVar s0(       "s0",     "s0",  100*s0_in     );
  RooConstVar s1(       "s1",     "s1",  1 );



  SiPMPdf pdf( "pdf", "pdf", x, ped, gain
             , s0, s1, mean, lambda, alpha, beta
             , dcfrac, eps );


  RooDataSet* data = pdf.generate( RooArgSet( x ), nevents );
  const double var = data->moment( x, 2 );
  const double m   = data->mean( x );
  delete data;
  const double enf = var / ( m*m ) * mean.getVal();

  // Outputting the results.
  std::ofstream f( output, std::ios::app );
  f << usr::fstr( "%8.6lf %8.6lf %8.6lf %8.6lf %8.6lf\n"
                , mean_in, lamb_in, alph_in,  s0_in, enf );
  std::cout << usr::fstr( "%8.6lf %8.6lf %8.6lf %8.6lf %8.6lf\n"
                        , mean_in, lamb_in, alph_in,  s0_in, enf );

  return 0;
}
