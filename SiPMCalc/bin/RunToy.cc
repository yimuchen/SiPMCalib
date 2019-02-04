#include <SiPMCalib/SiPMCalc/interface/SiPMPdf.hpp>
#include <UserUtils/Common/interface/ArgumentExtender.hpp>
#include <UserUtils/Common/interface/STLUtils/Filesystem.hpp>
#include <UserUtils/Common/interface/SystemUtils/Time.hpp>

#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooRandom.h"
#include "RooRealVar.h"

#include <boost/format.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>

int
main( int argc, char** argv )
{
  // Options defined in ToyRunCommon.cc
  usr::ArgumentExtender arg( "data/modelcfg.json" );
  arg.AddOptions( ToyOptions() );
  arg.ParseOptions( argc, argv );

  // Setting up default value and range for fitting
  RooRealVar x( "x", "Readout (mV#times ns)", Xmin( arg ), Xmax( arg ) );

  RooRealVar ped(    "ped",    "ped",
                     Pedestal( arg ),      -Gain( arg ),     Gain( arg ) );
  RooRealVar gain(   "gain",   "gain",
                     Gain( arg ),         Gain( arg )/2,     Gain( arg )*2 );
  RooRealVar s0(     "s0",     "s0",
                     ComNoise( arg ), ConNoise( arg )/4, ComNoise( arg )*4 );
  RooRealVar s1(     "s1",     "s1",
                     PixNoise( arg ), PixNoise( arg )/4, PixNoise( arg )*4 );
  RooRealVar mean(   "mean",   "mean",
                     Mean( arg ),         Mean( arg )/8,     Mean( arg )*8 );
  RooRealVar lambda( "lambda",    "lambda",
                     Lambda( arg ),       0, std::min( Lambda( arg )*4, 0.5 ) );
  RooRealVar dcfrac( "dcfrac", "dcfrac",
                     DCFrac( arg ), 0, std::min( DCFrac( arg )*4, 0.6 ) );
  RooRealVar alpha(  "alpha",  "alpha",
                     Alpha( arg ),  0, std::min( Alpha( arg )*4, 0.5 ) );
  RooRealVar beta(   "beta",   "beta",
                     Beta( arg ),        Beta( arg )/16,    Beta( arg )*16 );

  x.setBin( 16*( Xmax( arg ) - Xmin( arg ) )/( Gain( arg ) ) );

  // Defining PDF
  SiPMPdf full( "full", "full",
                x, ped, gain, s0, s1, mean, lambda, dcfrac, alpha, beta );
  SiPMPdf ndc( "ndc", "ndc",
               x, ped, gain, s0, s1, mean, lambda, alpha, beta );
  SiPMPdf nap( "nap", "nap",
               x, ped, gain, s0, s1, mean, lambda, dcfrac );
  SiPMPdf simp( "simp", "simp",
                x, ped, gain, s0, s1, mean, lambda  );
  SiPMDarkPdf dark( "dark", "dark",
                    x, ped, gain, s0, s1, dcfrac );

  RooAbsPdf& pdf = arg.Arg( "model" ) == "dark" ? dark :
                   arg.Arg( "model" ) == "simp" ? simp :
                   arg.Arg( "model" ) == "nap"  ? nap  :
                   arg.Arg( "model" ) == "ndc"  ? ndc  :
                   full;

  RooRandom::randomGenerator()->SetSeed( usr::CurrentTimeInNanSec() );

  std::ofstream fout;
  fout.open( filename( arg ), std::ios::app );

  // Begin generation loop
  for( int i = 0; i < arg.Arg<int>( "nToys" ); ++i ){
    ped    = Pedestal( arg );
    gain   = Gain( arg );
    s0     = ComNoise( arg );
    s1     = PixNoise( arg );
    mean   = Mean( arg );
    lambda = Lambda( arg );
    dcfrac = DCFrac( arg );
    alpha  = Alpha( arg );
    beta   = Beta( arg );

    RooAbsData* dat
      = arg.Arg( "fit" ) == "binned" ?
        pdf.generateBinned( RooArgSet( x ), arg.Arg<int>( "nEvents" ) ) :
        pdf.generate( RooArgSet( x ), arg.Arg<int>( "nEvents" ) );

    RooFitResult* fit = nullptr;
    unsigned iter     = 0;

    while( !fit && iter < 5 ){
      fit = oppdf.fitTo( *dat,
        RooFit::Verbose( kFALSE ),
        RooFit::PrintLevel( -1 ),
        RooFit::PrintEvalErrors( -1 ),
        RooFit::Warnings( kFALSE ),
        RooFit::Save()
        );

      if( fit->status() ){
        delete fit;
        fit = nullptr;
        ++iter;
      }
    }

    delete ans;
    delete dat;

    boost::format fmt( "%14.8lf  %14.9lf  " );

    fout << fmt % ped.getVal()    % pdf.getError()
         << fmt % gain.getVal()   % gain.getError()
         << fmt % s0.getVal()     % s0.getError()
         << fmt % s1.getVal()     % s1.getError()
         << fmt % mean.getVal()   % mean.getError()
         << fmt % lambda.getVal() % lambda.getError()
         << fmt % dcfrac.getVal() % dcfrac.getError()
         << fmt % alpha.getVal()  % alpha.getError()
         << fmt % beta.getVal()   % beta.getError()
         << std::endl;
  }

  fout.close();

  return 0;
}
