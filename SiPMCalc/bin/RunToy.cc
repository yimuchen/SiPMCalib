#include <SiPMCalib/SiPMCalc/interface/SiPMDarkPdf.hpp>
#include <SiPMCalib/SiPMCalc/interface/SiPMPdf.hpp>
#include <SiPMCalib/SiPMCalc/interface/ToyRunCommon.hpp>
#include <UserUtils/Common/interface/ArgumentExtender.hpp>
#include <UserUtils/Common/interface/STLUtils/Filesystem.hpp>
#include <UserUtils/Common/interface/SystemUtils/Time.hpp>

#include "RooConstVar.h"
#include "RooDataHist.h"
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
                     ComNoise( arg ), ComNoise( arg )/4, ComNoise( arg )*4 );
  RooRealVar s1(     "s1",     "s1",
                     PixNoise( arg ), PixNoise( arg )/4, PixNoise( arg )*4 );
  RooRealVar mean(   "mean",   "mean",
                     Mean( arg ),         Mean( arg )/8,     Mean( arg )*8 );
  RooRealVar lambda( "lambda",    "lambda",
                     Lambda( arg ),       0, std::min( Lambda( arg )*4, 0.5 ) );
  RooRealVar alpha(  "alpha",  "alpha",
                     Alpha( arg ),        0, std::min( Alpha( arg )*4, 0.5 ) );
  RooRealVar beta(   "beta",   "beta",
                     Beta( arg ),        Beta( arg )/16,    Beta( arg )*16 );
  RooRealVar dcfrac( "dcfrac", "dcfrac",
                     DCFrac( arg ),       0, std::min( DCFrac( arg )*4, 0.6 ) );
  RooRealVar epsilon( "epsilon", "epsilon",
                      0.01, 1e-6, 1e-1 );


  x.setBin( 20*( Xmax( arg ) - Xmin( arg ) )/( Gain( arg ) ) );

  // Defining PDF
  SiPMPdf full( "full", "full",
                x, ped, gain, s0, s1, mean, lambda,
                alpha, beta, dcfrac, epsilon );
  SiPMPdf ndc( "ndc", "ndc",
               x, ped, gain, s0, s1, mean, lambda, alpha, beta );
  SiPMPdf simp( "simp", "simp",
                x, ped, gain, s0, s1, mean, lambda  );
  SiPMDarkPdf dark( "dark", "dark",
                    x, ped, gain, s0, s1, dcfrac, RooFit::RooConst( 1e-3 ) );

  RooAbsPdf* pdf
    = arg.Arg( "model" ) == "dark" ? dynamic_cast<RooAbsPdf*>( &dark ) :
      arg.Arg( "model" ) == "simp" ? dynamic_cast<RooAbsPdf*>( &simp ) :
      arg.Arg( "model" ) == "ndc"  ? dynamic_cast<RooAbsPdf*>( &ndc ) :
      dynamic_cast<RooAbsPdf*>( &full );

  RooMsgService::instance().setSilentMode( true );
  RooRandom::randomGenerator()->SetSeed( usr::CurrentTimeInNanSec() );

  std::ofstream fout;
  fout.open( filename( arg ), std::ios::app );

  // Begin generation loop
  for( int i = 0; i < arg.Arg<int>( "nToys" ); ++i ){

    std::cout << boost::format( "\rRunning Toy %d/%d [%s]..." )
      % ( i+1 ) % arg.Arg<int>( "nToys" )
      % usr::CurrentTime() << std::flush;

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
        dynamic_cast<RooAbsData*>(
          pdf->generateBinned( RooArgSet( x ), arg.Arg<int>( "nEvents" ) ) ) :
        dynamic_cast<RooAbsData*>(
          pdf->generate( RooArgSet( x ), arg.Arg<int>( "nEvents" ) ) );

    RooFitResult* fit = nullptr;
    unsigned iter     = 0;

    const uint64_t tstart = usr::CurrentTimeInMuSec();

    while( iter <= 5 ){
      iter++;

      fit = pdf->fitTo( *dat,
        RooFit::Verbose( kFALSE ),
        RooFit::PrintLevel( -100 ),
        RooFit::PrintEvalErrors( -100 ),
        RooFit::Warnings( kFALSE ),
        RooFit::Save()
        );

      if( fit->status() && iter <= 5 ){
        delete fit;
      }
    }

    const uint64_t tend = usr::CurrentTimeInMuSec();

    boost::format fmt( "%14.8lf  %14.9lf  " );

    fout << fit->status() << "  ";
    fout << ( fmt % ped.getVal()    % ped.getError() );
    fout << ( fmt % gain.getVal()   % gain.getError() );
    fout << ( fmt % s0.getVal()     % s0.getError() );
    fout << ( fmt % s1.getVal()     % s1.getError() );
    fout << ( fmt % mean.getVal()   % mean.getError() );
    fout << ( fmt % lambda.getVal() % lambda.getError() );
    fout << ( fmt % dcfrac.getVal() % dcfrac.getError() );
    fout << ( fmt % alpha.getVal()  % alpha.getError() );
    fout << ( fmt % beta.getVal()   % beta.getError() );
    fout << tend - tstart << std::endl;

    std::cout << boost::format( "[%s] Done! (%ldsec)" )
      % usr::CurrentTime()
      % ( ( tend - tstart ) / 1e6 )
              << std::endl;

    delete fit;
    delete dat;
  }

  fout.close();

  return 0;
}
