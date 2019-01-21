#include <SiPMCalib/SiPMCalc/interface/SiPMPdf.hpp>
#include <UserUtils/Common/interface/ArgumentExtender.hpp>
#include <UserUtils/Common/interface/STLUtils/Filesystem.hpp>
#include <UserUtils/Common/interface/SystemUtils/Time.hpp>

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooRandom.h"
#include "RooFitResult.h"

#include <boost/format.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>

int
main( int argc, char** argv )
{
  usr::po::options_description desc( "Toy MC options" );
  desc.add_options()
    ( "nEvents,n", usr::po::value<int>(), "Number of events in a toy data set" )
    ( "nToys", usr::po::value<int>(), "Number of toy models to run" )
    ( "model,m", usr::po::value<std::string>(), "Model to use run the MC studies" )
  ;

  usr::ArgumentExtender arg( "data/modelcfg.json" );
  arg.AddOptions( desc );
  arg.ParseOptions( argc, argv );

  RooRealVar x( "x", "Readout (mV#times ns)", 300, 1500 );
  RooRealVar ped( "ped", "ped", 350, 100, 500 );
  RooRealVar gain( "gain", "gain", 200, 100, 500 );
  RooRealVar s0( "s0", "s0", 50, 15, 100 );
  RooRealVar s1( "s1", "s1", 10, 0.05, 100 );
  RooRealVar mean( "mean", "mean", 1.5, 0.01, 50 );
  RooRealVar lambda( "lambda", "lambda", 0.25, 0, 0.2 );
  RooRealVar acfrac("acfrac", "acfrac", 0 , 0 , 1 );
  RooRealVar acshift("acshift","acshift",0,0,1);
  RooRealVar acwidth("acwidth","acwidth",1,0.5,10);
  RooRealVar alpha( "alpha", "alpha", 0.1, 0, 1 );
  RooRealVar beta( "beta", "beta", 50, 20, 1000 );

  SiPMPdf full( "full", "full",
                x, ped, gain, s0, s1, mean, lambda,
                acfrac, acshift, acwidth,
                alpha, beta );
  SiPMPdf pgaus( "pgaus", "pgaus",
                 x, ped, gain, s0, s1, mean, lambda,
                 acfrac, acshift, acwidth );
  SiPMPdf gpgaus( "gpgaus", "gpgaus",
                  x, ped, gain, s0, s1, mean, lambda );

  RooAbsPdf& oppdf = arg.Arg<std::string>( "model" ) == "full" ? full :
                     arg.Arg<std::string>( "model" ) == "gpgaus" ? gpgaus :
                     arg.Arg<std::string>( "model" ) == "pgaus" ? pgaus :
                     full;

  const std::string filename = ( boost::format( "%s_%d.txt" )
                                 % arg.Arg<std::string>( "model" )
                                 % arg.Arg<int>( "nEvents" ) ).str();
  const usr::fs::path outputpath
    = usr::resultpath( "SiPMCalib", "SiPMCalc" ) / filename;

  if( !usr::fs::exists( outputpath ) ){
    std::cout << "Creating new file!" << std::endl;
    std::ofstream outst;
    outst.open( outputpath, std::ios::out );
    outst << std::setprecision( 6 )
          << ped.getVal() << " "
          << gain.getVal() << " "
          << s0.getVal() << " "
          << s1.getVal() << " "
          << mean.getVal() << " ";
    if( arg.Arg<std::string>( "model" ) == "gpgaus"
        || arg.Arg<std::string>( "model" ) == "full" ){
      outst << lambda.getVal() << " ";
    }
    if( arg.Arg<std::string>( "model" ) == "full" ){
      outst << alpha.getVal() << " ";
      outst << beta.getVal() << " ";
    }
    outst << std::endl;
    outst.close();
  }


  const double ped_c    = ped.getVal();
  const double gain_c   = gain.getVal();
  const double s0_c     = s0.getVal();
  const double s1_c     = s1.getVal();
  const double mean_c   = mean.getVal();
  const double lambda_c = lambda.getVal();
  const double alpha_c  = alpha.getVal();
  const double beta_c   = beta.getVal();

  RooRandom::randomGenerator()->SetSeed( usr::CurrentTimeInNanSec());

  std::ofstream outst;
  outst.open( outputpath, std::ios::app );

  for( int i = 0; i < arg.Arg<int>( "nToys" ); ++i ){
    ped    = ped_c;
    gain   = gain_c;
    s0     = s0_c;
    s1     = s1_c;
    mean   = mean_c;
    lambda = lambda_c;
    alpha  = alpha_c;
    beta   = beta_c;

    RooDataSet* dat = oppdf.generate(RooArgSet(x), arg.Arg<int>("nEvents") );

    RooFitResult* ans = nullptr;
    unsigned iter = 0;
    while( !ans && iter < 5 ){
      ans = oppdf.fitTo( *dat,
        RooFit::Verbose( kFALSE ),
        RooFit::PrintLevel( -1 ),
        RooFit::PrintEvalErrors( -1 ),
        RooFit::Warnings( kFALSE ),
        RooFit::Save()
        );

      if( ans->status() ){
        delete ans;
        ans = nullptr;
        ++iter;
      }
    }
    delete ans;

    outst << std::setprecision( 6 )
          << ped.getVal() << " " << ped.getError() << " "
          << gain.getVal() << " " << gain.getError() << " "
          << s0.getVal() << " " << s0.getError() << " "
          << s1.getVal() << " " << s1.getError() << " "
          << mean.getVal() << " " << mean.getError() << " ";
    if( arg.Arg<std::string>( "model" ) == "gpgaus"
        || arg.Arg<std::string>( "model" ) == "full" ){
      outst << lambda.getVal() << " " << lambda.getError() << " ";
    }
    if( arg.Arg<std::string>( "model" ) == "full" ){
      outst << alpha.getVal() << " " << alpha.getError() << " ";
      outst << beta.getVal() << " " << beta.getError() << " ";
    }
    outst << std::endl;

    delete dat;
  }

  outst.close();

  return 0;
}
