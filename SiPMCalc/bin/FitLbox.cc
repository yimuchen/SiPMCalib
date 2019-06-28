#include "SiPMCalib/SiPMCalc/interface/SiPMFormat.hpp"
#include "SiPMCalib/SiPMCalc/interface/SiPMPdf.hpp"
#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/MathUtils/interface/Measurement.hpp"
#include "UserUtils/PlotUtils/interface/Ratio1DCanvas.hpp"

#include "RooChi2Var.h"
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooRealVar.h"

#include <algorithm>
#include <boost/format.hpp>
#include <fstream>
#include <iostream>


/// Getting better results:
// * Pedestal subtraction.
// * Hardware issues (double check readout requirements)

int
main( int argc, char* argv[] )
{
  usr::po::options_description desc( "Options for fitting data to lower light SiPM spectrum model" );
  desc.add_options()
    ( "data", usr::po::value<std::string>(),
    "Input data .txt file" )
    ( "output", usr::po::value<std::string>(),
    "Primary output file name" )
    ( "estcheck", usr::po::value<std::string>(),
    "Estimation output files, leave blank for no output" )
    ( "baseline", usr::po::value<std::string>(),
    "Baseline line .txt file" )
    ( "adcbin", usr::po::value<int>()->default_value( 16 ),
    "Area used for for binning" )
    ( "start", usr::po::value<int>(),
    "Sample to begin integration, leave blank for begining of sample" )
    ( "end", usr::po::value<int>(),
    "Sample to stop integration, leave blank for end of sample" )
    ( "maxarea", usr::po::value<double>(),
    "Maximum area for perform fit on" )
    ( "dcfrac", usr::po::value<double>(),
    "Externally obtained DC fraction to use, leave floating if not provided" )
    ( "epsilon", usr::po::value<double>(),
    "External factor to be used for dark current contribution cutoff, leave floating if not provided" )
    ( "biasvoltage", usr::po::value<double>(),
    "Bias voltage to be displayed on plot, not entry if left empty" )
    ( "lumitype", usr::po::value<std::string>()->default_value( "Laser setup" ),
    "Set up to be written on top right of the plot" )
  ;
  usr::ArgumentExtender arg;
  arg.AddOptions( desc );
  arg.ParseOptions( argc, argv );

  if( !arg.CheckArg( "data" ) ){
    std::cerr << "Please include input data" << std::endl;
    return 1;
  }
  if( !arg.CheckArg( "output" ) ){
    std::cerr << "Please specify output" << std::endl;
    return 1;
  }

  arg.SetNameScheme( {{"output", ""}} );

  const std::string input = arg.Arg( "data" );
  const std::string base  = arg.ArgOpt<std::string>( "baseline", "" );
  const std::string est   = arg.ArgOpt<std::string>( "estcheck", "" );
  const double adcbin     = arg.Arg<int>( "adcbin" );
  const unsigned start    = arg.ArgOpt<int>( "start", 0 );
  const unsigned end      = arg.ArgOpt<int>( "end", -1 );

  // Making objects
  SiPMFormat fmt( input, adcbin, start, end, base  );
  std::cout << base << std::endl;
  fmt.RunLumiEstimate( est );

  const unsigned opend = std::min( end, fmt.PreSamples() + fmt.PostSamples() );
  const unsigned npeak = std::min( 10.0,
    std::max( 5.0, fmt.estmean + 2*sqrt( fmt.estmean ) ) );

  const double areamax =
    arg.ArgOpt<double>( "maxarea",
      fmt.estped + npeak*fmt.estgain - 2* fmt.ests0 );

  fmt.MakeDataSet( areamax );


  const double min = std::numeric_limits<double>::min();

  RooRealVar ped( "ped", "ped", -300, 300 );
  RooRealVar gain( "gain", "gain", 0, 1000 );
  RooRealVar s0( "s0", "s0", min, 1000 );
  RooRealVar s1( "s1", "s1", 0, 100 );
  RooRealVar mean( "mean", "mean", min, 50 );
  RooRealVar lambda( "lambda", "lambda", 0, 0.50 );
  RooRealVar alpha( "alpha", "alpha", 0, 0.5 );
  RooRealVar beta( "beta", "beta", 10, 20000 );
  RooRealVar dcfrac( "dcfrac", "dcfrac", 0, 0.4 );
  RooRealVar eps( "eps", "eps", 1e-5, 1e-1 );

  // First pass, setting secondary effects to be constant.
  SiPMPdf p( "p", "p", fmt.x(), ped, gain, s0, s1,
             mean, lambda,
             alpha, beta,
             dcfrac, eps        );

  ped    = fmt.estped;
  gain   = fmt.estgain;
  s0     = fmt.ests0;
  s1     = fmt.ests1;
  mean   = fmt.estmean;
  lambda = fmt.estlambda;
  alpha  = 0;
  beta   = 1000;
  dcfrac = arg.ArgOpt<double>( "dcfrac", 0 );
  eps    = arg.ArgOpt<double>( "epsilon", 1e-2 );

  ped.setConstant( true );
  gain.setConstant( true );
  s0.setConstant( true );
  mean.setConstant( true );
  dcfrac.setConstant( true );
  eps.setConstant( true );
  alpha = 0.05;

  p.fitTo( fmt.data() );

  ped.setConstant( false );
  gain.setConstant( false );
  s0.setConstant( false );
  mean.setConstant( false );
  dcfrac.setConstant( arg.CheckArg( "dcfrac" ) );
  eps.setConstant( arg.CheckArg( "epsilon" ) );

  p.fitTo( fmt.data() );
  p.fitTo( fmt.data() );


  usr::plt::Ratio1DCanvas c( fmt.x() );

  auto& fitgraph = c.PlotPdf( p,
    RooFit::Normalization( fmt.data().sumEntries() ),
    usr::plt::EntryText( "Fit" ) );

  auto& datgraph = c.PlotData( fmt.data(),
    usr::plt::EntryText( "Data" ) );

  fitgraph.SetLineColor( kBlue );
  fitgraph.SetFillColor( kCyan );
  datgraph.SetMarkerSize( 0.2 );


  c.PlotScale( fitgraph, fitgraph,
    usr::plt::PlotType( usr::plt::simplefunc ) );

  c.PlotScale( datgraph, fitgraph,
    usr::plt::PlotType( usr::plt::scatter ) );

  // Constant calculations
  const double window    = ( opend-start ) * fmt.TimeInterval();
  const double lval      = lambda.getVal();
  const double lerr      = lambda.getError();
  const double xtalk_cen = 1 - TMath::Exp( -lval );
  const double xtalk_hi  = 1 - TMath::Exp( -( lval + lerr ) );
  const double xtalk_lo  = 1 - TMath::Exp( -( lval - lerr ) );
  const double xtalk_err = ( xtalk_hi - xtalk_lo )/2;

  const usr::Measurement pdc( dcfrac.getVal(),
                              dcfrac.getError(), dcfrac.getError() );
  const usr::Measurement tdc = window / 1000. / pdc;
  const usr::Measurement mb( beta.getVal(), beta.getError(), beta.getError() );
  const usr::Measurement mg( gain.getVal(), gain.getError(), gain.getError() );
  const usr::Measurement tap = ( 30. / mg ) * mb;

  c.TopPad().DrawLuminosity( arg.Arg( "lumitype" ) );
  c.TopPad().DrawCMSLabel( "", "Spectral fit" );
  if( arg.CheckArg( "biasvoltage" ) ){
    c.TopPad()
    .WriteLine( ( boost::format( "SiPM Bias: %.2lf V" )
                  % arg.Arg<double>( "biasvoltage" ) ).str() );
  }

  c.BottomPad().Yaxis().SetTitle( "Data/Global fit" );
  c.SetLogy( true );
  c.TopPad().SetYaxisMax( c.TopPad().GetYaxisMax() );


  c.SaveAsPDF( arg.MakePDFFile( "Spectralfit" ) );

  boost::format linefmt( "%s & = & %s & [\\text{%s}] \\\\" );
  boost::format p1fmt( "%.1lf & \\pm%.2lf" );
  boost::format p2fmt( "%.2lf & \\pm%.3lf" );

  std::cout
    << ( linefmt
         % "\\left\\langle N_{\\gamma}\\right\\rangle"
         % ( p2fmt % mean.getVal() % mean.getError() ).str()
         % "Photons"  )<< std::endl
    << ( linefmt
       % "\\text{Gain}"
       % ( p2fmt
           % ( gain.getVal()  * fmt.ADCConversion() )
           % ( gain.getError() * fmt.ADCConversion() ) ).str()
       % "mV-ns" )  << std::endl
    << ( linefmt
       % "\\sigma_\\text{com}"
       % ( p2fmt
           % ( s0.getVal()  * fmt.ADCConversion() )
           % ( s0.getError() * fmt.ADCConversion() ) ).str()
       % "mV-ns" )  << std::endl
    << ( linefmt
       % "\\sigma_\\text{pix}"
       % ( p2fmt
           % ( s1.getVal()  * fmt.ADCConversion() )
           % ( s1.getError() * fmt.ADCConversion() ) ).str()
       % "mV-ns" ) << std::endl
    << ( linefmt
       % "P_\\text{ct}"
       % ( p2fmt % ( xtalk_cen*100 ) % ( xtalk_err*100 ) ).str()
       % "\\%" )  << std::endl
    << ( linefmt
       % "P_\\text{ap}"
       % ( p2fmt % ( alpha.getVal()*100 ) % ( alpha.getError()*100 ) ).str()
       % "\\%" )  << std::endl
    << ( linefmt
       % "\\tau_\\text{ap}"
       % ( p2fmt % ( tap.CentralValue() ) % tap.AbsAvgError() ).str()
       % "ns" )  << std::endl
    << ( linefmt
       % "\\tau_\\text{dc}"
       % ( p2fmt % ( tdc.CentralValue() ) % tdc.AbsAvgError() ).str()
       % "ns" )  << std::endl
  ;

  std::cout
    << "ped " << ped.getVal() <<  "  " << ped.getError() << std::endl
    << "gain " << gain.getVal() <<  "  " << gain.getError() << std::endl
    << "s0 " << s0.getVal() <<  "  " << s0.getError() << std::endl
    << "s1 " << s1.getVal() <<  "  " << s1.getError() << std::endl
    << "mean " << mean.getVal() <<  "  " << mean.getError() << std::endl
    << "lambda " << lambda.getVal() <<  "  " << lambda.getError() << std::endl
    << "alpha " << alpha.getVal() <<  "  " << alpha.getError() << std::endl
    << "beta " << beta.getVal() << "  " << beta.getError() << std::endl
    << "dcfrac " << dcfrac.getVal() <<  "  " << dcfrac.getError() << std::endl
    << "epsilon " << eps.getVal() <<  "  " << eps.getError() << std::endl
  ;

  return 0;
}
