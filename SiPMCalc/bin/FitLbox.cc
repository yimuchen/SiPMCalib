#include "SiPMCalib/SiPMCalc/interface/SiPMFormat.hpp"
#include "SiPMCalib/SiPMCalc/interface/SiPMPdf.hpp"
#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/Common/interface/STLUtils/OStreamUtils.hpp"
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"
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

  // Declaring fit variables
  SiPMFormat fmt( input, adcbin, start, end, base  );
  const double min = std::numeric_limits<double>::min();

  RooRealVar ped( "ped", "ped", -300, 300 );
  RooRealVar gain( "gain", "gain", 0, 1000 );
  RooRealVar s0( "s0", "s0", min, 1000 );
  RooRealVar s1( "s1", "s1", 0, 100 );
  RooRealVar mean( "mean", "mean", min, 50 );
  RooRealVar lambda( "lambda", "lambda", 0.01, 0.50 );
  RooRealVar alpha( "alpha", "alpha", 0, 0.5 );
  RooRealVar beta( "beta", "beta", 10, 20000 );
  RooRealVar dcfrac( "dcfrac", "dcfrac", 0, 0.4 );
  RooRealVar eps( "eps", "eps", 1e-5, 1e-1 );

  SiPMPdf p( "p", "p", fmt.x(), ped, gain, s0, s1,
             mean, lambda,
             alpha, beta,
             dcfrac, eps  );

  p.RunEstimate( fmt.data(), est );// Running initial estimation.

  const unsigned npeak = std::min( 10.0,
    std::max( 5.0, mean.getVal() + 2*sqrt( mean.getVal() ) ) );

  const double areamax =
    arg.ArgOpt<double>( "maxarea",
      ped.getVal() + npeak*gain.getVal() + 2* s0.getVal() );

  // Not running over all fit results.
  fmt.TruncateDataSet( areamax );
  alpha  = 0.05;
  beta   = 1000;
  dcfrac = arg.ArgOpt<double>( "dcfrac", 0 );
  eps    = arg.ArgOpt<double>( "epsilon", 1e-2 );
  dcfrac.setConstant( arg.CheckArg( "dcfrac" ) );
  eps.setConstant( arg.CheckArg( "epsilon" ) );

  p.fitTo( fmt.data() );
  p.fitTo( fmt.data() );


  usr::plt::Ratio1DCanvas c( fmt.x() );

  auto& fitgraph = c.PlotPdf( p,
    RooFit::Normalization( fmt.data().sumEntries(), RooAbsReal::NumEvent ),
    usr::plt::EntryText( "Fit" ),
    usr::plt::LineColor( usr::plt::col::blue ) );

  auto& datgraph = c.PlotData( fmt.data(),
    usr::plt::EntryText( "Data" ),
    usr::plt::MarkerSize( 0.2 ) );

  c.PlotScale( fitgraph, fitgraph,
    usr::plt::PlotType( usr::plt::simplefunc ) );

  c.PlotScale( datgraph, fitgraph,
    usr::plt::PlotType( usr::plt::scatter ) );

  // Constant calculations
  const unsigned opend   = std::min( end, fmt.PreSamples() + fmt.PostSamples() );
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
    .WriteLine( usr::fstr( "SiPM Bias: %.2lf V"
                         , arg.Arg<double>( "biasvoltage" ) ) );
  }

  c.BottomPad().Yaxis().SetTitle( "Data/Global fit" );
  c.SetLogy( true );
  c.TopPad().SetYaxisMax( c.TopPad().GetYaxisMax() );


  c.SaveAsPDF( arg.MakePDFFile( "Spectralfit" ) );

  const std::string rawfmt( "%10s %10.5lf %10.5lf\n" );
  usr::fout( "Raw Fit parameter results" );
  usr::fout( rawfmt, "ped",     ped.getVal(),    ped.getError() );
  usr::fout( rawfmt, "gain",    gain.getVal(),   gain.getError() );
  usr::fout( rawfmt, "s0",      s0.getVal(),     s0.getError() );
  usr::fout( rawfmt, "s1",      s1.getVal(),     s1.getError() );
  usr::fout( rawfmt, "mean",    mean.getVal(),   mean.getError() );
  usr::fout( rawfmt, "lambda",  lambda.getVal(), lambda.getError() );
  usr::fout( rawfmt, "alpha",   alpha.getVal(),  alpha.getError() );
  usr::fout( rawfmt, "beta",    beta.getVal(),   beta.getError() );
  usr::fout( rawfmt, "dcfrac",  dcfrac.getVal(), dcfrac.getError() );
  usr::fout( rawfmt, "epsilon", eps.getVal(),    eps.getError() );


  std::cout << "\n\n\n" << usr::separator( '-' ) << "\n\n\n" << std::endl;
  const std::string linefmt = "%.30s & = & %.20s & [\\text{%s}] \\\\\n";
  const std::string p1fmt   = "%.1lf & \\pm%.2lf";
  const std::string p2fmt   = "%.2lf & \\pm%.3lf";

  // Generating strings.
  const std::string mean_s = usr::fstr( p2fmt, mean.getVal(), mean.getError() );
  const std::string gain_s = usr::fstr( p2fmt
                                      , gain.getVal()  * fmt.ADCConversion()
                                      , gain.getError() * fmt.ADCConversion() );
  const std::string s0_s = usr::fstr( p2fmt
                                    , s0.getVal()  * fmt.ADCConversion()
                                    , s0.getError() * fmt.ADCConversion() );
  const std::string s1_s = usr::fstr( p2fmt
                                    , s1.getVal()  * fmt.ADCConversion()
                                    , s1.getError() * fmt.ADCConversion() );
  const std::string ct_s = usr::fstr( p2fmt,  xtalk_cen*100,  xtalk_err*100 );
  const std::string ap_s = usr::fstr( p2fmt
                                    , alpha.getVal()*100, alpha.getError()*100 );
  const std::string tap_s = usr::fstr( p2fmt
                                     , tap.CentralValue(), tap.AbsAvgError() );
  const std::string tdc_s = usr::fstr( p2fmt
                                     , tdc.CentralValue(), tdc.AbsAvgError() );
  const std::string mean_title = "\\left\\langle N_{\\gamma}\\right\\rangle";

  usr::fout( "Human readable fit results" );
  usr::fout( linefmt, mean_title,            mean_s, "Photons" );
  usr::fout( linefmt, "\\text{Gain}",        gain_s, "mV-ns" );
  usr::fout( linefmt, "\\sigma_\\text{com}", s0_s,   "mV-ns" );
  usr::fout( linefmt, "\\sigma_\\text{pix}", s1_s,   "mV-ns" );
  usr::fout( linefmt, "P_\\text{ct}",        ct_s,   "\\%" );
  usr::fout( linefmt, "P_\\text{ap}",        ap_s,   "\\%" );
  usr::fout( linefmt, "\\tau_\\text{ap}",    tap_s,  "ns" );
  usr::fout( linefmt, "\\tau_\\text{dc}",    tdc_s,  "\\mus" );


  return 0;
}
