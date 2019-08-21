#include "SiPMCalib/SiPMCalc/interface/SiPMDarkPdf.hpp"
#include "SiPMCalib/SiPMCalc/interface/SiPMFormat.hpp"

#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/Common/interface/Maths.hpp"
#include "UserUtils/Common/interface/STLUtils/OStreamUtils.hpp"
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"
#include "UserUtils/MathUtils/interface/Measurement/Measurement.hpp"
#include "UserUtils/PlotUtils/interface/Ratio1DCanvas.hpp"

#include "RooChi2Var.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooRealVar.h"

#include <algorithm>
#include <fstream>

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
    ( "adcbin", usr::po::value<int>()->default_value( 8 ),
    "Area used for for binning" )
    ( "start", usr::po::value<int>(),
    "Sample to begin integration, leave blank for begining of sample" )
    ( "end", usr::po::value<int>(),
    "Sample to stop integration, leave blank for end of sample" )
    ( "maxarea", usr::po::value<double>(),
    "Maximum area for perform fit on, leave blank for auto determination" )
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
  const std::string est   = arg.ArgOpt<std::string>( "estcheck", "" );
  const double adcbin     = arg.Arg<int>( "adcbin" );
  const unsigned start    = arg.ArgOpt<int>( "start", 0 );
  const unsigned end      = arg.ArgOpt<int>( "end", 60 );

  SiPMFormat fmt( input, adcbin, start, end );
  fmt.RunDarkEstimate( est );

  const unsigned opend = std::min( end, fmt.PreSamples() + fmt.PostSamples() );

  const double max = arg.ArgOpt<double>( "maxarea",
    fmt.estped + 1.5* fmt.estgain );
  fmt.MakeDataSet( max );

  RooRealVar ped( "ped", "ped", -200, 200 );
  RooRealVar gain( "gain", "gain", 1, 500 );
  RooRealVar s0( "s0", "s0", 0.01, 100 );
  RooRealVar s1( "s1", "s1", 0.01, 50 );
  RooRealVar dcfrac( "dcfrac", "dcfrac", 0, 0.2 );
  RooRealVar epsilon( "epslion", "epsilon", 1e-5, 1e-1 );

  SiPMDarkPdf pdf( "dark", "dark",
                   fmt.x(), ped, gain, s0, s1, dcfrac, epsilon );

  // Running inital estimation
  ped    = fmt.estped;
  gain   = fmt.estgain;
  s0     = fmt.ests0;
  s1     = 5;
  dcfrac = fmt.estdcfrac;

  pdf.fitTo( fmt.data(),
    RooFit::Range( fmt.estped - 3*fmt.ests0, fmt.estped +  1.2*fmt.estgain ) );
  pdf.fitTo( fmt.data(),
    RooFit::Range( fmt.estped - 3*fmt.ests0, fmt.estped +  1.2*fmt.estgain ) );
  pdf.fitTo( fmt.data(),
    RooFit::Range( fmt.estped - 3*fmt.ests0, fmt.estped +  1.2*fmt.estgain ) );

  usr::plt::Ratio1DCanvas c( fmt.x() );
  auto& fitgraph = c.PlotPdf( pdf,
    RooFit::Normalization( fmt.data().sumEntries() ),
    usr::plt::EntryText( "Model Fit" ) );
  auto& datgraph = c.PlotData( fmt.data(),
    usr::plt::EntryText( "SiPM readout" ) );

  fitgraph.SetLineColor( kBlue );
  fitgraph.SetFillColor( kCyan );
  datgraph.SetMarkerSize( 0.2 );

  c.PlotScale( fitgraph, fitgraph,
    usr::plt::PlotType( usr::plt::scatter ) );
  c.PlotScale( datgraph, fitgraph,
    usr::plt::PlotType( usr::plt::scatter ) );

  // More information from fit parameters values
  const double window = ( opend-start ) * fmt.TimeInterval();
  const usr::Measurement pdc( dcfrac.getVal(),
                              dcfrac.getError(), dcfrac.getError() );
  const usr::Measurement tdc = window / 1000. / pdc;


  c.DrawLuminosity( "Random trigger" );
  c.DrawCMSLabel( "", "Spectral Fit" );
  c.TopPad()
  .WriteLine( usr::fstr( "Int. window = %d[ns]", window ) )
  .WriteLine( usr::fstr( "#tau_{DC} = %.2lf_{%.2lf}[us]",
    tdc.CentralValue(), tdc.AbsAvgError() ) );


  c.BottomPad().Yaxis().SetTitle( "Data/fit" );
  c.SetLogy( true );
  c.TopPad().SetYaxisMax( c.TopPad().GetYaxisMax() * 300 );
  c.SaveAsPDF( arg.MakePDFFile( "DarkSpectralFit" ) );

  const std::string sf = "%10s  %10.2lf %10.2lf %10.2lf\n";
  usr::fout( sf, "ped",     fmt.estped,    ped.getVal(),     ped.getError() );
  usr::fout( sf, "gain",    fmt.estgain,   gain.getVal(),    gain.getError() );
  usr::fout( sf, "s0 ",     fmt.ests0,     s0.getVal(),      s0.getError() );
  usr::fout( sf, "s1 ",     fmt.ests1,     s1.getVal(),      s1.getError() );
  usr::fout( sf, "dcfrac",  fmt.estdcfrac, dcfrac.getVal(),  dcfrac.getError() );
  usr::fout( sf, "epslion", 0,             epsilon.getVal(), epsilon.getError() );

  return 0;
}
