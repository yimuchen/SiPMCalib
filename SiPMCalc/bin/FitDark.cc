#include "SiPMCalib/Common/interface/MakeRooData.hpp"
#include "SiPMCalib/Common/interface/WaveFormat.hpp"
#include "SiPMCalib/SiPMCalc/interface/SiPMDarkPdf.hpp"

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
    ( "data", usr::po::reqvalue<std::string>(),
    "Input data .txt file" )
    ( "output", usr::po::reqvalue<std::string>(),
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

  arg.SetNameScheme( {{"output", ""}} );

  const std::string input = arg.Arg( "data" );
  const std::string est   = arg.ArgOpt<std::string>( "estcheck", "" );
  const double adcbin     = arg.Arg<int>( "adcbin" );
  const unsigned start    = arg.ArgOpt<int>( "start", 0 );
  const unsigned end      = arg.ArgOpt<int>( "end", 60 );

  WaveFormat wformat( input );

  RooRealVar x( "x", "x", -1000, 10000 );
  RooRealVar ped( "ped", "ped", -200, 200 );
  RooRealVar gain( "gain", "gain", 1, 500 );
  RooRealVar s0( "s0", "s0", 0.01, 100 );
  RooRealVar s1( "s1", "s1", 0.01, 50 );
  RooRealVar dcfrac( "dcfrac", "dcfrac", 0, 0.2 );
  RooRealVar epsilon( "epslion", "epsilon", 1e-5, 1e-1 );
  SiPMDarkPdf pdf( "dark", "dark",
                   x, ped, gain, s0, s1, dcfrac, epsilon );

  const auto list = wformat.SumList( start, end );

  SetRange( x, adcbin, -1, list );
  std::unique_ptr<RooDataHist> data( MakeData( x, list, -1 ) );
  pdf.RunEstimate( *data );

  const double max = arg.ArgOpt<double>( "maxarea",
    ped.getVal() + 1.5*gain.getVal() );
  data.reset( MakeData( x, list, max ) );

  pdf.fitTo( *data,
    RooFit::Range( ped.getVal() - 3*s0.getVal()
                 , ped.getVal() + 1.2*gain.getVal() ) );
  pdf.fitTo( *data,
    RooFit::Range( ped.getVal() - 3*s0.getVal()
                 , ped.getVal() + 1.2*gain.getVal() ) );
  pdf.fitTo( *data,
    RooFit::Range( ped.getVal() - 3*s0.getVal()
                 , ped.getVal() + 1.2*gain.getVal() ) );

  usr::plt::Ratio1DCanvas c( x );
  auto& fitgraph = c.PlotPdf( pdf,
    RooFit::Normalization( data->sumEntries() ),
    usr::plt::LineColor( usr::plt::col::blue ),
    usr::plt::FillColor( usr::plt::col::cyan ),
    usr::plt::EntryText( "Model Fit" ) );
  auto& datgraph = c.PlotData( *data,
    usr::plt::EntryText( "SiPM readout" ),
    usr::plt::MarkerSize( 0.2 ) );

  c.PlotScale( fitgraph, fitgraph,
    usr::plt::PlotType( usr::plt::scatter ) );
  c.PlotScale( datgraph, fitgraph,
    usr::plt::PlotType( usr::plt::scatter ) );

  // More information from fit parameters values
  const double window = ( end-start ) * wformat.Time();
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

  const std::string sf = "%10s %10.2lf %10.2lf\n";
  usr::fout( sf, "ped",     ped.getVal(),     ped.getError() );
  usr::fout( sf, "gain",    gain.getVal(),    gain.getError() );
  usr::fout( sf, "s0 ",     s0.getVal(),      s0.getError() );
  usr::fout( sf, "s1 ",     s1.getVal(),      s1.getError() );
  usr::fout( sf, "dcfrac",  dcfrac.getVal(),  dcfrac.getError() );
  usr::fout( sf, "epslion", epsilon.getVal(), epsilon.getError() );
  return 0;
}
