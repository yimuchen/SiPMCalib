#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/PlotUtils/interface/PlotCommon.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

#include <fstream>
#include <sstream>
#include <string>

#include "TF1.h"
#include "TH1D.h"
#include "TSpectrum.h"

int
main( int argc, char* argv[] )
{
  usr::po::options_description desc( "Options for Pulser shape plotting" );
  desc.add_options()
    ( "data", usr::po::value<std::string>(),   "Input data .txt file" )
    ( "output", usr::po::value<std::string>(), "Primary output file name" )
    ( "min", usr::po::value<double>(), "X axis minimum" )
    ( "max", usr::po::value<double>(), "X axis maximum" )
    ( "nbins", usr::po::value<unsigned>(), "number of bins" )
  ;
  usr::ArgumentExtender arg;
  arg.AddOptions( desc );
  arg.ParseOptions( argc, argv );

  if( !arg.CheckArg( "data" ) ){ return 0; }
  if( !arg.CheckArg( "output" ) ){ return 0; }
  if( !arg.CheckArg( "min" ) ){ return 0; }
  if( !arg.CheckArg( "max" ) ){ return 0; }
  if( !arg.CheckArg( "nbins" ) ){ return 0; }

  std::cout << "Starting" << std::endl;

  std::fstream infile( arg.Arg<std::string>( "data" ), std::ios::in );
  std::string line;

  const double xmin    = arg.Arg<double>( "min" );
  const double xmax    = arg.Arg<double>( "max" );
  const unsigned nbins = arg.Arg<unsigned>( "nbins" );
  const double binwidth = (xmax - xmin )/nbins;

  TH1D hist( "data", "", nbins, xmin, xmax );
  TF1 f( "f", "gaus", xmin, xmax );

  // Getting ride of the first four line;
  std::getline( infile, line );
  std::getline( infile, line );
  std::getline( infile, line );
  std::getline( infile, line );

  while( std::getline( infile, line ) ){
    std::istringstream iss( line );
    double input;
    iss >> input;// Getting first token only
    if( input < xmin || input > xmax ){ continue; }
    hist.Fill( input );
  }

  TSpectrum spec(5); // 5 Peaks is enough.
  spec.Search( &hist, 1, "nobackground" );
  // spec.SetMarkerColorAlpha( usr::plt::col::white, 0 );

  const double x0  = spec.GetPositionX()[0];
  const double sig = 3*binwidth ;
  f.SetRange( x0 - 1.5*sig, x0 + 1.5*sig );
  auto fit = hist.Fit( &f, "QN0 L S R" );

  usr::plt::Simple1DCanvas c;

  c.PlotHist( hist,
    usr::plt::PlotType( usr::plt::scatter ) );

  c.PlotFunc( f,
    usr::plt::EntryText( "Gaussian Fit" ),
    usr::plt::VisualizeError( fit ),
    usr::plt::PlotType( usr::plt::fittedfunc ),
    usr::plt::LineColor( usr::plt::col::blue ),
    usr::plt::FillColor( usr::plt::col::cyan )
    );

  c.PlotHist( hist,
    usr::plt::EntryText( "Data" ),
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
    usr::plt::MarkerSize( 0.5 ),
    usr::plt::MarkerColor( usr::plt::col::black ),
    usr::plt::LineColor( usr::plt::col::black )
    );
  c.Xaxis().SetTitle( "Time measured [ns]" );
  c.Yaxis().SetTitle( "Events" );
  c.Xaxis().SetNdivisions( 505 );

  const std::string toplabel =
    arg.Arg( "data" ).find( "jitter" ) != std::string::npos ? "Jitter test" :
    arg.Arg( "data" ).find( "Laser" ) != std::string::npos ? "Laser setup" :
    arg.Arg( "data" ).find( "Pulse" ) != std::string::npos ? "LED Pulser" :
    "";

  c.DrawLuminosity( toplabel );
  c.DrawCMSLabel( "", "Delay width" );
  c.Pad().WriteLine( usr::fstr( "Gauss width=%.0lf#pm%.0lf [ps]",
    f.GetParameter( 2 ) * 1000, f.GetParError( 2 ) * 1000 ) );
  std::cout << f.GetParameter( 0 ) << " "
            << f.GetParameter( 1 ) << " "
            << f.GetParameter( 2 ) << std::endl;

  c.SaveAsPDF( arg.Arg( "output" ) );

  return 0;
}
