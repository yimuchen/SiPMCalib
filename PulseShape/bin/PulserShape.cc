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
    ( "fittype", usr::po::value<std::string>(), "Type of function for fitting" )
    ( "fitmin", usr::po::defvalue<double>( -1 ), "Fit range minimum" )
    ( "fitmax", usr::po::defvalue<double>( -1 ), "Fit range maximum" )
  ;
  usr::ArgumentExtender arg;
  arg.AddOptions( desc );
  arg.ParseOptions( argc, argv );

  const std::string data   = arg.Arg( "data" );
  const std::string output = arg.Arg( "output" );
  const double xmin        = arg.Arg<double>( "min" );
  const double xmax        = arg.Arg<double>( "max" );
  const unsigned nbins     = arg.Arg<unsigned>( "nbins" );
  const double binwidth    = ( xmax - xmin )/nbins;

  const std::string fittype
    = arg.CheckArg( "fittype" ) && arg.Arg( "fittype" ) == "expo" ? "expo" :
      "gaus";
  const std::string fitentry = fittype == "expo" ? "Exponential Fit" :
                               "Gaussian Fit";

  std::cout << "Starting" << std::endl;

  std::fstream infile( data, std::ios::in );
  std::string line;


  TH1D hist( "data", "", nbins, xmin, xmax );
  TF1 f( "f", fittype.c_str(), xmin, xmax );

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

  TSpectrum spec( 5 );// 5 Peaks is enough.
  spec.Search( &hist, 1, "nobackground nodraw goff" );
  // spec.SetMarkerColorAlpha( usr::plt::col::white, 0 );

  const double x0 = spec.GetPositionX()[0];
  const double fitmin
    = arg.Arg<double>( "fitmin" ) > 0 ? arg.Arg<double>( "fitmin" ) :
      fittype == "expo"               ? x0 + 3*binwidth :
      x0 - 6*binwidth;
  const double fitmax
    = arg.Arg<double>( "fitmax" ) > 0 ? arg.Arg<double>( "fitmax" ) :
      fittype == "expo"               ? x0 + ( nbins/2 )*binwidth :
      x0 + 6*binwidth;
  f.SetRange( fitmin, fitmax );
  auto fit = hist.Fit( &f, "QN0 L S R" );

  usr::plt::Simple1DCanvas c;

  c.PlotHist( hist,
    usr::plt::PlotType( usr::plt::scatter ) );

  c.PlotFunc( f,
    usr::plt::EntryText( fitentry ),
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
    data.find( "jitter" ) != std::string::npos ? "Jitter test" :
    data.find( "Laser"  ) != std::string::npos ? "Laser setup" :
    data.find( "LED"    ) != std::string::npos ? "Driven LED" :
    data.find( "Pulse"  ) != std::string::npos ? "LED Pulser" :
    "";

  c.DrawLuminosity( toplabel );
  c.DrawCMSLabel( "", "Delay width" );
  if( fittype == "expo" ){
    const double tau = -1/f.GetParameter( 1 );
    const double err = fabs( tau * f.GetParError( 1 ) / f.GetParameter( 1 ) );
    std::cout << tau << " " << err << std::endl;
    if( tau < 1 ){
      c.Pad().WriteLine( usr::fstr( "Pulse Tail=%.0lf#pm%.0lf ps",
        tau * 1000, err * 1000 ) );
    } else {
      c.Pad().WriteLine( usr::fstr( "Pulse Tail=%.2lf#pm%.2lf ns",
        tau, err ) );
    }
  } else {
    const double t = f.GetParameter( 2 );
    const double e = f.GetParError( 2 );
    if( t < 2 ){
      c.Pad().WriteLine( usr::fstr( "Pulse Width=%.0lf#pm%.0lf ps",
        t * 1000, e * 1000 ) );
    } else {
      c.Pad().WriteLine( usr::fstr( "Pulse Width=%.2lf#pm%.2lf ns",
        t, e ) );
    }
  }
  std::cout << f.GetParameter( 0 ) << " "
            << f.GetParameter( 1 ) << " "
            << f.GetParameter( 2 ) << std::endl;

  c.SaveAsPDF( arg.Arg( "output" ) );

  return 0;
}
