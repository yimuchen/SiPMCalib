#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"
#include "UserUtils/PlotUtils/interface/PlotCommon.hpp"
#include "UserUtils/Common/interface/ArgumentExtender.hpp"

#include <fstream>
#include <sstream>
#include <string>

#include "TH1D.h"
#include "TF1.h"

int main(int argc, char* argv[])
{
  usr::po::options_description desc( "Options for Pulser shape plotting" );
  desc.add_options()
    ( "data", usr::po::value<std::string>(),   "Input data .txt file" )
    ( "min", usr::po::value<double>(), "X axis minimum" )
    ( "max", usr::po::value<double>(), "X axis maximum" )
    ( "nbins", usr::po::value<unsigned>(), "number of bins" )
  ;
  usr::ArgumentExtender arg;
  arg.AddOptions( desc );
  arg.ParseOptions( argc, argv );

  if( !arg.CheckArg("data") ){ return 0 ; }
  if( !arg.CheckArg("min") ){ return 0 ; }
  if( !arg.CheckArg("max") ){ return 0 ; }
  if( !arg.CheckArg("nbins") ){ return 0 ; }

  std::fstream infile( arg.Arg<std::string>("data") , std::ios::in );
  std::string  line;

  const double xmin = arg.Arg<double>("min");
  const double xmax = arg.Arg<double>("max");
  const unsigned nbins = arg.Arg<unsigned>("nbins");

  TH1D hist("data", "", nbins, xmin, xmax );
  TF1  f("f", "expo", xmin, xmax );

  // Getting ride of the first four line;
  std::getline( infile, line );
  std::getline( infile, line );
  std::getline( infile, line );
  std::getline( infile, line );

  while( std::getline( infile, line ) ){
    std::istringstream iss(line);
    double input;
    iss >> input; // Getting first token only
    input /= 1000.;
    if( input < xmin || input > xmax ){ continue; }
    hist.Fill( input );
  }

  auto fit = hist.Fit( &f, "QN0 L S R" );

  usr::plt::Simple1DCanvas c;


  c.PlotFunc( f,
    usr::plt::EntryText( "Exponential Fit"),
    usr::plt::VisualizeError( fit ),
    usr::plt::PlotType( usr::plt::fittedfunc ),
    usr::plt::LineColor( usr::plt::col::blue ),
    usr::plt::FillColor( usr::plt::col::cyan )
    );

  c.PlotHist( hist,
    usr::plt::EntryText("Data"),
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
    usr::plt::MarkerSize( 0.5 ),
    usr::plt::MarkerColor( usr::plt::col::black ),
    usr::plt::LineColor( usr::plt::col::black )
     );
  c.Xaxis().SetTitle("Time measured [#mu s]" );
  c.Yaxis().SetTitle("Events");
  c.Xaxis().SetNdivisions( 205 );
  c.DrawLuminosity("Random Trigger");
  c.DrawCMSLabel("","Light Timing");
  c.Pad().WriteLine( usr::fstr("#tau_{Dark} = %.1f#pm%.2lf [#mu s]",
    1/f.GetParameter(1) * (-1.) ,
    1/f.GetParError(1)
  ));

  c.SaveAsPDF( "RandomTrigger.pdf" );

  return 0;
}

