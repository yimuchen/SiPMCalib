#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/PlotUtils/interface/PlotCommon.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

#include <fstream>
#include <sstream>
#include <string>

#include "TF1.h"
#include "TH1D.h"

static TH1D        MakeHistogram( const usr::ArgumentExtender& args );
static TF1         MakeFunction( const usr::ArgumentExtender& args );
static TFitResult* RunFit( TH1D&                        hist,
                           TF1&                         func,
                           const usr::ArgumentExtender& args  );
static void MakePlot( TH1D&                        hist,
                      TF1&                         func,
                      TFitResult*                  fit,
                      const usr::ArgumentExtender& args );

int
main( int argc, char*argv[] )
{
  usr::po::options_description iodesc( "Input/Output options" );
  iodesc.add_options()
    ( "data", usr::po::reqvalue<std::string>(), "Input data .txt file" )
    ( "output", usr::po::reqvalue<std::string>(), "Primary output file name" )
  ;
  usr::po::options_description fitdesc( "Fitting options" );
  fitdesc.add_options()
    ( "min", usr::po::reqvalue<double>(), "X axis minimum" )
    ( "max", usr::po::reqvalue<double>(), "X axis maximum value" )
    ( "nbins",
    usr::po::reqvalue<unsigned>(),
    "number of bins to split the x axis" )
    ( "fittype",
    usr::po::defvalue<std::string>( "gaus" ),
    "Type of function for fitting, must be a string representing a function name"
    " define in ROOT (gaus/expo... and such)" )
    ( "fitmin",
    usr::po::defvalue<double>( -1 ),
    "Fit range minimum, leave negative for automatic determination" )
    ( "fitmax",
    usr::po::defvalue<double>( -1 ),
    "Fit range maximum, leave negative for automatic determination" )
  ;
  usr::po::options_description plotdesc( "Plot options" );
  plotdesc.add_options()
    ( "plotlabel",
    usr::po::defvalue<std::string>( "" ),
    "Additional label to display on the plot" )
  ;

  usr::ArgumentExtender args;
  args.AddOptions( iodesc );
  args.AddOptions( fitdesc );
  args.AddOptions( plotdesc );
  args.ParseOptions( argc, argv );

  TH1D       hist = MakeHistogram( args );
  TF1        func = MakeFunction( args );
  TFitResult*fit  = RunFit( hist, func, args );
  MakePlot( hist, func, fit, args );

  return 0;
}


static TH1D
MakeHistogram( const usr::ArgumentExtender& args )
{
  std::fstream   infile( args.Arg<std::string>( "data" ), std::ios::in );
  const unsigned nbins = args.Arg<unsigned>( "nbins" );
  const double   xmin  = args.Arg<double>( "xmin" );
  const double   xmax  = args.Arg<double>( "xmax" );
  TH1D           hist( "data", "", nbins, xmin, xmax );
  std::string    line;

  // Getting ride of the first four line;
  std::getline( infile, line );
  std::getline( infile, line );
  std::getline( infile, line );
  std::getline( infile, line );

  while( std::getline( infile, line ) ){
    std::istringstream iss( line );
    double             input;
    iss >> input;// Getting first token only
    if( input < xmin || input > xmax ){ continue; }
    hist.Fill( input );
  }

  return hist;
}


static TF1
MakeFunction( const usr::ArgumentExtender& args )
{
  const std::string fittype = args.Arg<std::string>( "fittype" );
  const double      xmin    = args.Arg<double>( "xmin" );
  const double      xmax    = args.Arg<double>( "xmax" );

  if( fittype != "expo" && fittype != "gaus" ){
    usr::log::PrintLog( usr::log::WARNING,
                        "Function is not officially supported, program may misbehave" );
  }

  TF1 f( "f", fittype.c_str(), xmin, xmax );
  return f;
}


static TFitResult*
RunFit( TH1D& hist, TF1& func, const usr::ArgumentExtender& args  )
{
  // Additional parsing of range
  const double   xmin     = args.Arg<double>( "xmin" );
  const double   xmax     = args.Arg<double>( "xmax" );
  const unsigned nbins    = args.Arg<unsigned>( "nbins" );
  const double   binwidth = ( xmax-xmin ) / nbins;

  // Getting the central position
  const double x0 = hist.GetBinCenter( hist.GetMaximumBin() );

  // Additional parsing for fit range.
  const std::string ftype = args.Arg<std::string>( "fittype" );
  const double      fmin  = args.Arg<double>( "fitmin" );
  const double      fmax  = args.Arg<double>( "fitmax" );

  // Getting the actual fit range
  const double fitmin = fmin > 0 ?
                        fmin :
                        ftype == "expo" ?
                        x0+3 * binwidth :
                        x0-6 * binwidth;
  const double fitmax = fmax > 0 ?
                        fmax :
                        ftype == "expo" ?
                        x0+( nbins / 2 ) * binwidth :
                        x0+6 * binwidth;
  func.SetRange( fitmin, fitmax );
  return hist.Fit( &func, "QN0 L S R" ).Get();
}


static void
MakePlot( TH1D&                        hist,
          TF1&                         func,
          TFitResult*                  fit,
          const usr::ArgumentExtender& args )
{
  // Additional parsing
  const std::string ftype    = args.Arg<std::string>( "fittype" );
  const std::string fitlabel = ftype == "expo" ?
                               "Exponential Fit" :
                               ftype == "gaus" ?
                               "Gaussian Fit" :
                               "Custom fit";
  const std::string label = args.Arg<std::string>( "plotlabel" );

  usr::plt::Simple1DCanvas c;

  c.PlotHist( hist, usr::plt::PlotType( usr::plt::scatter ) );

  c.PlotFunc( func,
              usr::plt::EntryText( fitlabel ),
              usr::plt::VisualizeError( fit ),
              usr::plt::PlotType(
                usr::plt::fittedfunc ),
              usr::plt::LineColor(
                usr::plt::col::blue ),
              usr::plt::FillColor( usr::plt::col::cyan ));

  c.PlotHist( hist,
              usr::plt::EntryText( "Data" ),
              usr::plt::PlotType( usr::plt::scatter ),
              usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
              usr::plt::MarkerSize(
                0.3 ),
              usr::plt::MarkerColor( usr::plt::col::black ),
              usr::plt::LineColor( usr::plt::col::black ));
  c.Xaxis().SetTitle( "Time measured [ns]" );
  c.Yaxis().SetTitle( "Events" );
  c.Xaxis().SetNdivisions( 505 );

  c.DrawLuminosity( label );
  c.DrawCMSLabel( "", "Delay width" );
  if( ftype == "expo" ){
    const double tau = -1 / func.GetParameter( 1 );
    const double err = fabs( tau * func.GetParError( 1 ) / func.GetParameter(
                               1 ) );
    if( tau < 1 ){
      c.Pad().WriteLine( usr::fstr( "Pulse Tail=%.0lf#pm%.0lf ps",
                                    tau * 1000,
                                    err * 1000 ) );
    } else {
      c.Pad().WriteLine( usr::fstr( "Pulse Tail=%.2lf#pm%.2lf ns", tau, err ) );
    }
  } else {
    const double t = func.GetParameter( 2 );
    const double e = func.GetParError( 2 );
    if( t < 2 ){
      c.Pad().WriteLine( usr::fstr( "Pulse Width=%.0lf#pm%.0lf ps",
                                    t * 1000,
                                    e * 1000 ) );
    } else {
      c.Pad().WriteLine( usr::fstr( "Pulse Width=%.2lf#pm%.2lf ns", t, e ) );
    }
  }
  c.SaveAsPDF( args.Arg( "output" ) );
}
