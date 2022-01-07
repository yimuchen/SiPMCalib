#include "SiPMCalib/Common/interface/StdFormat.hpp"
#include "SiPMCalib/SiPMCalc/interface/NonLinearModel.hpp"

#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/MathUtils/interface/RootMathTools.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"


int
main( int argc, char* argv[] )
{
  usr::po::options_description desc( "Power correction generator" );
  desc.add_options()
    ( "powerdata",
    usr::po::reqvalue<std::string>(),
    "Data file for various power scans" )
    ( "corrfile",
    usr::po::reqvalue<std::string>(),
    "Output .txt file to store correction results." )
    ( "fitplot",
    usr::po::value<std::string>(),
    "Output .pdf file to store fit result plot." )
    ( "ped",
    usr::po::value<double>(),
    "Specifying the pedestal, if not specified the the constant component of the fit results will be taken as the pedestal" )
  ;

  usr::ArgumentExtender args;
  args.AddVerboseOpt();
  args.AddOptions( desc );
  args.ParseOptions( argc, argv );

  const StdFormat           input( args.Arg<std::string>( "powerdata" ) );
  const std::vector<double> bias    = input.Bias();
  const std::vector<double> readout = input.DataCol( 0 );
  const std::vector<double> unc     = input.DataCol( 1 );
  const std::vector<double> zero( 0, bias.size() );
  const double              bmin = usr::GetMinimum( bias );
  const double              bmax = usr::GetMaximum( bias );

  TGraphErrors g( bias.size(), bias.data(), readout.data(), zero.data(),
                  unc.data() );
  TF1 func( "func", &BiasModel, bmin, bmax, 2 );

  // Saving fit result
  TFitResult fit = usr::fit::FitGraph( g, func );

  // Saving the results to a corrector file.
  const double ped = args.CheckArg( "ped" ) ?
                     args.Arg<double>( "ped" ) :
                     func.GetParameter( 1 );
  std::ofstream outfile( args.Arg<std::string>( "corrfile" ) );
  for( const auto b : bias ){
    outfile << b << " " << func.Eval( b )-ped << std::endl;
  }

  if( args.CheckArg( "fitplot" )){
    // Saving the plot file
    usr::plt::Simple1DCanvas c;
    c.PlotFunc( func,
                usr::plt::VisualizeError( fit ),
                usr::plt::EntryText(
                  "Fit" ),
                usr::plt::TrackY( usr::plt::tracky::both ),
                usr::plt::PlotType( usr::plt::fittedfunc ),
                usr::plt::LineColor( usr::plt::col::blue ),
                usr::plt::FillColor( usr::plt::col::cyan ),
                usr::plt::MarkerSize( 0.2 ) );
    c.PlotGraph( g,
                 usr::plt::TrackY( usr::plt::tracky::both ),
                 usr::plt::PlotType(
                   usr::plt::scatter ),
                 usr::plt::MarkerSize( 0.2 ) );
    c.Xaxis().SetTitle( "Bias Voltage [mV]" );
    c.Yaxis().SetTitle( "Readout [mV-ns]" );
    c.Pad().SetLogy( kTRUE );

    c.SaveAsPDF( args.Arg<std::string>( "fitplot" ) );
  }
  return 0;
}
