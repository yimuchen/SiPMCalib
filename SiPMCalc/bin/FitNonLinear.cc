#include "SiPMCalib/SiPMCalc/interface/SiPMNonLinearFit.hpp"
#include "UserUtils/Common/interface/ArgumentExtender.hpp"

// Helper functions

int
main( int argc, char*argv[] )
{
  usr::po::options_description desc(
    "Additional input output options for the SiPM analysis chain" );
  desc.add_options()
    ( "zscandata",
    usr::po::reqvalue<std::string>(),
    "Data file to the luminosity zscan data" )
    ( "refpoint",
    usr::po::multivalue<double>(),
    "Reference point to use as the number of photons. This option expects at least 3 doubles"
    "numbers to specify the referene point"
    "- The number of estimated photons at the point of reference"
    "- The z coordinate value of the estimate photons at the point of reference [mm]"
    "- The estimated bias voltage [mV]" )
    ( "gain",
    usr::po::defvalue<double>( 1.0 ),
    "Whether or not to specify the gain of the the SiPM" )
    ( "ped",
    usr::po::defvalue<double>( 0.0 ),
    "Whether or not to specify the pedestal, this is automatically evaluated if not specified." )
    ( "linplot",
    usr::po::value<std::string>(),
    "Output PDF file to place linearity fit result" )
  ;

  usr::ArgumentExtender args;
  args.AddVerboseOpt();
  args.AddOptions( SiPMNonLinearFit::FitArguments() );
  args.AddOptions( SiPMNonLinearFit::OutputArguments() );
  args.AddOptions( desc );
  args.ParseOptions( argc, argv );

  std::unique_ptr<SiPMNonLinearFit> fitter( new SiPMNonLinearFit() );

  // Getting the reference point
  const auto refargs = args.ArgList<double>( "refpoint" );
  assert( refargs.size() >= 3 );
  const double ref_nphotons = refargs.at( 0 );
  const double ref_z        = refargs.at( 1 );
  const double ref_bias     = refargs.at( 2 );
  const double gain         = args.Arg<double>( "gain" );
  const double ped          = args.Arg<double>( "ped" );

  // Updating the arguments
  usr::log::PrintLog( usr::log::DEBUG, "Updating settings" );
  fitter->UpdateSettings( args );

  // Running the main control flow
  usr::log::PrintLog( usr::log::DEBUG, "Making the fast component" );
  fitter->ReadFiles( args.Arg<std::string>( "zscandata" ) );
  fitter->MakeLinearGraph();

  usr::log::PrintLog( usr::log::DEBUG, "Making the nonlinear graph" );
  fitter->MakeNonLinearGraph( ref_nphotons, ref_z, ref_bias, gain, ped );

  // Running the fit
  usr::log::PrintLog( usr::log::DEBUG, "Running the non-linear fit" );
  fitter->RunNonLinearFit();


  usr::log::PrintLog( usr::log::DEBUG, "Making the plot" );
  fitter->PlotNonLinearity( "mytest.pdf" );

  // Running additional objects

  return 0;
}
