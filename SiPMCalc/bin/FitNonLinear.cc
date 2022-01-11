/**
 * @file FitNonLinear.cc
 * @ingroup SiPMCalc
 * @brief Script for running the Fit of data of a zscan result to the nonlinear
 * model.
 *
 * @details Basic scripts for running the fitting routine defined
 * SiPMNonLinearFit. Asside from the options diefined the SiPMNonLinearFit, the
 * scripts provides arguments to all the input arguments defined by the various
 * SiPMNonLinearFit function calls:
 *
 * - `zscandata`: The path to the data file collected via a z scaon command
 * - `corrdata`: The path to the LED bias and temperature compensation file.
 * - `refpoint`: The 3 inputs required for the "number of photons" reference
 *   point, the input expects 3 numbers:
 *   - The number of estimated photons are the specified reference.
 *   - The gantry z coordinates of the reference points in [mm]
 *   - The estimated bias voltage of the reference points
 * - The `gain` and `ped`estal subtraction required for the dat morphing.
 *
 * The various outputs plots:
 * - `linplot`: Where to store the plot for the inverse square law fit plots.
 * - `origplot`: Where to store the plot for all data points as a function of
 *   fitted z coordinates.
 * - `nlplot`: Where to store the plots for the nonlinear fit results plot. *
 *
 */

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
    ( "corrdata", usr::po::value<std::string>(), "Correction data to be used" )
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
    ( "origplot",
    usr::po::value<std::string>(),
    "Output PDF file to plot original data scatters" )
    ( "linplot",
    usr::po::value<std::string>(),
    "Output PDF file to place linearity fit result" )
    ( "nlplot",
    usr::po::reqvalue<std::string>(),
    "Output PDF file to place non-linearity fit result" )
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
  if( args.CheckArg( "corrdata" ) ){
    fitter->ReadFiles( args.Arg<std::string>( "zscandata" ),
                       args.Arg<std::string>( "corrdata" ));
  } else {
    fitter->ReadFiles( args.Arg<std::string>( "zscandata" ) );
  }

  // Running the linear fit
  fitter->MakeLinearGraph();
  usr::log::PrintLog( usr::log::DEBUG, "Making the nonlinear graph" );
  fitter->MakeNonLinearGraph( ref_nphotons, ref_z, ref_bias, gain, ped );

  // Plotting the original graph
  if( args.CheckArg( "origplot" ) ){
    fitter->PlotOriginal( args.Arg<std::string>( "origplot" ) );
  }

  // Running the fit
  usr::log::PrintLog( usr::log::DEBUG, "Running the non-linear fit" );
  fitter->RunNonLinearFit();

  usr::log::PrintLog( usr::log::DEBUG, "Making the plot" );
  fitter->PlotNonLinearity( args.Arg<std::string>( "nlplot" ) );

  // Running additional objects

  return 0;
}
