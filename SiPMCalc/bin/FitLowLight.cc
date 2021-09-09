#include "SiPMCalib/SiPMCalc/interface/SiPMLowLightFit.hpp"
#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/Common/interface/STLUtils/OStreamUtils.hpp"

#include <fstream>

int
main( int argc, char* argv[] )
{
  usr::po::options_description desc(
    "Single run fit options for default low light fit" );
  desc.add_options()
    ( "configfile,c", usr::po::value<std::string>(),
    "Configuration file for overloading default settings in json format. "
    "Leave empty to load all-default fitting" )
    ( "outputdir,o", usr::po::defvalue<std::string>( "results/" ),
    "Output directory of the various objects" )
    ( "commonpostfix,p", usr::po::defvalue<std::string>( "" ),
    "Output prefix for the all files" )
  ;

  usr::po::options_description savedesc(
    "Options for saving the results, leave blank to ignore save, notice that the"
    "savefile will be augmented by the common file parsings, Notice that if all "
    "of savefit, savelatex and savetex are not set, then the full full will not "
    "be run " );
  savedesc.add_options()
    ( "savefit", usr::po::value<std::string>(),
    "Saving the grand spectrum fit results" )
    ( "savelatex", usr::po::value<std::string>(),
    "Saving the fit results as a latex table" )
    ( "savetxt", usr::po::value<std::string>(),
    "Saving the fit results as a raw .txt file" )
    ( "saveestpeak", usr::po::value<std::string>(),
    "Saving the peak estimation plot" )
    ( "saveestgain", usr::po::value<std::string>(),
    "Saving the gain estimation fit plot" )
    ( "saveestwidth", usr::po::value<std::string>(),
    "Saving the width estimation fit plot" )
    ( "saveestpoisson", usr::po::value<std::string>(),
    "Saving the poisson estimation plot" )
  ;

  usr::po::options_description vdesc( "verbose level settings" );
  vdesc.add_options()
    ( "verbose", usr::po::defvalue<int>( usr::log::WARNING ), "Print level for routine" )
  ;

  usr::ArgumentExtender args;
  args.AddOptions( desc );
  args.AddOptions( savedesc );
  args.AddOptions( SiPMLowLightFit::DataArguments() );
  args.AddOptions( SiPMLowLightFit::FitArguments() );
  args.AddOptions( SiPMLowLightFit::OperationArguments() );
  args.AddOptions( SiPMLowLightFit::EstArguments() );
  args.AddOptions( vdesc );
  args.ParseOptions( argc, argv );

  usr::log::SetLogLevel( args.Arg<int>( "verbose" ) );
  args.AddDirScheme( usr::ArgumentExtender::ArgPathScheme( "outputdir", "" ) );
  args.AddNameScheme( usr::ArgumentExtender::ArgPathScheme( "commonpostfix", "" ) );
  SiPMLowLightFit* mgr = args.CheckArg( "configfile" ) ?
                         new SiPMLowLightFit( args.Arg<std::string>( "configfile" ) ) :
                         new SiPMLowLightFit();
  mgr->UpdateSettings( args );

  // Parsing the data
  usr::log::PrintLog( usr::log::INFO, "Parsing the data file" );
  mgr->MakeBinnedData();

  // Running the estimations
  usr::log::PrintLog( usr::log::INFO, "Running the PDF estimation" );
  mgr->RunPDFEstimation();

  // Saving the estimation plots if requested.
  if( args.CheckArg( "saveestpeak" ) ){
    mgr->PlotPeakFind( args.MakePDFFile( args.Arg<std::string>( "saveestpeak" ) ) );
  }
  if( args.CheckArg( "saveestgain" ) ){
    mgr->PlotGainFit( args.MakePDFFile( args.Arg<std::string>( "saveestgain" ) ) );
  }
  if( args.CheckArg( "saveestwidth" ) ){
    mgr->PlotWidthFit( args.MakePDFFile( args.Arg<std::string>( "saveestwidth" ) ) );
  }
  if( args.CheckArg( "saveestpoisson" ) ){
    mgr->PlotPoissonFit( args.MakePDFFile( args.Arg<std::string>( "saveestpoisson" ) ) );
  }

  const bool runfit = args.CheckArg( "savefit" ) || args.CheckArg( "savelatex" )
                      || args.CheckArg( "savetxt" );
  if( !runfit ){
    usr::log::PrintLog( usr::log::INFO, "Early exit!" );
    return 0;
  }
  // Running the fit
  usr::log::PrintLog( usr::log::INFO, "Running the PDF fitting routine" );
  mgr->RunFit();

  // Saving the requested outputs
  usr::log::PrintLog( usr::log::INFO, "Saving the fit result" );
  if( args.CheckArg( "savefit" ) ){
    mgr->PlotSpectrumFit( args.MakePDFFile( args.Arg<std::string>( "savefit" ) ) );
  }
  if( args.CheckArg( "savelatex" ) ){
    std::ofstream table( args.MakeTEXFile( args.Arg<std::string>( "savelatex" ) ) );
    mgr->PrintTable( table );
  }

  if( args.CheckArg( "savetxt" ) ){
    std::ofstream raw( args.MakeTXTFile( args.Arg<std::string>( "savetxt" ) ) );
    mgr->PrintRaw( raw );
  }

  delete mgr;
  return 0;
}
