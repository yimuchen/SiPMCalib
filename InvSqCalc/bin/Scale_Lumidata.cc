#include "UserUtils/Common/interface/ArgumentExtender.hpp"

#include <fstream>
#include <iostream>
#include <sstream>

int
main( int argc, char* argv[] )
{
  usr::po::options_description desc( "Options for shifting Lumi data" );
  desc.add_options()
    ( "data", usr::po::value<std::string>(),
    "Input data .txt file" )
    ( "output", usr::po::value<std::string>(),
    "Primary output file name" )
    ( "scale,s", usr::po::defvalue<double>( 1 ),
    "Scaling the luminosity and uncertainty by common value" )
    ( "uncscale,u", usr::po::defvalue<double>( 1 ),
    "Scaling the uncertainty by a value" )
    ( "pedestal,p", usr::po::defvalue<double>( 0 ),
    "Shifting the luminosity value by a certain amount, "
    "will take place before scaling" )
  ;

  usr::ArgumentExtender arg;
  arg.AddOptions( desc );
  arg.ParseOptions( argc, argv );

  const std::string inputfile  = arg.Arg<std::string>( "data" );
  const std::string outputfile = arg.Arg<std::string>( "output" );
  const double scale           = arg.Arg<double>( "scale" );
  const double unc_scale       = arg.Arg<double>( "uncscale" );
  const double ped             = arg.Arg<double>( "pedestal" );

  std::ifstream in( inputfile );
  std::ofstream out( outputfile );

  std::string line;

  while( getline( in, line ) ){
    std::stringstream ss( line );
    double x, y, z, lumi, unc, capval;
    ss >> x >> y >> z >> lumi >> unc >> capval;

    // Scaling the data
    lumi += ped;
    lumi *= scale;
    unc  *= scale * unc_scale;

    out << usr::fstr( "%.1lf %.1lf %.1lf %lf %lf %lf\n"
                    , x, y, z, lumi, unc, capval ) << std::flush;
  }

  return 0;

}
