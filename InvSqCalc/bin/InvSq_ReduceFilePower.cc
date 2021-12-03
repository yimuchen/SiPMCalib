/**
 * @file InvSq_ReduceFilePower.cc
 * @brief Reducing an input file to contain just rows with power levels within a
 * certain range.
 */

#include "SiPMCalib/Common/interface/StdFormat.hpp"
#include "UserUtils/Common/interface/ArgumentExtender.hpp"

int
main( int argc, char** argv )
{
  usr::po::options_description desc(
    "Options for reducing the file via power output" );
  desc.add_options()
    ( "input,i", usr::po::reqvalue<std::string>(), "Input data file" )
    ( "output,o", usr::po::reqvalue<std::string>(), "Output data file" )
    ( "min", usr::po::reqvalue<double>(), "Minimum power value" )
    ( "max", usr::po::reqvalue<double>(), "Maximum power value" )
  ;

  usr::ArgumentExtender arg;
  arg.AddOptions( desc );
  arg.ParseOptions( argc, argv );

  const double pmin = arg.Arg<double>( "min" );
  const double pmax = arg.Arg<double>( "max" );

  StdFormat input( arg.Arg( "input" ) );
  auto      reduce = [pmin, pmax]( const StdFormat::RowFormat& x )->bool {
                       return x.bias > pmin && x.bias < pmax;
                     };
  StdFormat output = input.MakeReduced( reduce );
  output.WriteToFile( arg.Arg( "output" ) );
  return 0;
}
