#include "SiPMCalib/SiPMCalc/interface/ToyRunCommon.hpp"

#include <fstream>

int
main( int argc, char** argv )
{
  usr::ArgumentExtender arg( "data/modelcfg.json" );
  arg.AddOptions( ToyOptions() );
  arg.ParseOptions( argc, argv );

  std::ifstream fin( pullfilename( arg ), std::ifstream::in );

  // Variables for pull result reading
  std::string line;
  std::string varname;
  double pullmean, pullmeanerr;
  double pullsig, pullsigerr;
  double npullmean, npullmeanerr;
  double npullsig, npullsigerr;

  while( std::getline( fin, line ) ){

    std::istringstream linestream( line );
    linestream  >> varname
    >> pullmean >> pullmeanerr
    >> pullsig >> pullsigerr
    >> npullmean >> npullmeanerr
    >> npullsig >> npullsigerr;

    // return 0;
  }
  return 0;
}
