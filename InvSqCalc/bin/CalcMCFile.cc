#include "SiPMCalib/InvSqCalc/interface/MCFormat.hpp"

#include <boost/format.hpp>
#include <iostream>

int main( int argc, char** argv )
{
  MCManager mgr(argv[1]);

  for( const auto& m : mgr._modellist ){
    const auto lumi = m.Lumi(); 
    std::cout << 
      boost::format( "%8.1lf %8.1lf %8.1lf %8.1lf %10.6lf %10.6lf %10.6lf" )
      % m.r % m.l % m.o % m.z 
      % lumi.CentralValue() % lumi.AbsUpperError() % lumi.AbsLowerError() 
      << std::endl;
  }

  return 0 ;
}