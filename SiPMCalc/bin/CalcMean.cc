#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/Common/interface/Format.hpp"
#include "UserUtils/Common/interface/Maths.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

int
main( int argc, char* argv[] )
{
  usr::po::options_description desc( "Options for calculating mean of area" );
  desc.add_options()
    ( "data", usr::po::value<std::string>(),
    "Input data .txt file" )
    ( "start", usr::po::value<int>()->default_value(-1000),
    "starting time slice" )
    ( "end", usr::po::value<int>()->default_value(2147483647),
    "ending time slice" )
  ;
  usr::ArgumentExtender arg;
  arg.AddOptions( desc );
  arg.ParseOptions( argc, argv );

  if( !arg.CheckArg( "data" ) ){
    std::cerr << "Please include input data" << std::endl;
    return 1;
  }

  const std::string file( arg.Arg( "data" ) );
  unsigned time, ncapture, presample, postsample;
  double convfactor;

  std::string line;
  std::ifstream fin( file, std::ios::in );

  unsigned linecount = 0;

  // Getting first line
  std::getline( fin, line );
  std::istringstream linestream( line );
  linestream >> time >> ncapture
  >> presample >> postsample
  >> convfactor;


  std::vector<double> arealist;


  // Getting all other lines
  while( std::getline( fin, line ) ){
    ++linecount;

    if( line.find( ' ' ) == std::string::npos ){
      double area = 0;

      const int op_start = std::max( -((int)presample), arg.Arg<int>("start") );
      const int op_end = std::min( (int)postsample, arg.Arg<int>("end") );

      for( int i = op_start; i < op_end; ++i ){
        const int idx = i - op_start;
        const int8_t v1 = line[2*idx] >= 'a' ? 10 + line[2*idx] - 'a' :
                          line[2*idx] >= 'A' ? 10 + line[2*idx] - 'A' :
                          line[2*idx] - '0';
        const int8_t v2 = line[2*idx+1] >= 'a' ? 10 + line[2*idx+1] - 'a' :
                          line[2*idx+1] >= 'A' ? 10 + line[2*idx+1] - 'A' :
                          line[2*idx+1] - '0';
        const int8_t v = v1 << 4 | v2;
        area += (double)v  * (double)time * ( -1 ) * ( convfactor * 256.0 );
      }

      arealist.push_back( area );

    } else {
      std::istringstream stream( line );
      double area;

      std::string tmp = line.substr( 0, 12 );

      while( stream.good() ){
        stream >> area;
        area *= (double)time;
        area *= (double)-1;
        area *= convfactor;
        arealist.push_back( area );
      }
    }

  }

  std::cout << arg.Arg("data") << " "
            << usr::fmt::base::decimal( usr::Mean( arealist ), 3 ) << ", "
            << usr::fmt::base::decimal( usr::StdDev( arealist ), 4 )
            << std::endl;


  return 0;
}
