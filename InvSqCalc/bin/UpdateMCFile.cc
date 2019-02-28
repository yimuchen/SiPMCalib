#include "SiPMCalib/InvSqCalc/interface/MCFormat.hpp"

#include "UserUtils/Common/interface/ArgumentExtender.hpp"

#include <iostream>
#include <boost/format.hpp>

int main( int argc, char** argv )
{
  usr::po::options_description desc("Update MC Files options");
  desc.add_options()
    ("file,f", usr::po::value<std::string>(), "File for input/output" )
    ("r", usr::po::value<std::vector<double>>()->multitoken(), "Aperture radius list" )
    ("l", usr::po::value<std::vector<double>>()->multitoken(), "SiPM size")
    ("o", usr::po::value<std::vector<double>>()->multitoken(), "offset")
    ("z", usr::po::value<std::vector<double>>()->multitoken(), "z value")
    ("n", usr::po::value<int>(), "number of runs to perform")
  ;

  usr::ArgumentExtender arg;
  arg.AddOptions( desc );

  arg.ParseOptions( argc, argv ); 

  MCManager* mgr =  usr::fs::exists(arg.Arg("file")) ? 
                    new MCManager( arg.Arg("file") ) :
                    new MCManager() ; 

  for( const auto rv : arg.ArgList<double>("r") ){
    for( const auto lv : arg.ArgList<double>("l") ){
      for( const auto ov : arg.ArgList<double>("o") ){
        for( const auto zv : arg.ArgList<double>("z") ){
          std::cout << boost::format("%lf %lf %lf %lf")
            % rv % lv %ov %zv  << std::endl;
          mgr->GetModel( rv, lv, ov, zv ).Run( arg.Arg<int>("n") );
        }
      }
    }
  }

  mgr->SaveToTXT( arg.Arg("file") ) ; 
  delete mgr ; 

  return 0 ; 
}