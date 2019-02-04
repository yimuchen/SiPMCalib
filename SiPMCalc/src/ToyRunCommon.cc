#include "SiPMCalib/SiPMCalc/interface/ToyRunCommon.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <cmath>

usr::po::options_description
ToyOptions()
{
  usr::po::options_description desc( "Toy MC options" );
  desc.add_options()
    ( "nEvents,n",   usr::po::value<int>(), "Number of events in a toy set" )
    ( "nToys",       usr::po::value<int>(), "Number of toy models to run" )
    ( "model,m",     usr::po::value<std::string>(), "Model to run MC studies" )
    ( "fit,f",       usr::po::value<std::string>(), "binned or unbinned" )
    ( "rate,r",      usr::po::value<float>(), "Average rate of geiger charges" )
    ( "xtalkrate,x", usr::po::value<float>(), "specify cross talk rate" )
    ( "darkrate,d",  usr::po::value<float>(), "Dark current rate" )
    ( "alpha,a",     usr::po::value<float>(), "after pulsing rate" )
  ;
  return desc;
}

double
Pedestal( const usr::ArgumentExtender& arg )
{ return 0; }

double
Gain( const usr::ArgumentExtender& arg )
{ return 200; }

double
ComNoise( const usr::ArgumentExtender& )
{ return 40; }

double
PixNoise( const usr::ArgumentExtender& )
{ return 15; }

double
Mean( const usr::ArgumentExtender& arg )
{ return arg.ArgOpt<float>( "rate", 4 ); }

double
Lambda( const usr::ArgumentExtender& arg )
{ return arg.ArgOpt<float>( "xtalkrate", 0.2 ); }

double
DCFrac( const usr::ArgumentExtender& arg )
{ return arg.ArgOpt<float>( "darkrate", 0.1 ); }

double
Alpha( const usr::ArgumentExtender& arg )
{ return arg.ArgOpt<float>( "alpha", 0.1 ); }

double
Beta( const usr::ArgumentExtender& arg )
{ return Gain( arg )*1.1; }

double
Xmin( const usr::ArgumentExtender& arg )
{ return Pedestal( arg ) - 1.5*ComNoise( arg ); }

double
Xmax( const usr::ArgumentExtender& arg )
{
  return arg.Arg( "model" ) == "dark" ?
         Pedestal( arg ) + Gain( arg )
         + std::sqrt( ComNoise( arg )*ComNoise( arg ) + PixNoise( arg ) * PixNoise( arg ) ) :
         Pedestal( arg ) + Gain( arg ) * std::max( Mean( arg ) + 1.5*sqrt( Mean( arg ) ), 5.0 );
}

std::fs::path
filename( const usr::ArgumentExtender& arg )
{
  std::string ans = ( boost::format( "%s_nEvt%d_%s_r%lg_x%lg_dc%lg_a%lg" )
                      % arg.Arg<std::string>( "model" )
                      % arg.Arg<int>( "nEvents" )
                      % arg.Arg<std::string>( "fit" )
                      % Mean()
                      % Lambda()
                      % DCFrac()
                      % Alpha()
                      ).str();
  boost::replace_all( ans, ".", "p" );
  return usr::resultpath( "SiPMCalib", "SiPMCalc" ) / ( ans + ".txt" );
}
