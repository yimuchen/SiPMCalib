#include "SiPMCalib/SiPMCalc/interface/ToyRunCommon.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <cmath>

usr::po::options_description ToyOptions()
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
{ return 20; }

double
PixNoise( const usr::ArgumentExtender& )
{ return 5; }

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

usr::fs::path
filename( const usr::ArgumentExtender& arg )
{
  std::string ans = ( boost::format( "%s_nEvt%d_%s" )
                      % arg.Arg( "model" )
                      % arg.Arg<int>( "nEvents" )
                      % arg.Arg( "fit" ) ).str();

  if( arg.Arg("model") != "dark" ){
    ans += (boost::format("_r%lg_x%lg") % Mean(arg) % Lambda(arg)).str();
  }
  if( arg.Arg("model") != "ndc" && arg.Arg("model") != "simp" ){
    ans += (boost::format( "_dc%lg") % DCFrac(arg) ).str();
  }
  if( arg.Arg("model") != "nap" && arg.Arg("model") != "simp" ){
    ans += (boost::format("_a%lg") % Alpha(arg) ).str();
  }

  boost::replace_all( ans, ".", "p" );
  return usr::resultpath( "SiPMCalib", "SiPMCalc" ) / ( ans + ".txt" );
}

usr::fs::path
pullfilename( const usr::ArgumentExtender& arg )
{
  std::string ans = ( boost::format( "pullsummary_%s_nEvt%d_%s" )
                      % arg.Arg( "model" )
                      % arg.Arg<int>( "nEvents" )
                      % arg.Arg( "fit" ) ).str();

  if( arg.Arg("model") != "dark" ){
    ans += (boost::format("_r%lg_x%lg") % Mean(arg) % Lambda(arg)).str();
  }
  if( arg.Arg("model") != "ndc" && arg.Arg("model") != "simp" ){
    ans += (boost::format( "_dc%lg") % DCFrac(arg) ).str();
  }
  if( arg.Arg("model") != "nap" && arg.Arg("model") != "simp" ){
    ans += (boost::format("_a%lg") % Alpha(arg) ).str();
  }

  boost::replace_all( ans, ".", "p" );
  return usr::resultpath( "SiPMCalib", "SiPMCalc" ) / ( ans + ".txt" );
}

usr::fs::path
pullplotfilename( const usr::ArgumentExtender& arg, const std::string& name )
{
  std::string ans = ( boost::format( "%s_%s_nEvt%d_%s" )
                      % name
                      % arg.Arg( "model" )
                      % arg.Arg<int>( "nEvents" )
                      % arg.Arg( "fit" ) ).str();

  if( arg.Arg("model") != "dark" ){
    ans += (boost::format("_r%lg_x%lg") % Mean(arg) % Lambda(arg)).str();
  }
  if( arg.Arg("model") != "ndc" && arg.Arg("model") != "simp" ){
    ans += (boost::format( "_dc%lg") % DCFrac(arg) ).str();
  }
  if( arg.Arg("model") != "nap" && arg.Arg("model") != "simp" ){
    ans += (boost::format("_a%lg") % Alpha(arg) ).str();
  }

  boost::replace_all( ans, ".", "p" );
  return usr::resultpath( "SiPMCalib", "SiPMCalc" ) / "pullplot" / ( ans + ".pdf" );
}