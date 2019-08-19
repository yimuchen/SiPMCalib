#include "UserUtils/Common/interface/BoostUtils/PTreeUtils.hpp"
#include "UserUtils/Common/interface/Maths.hpp"
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

#include "cmath"
#include "TF1.h"
#include "TH1D.h"

int
main( int argc, char* argv[] )
{
  usr::pt::ptree dataj =
    usr::GetSubTree( usr::FromJsonFile( "data/viscalib_test.cfg" ),
      "FOV scan calibration" );

  const std::string maintag = "FOV scan calibration";
  const std::string calibz  = usr::fstr( "%.1f", 12.3 );

  std::map<int, double> calibx;
  std::map<int, double> caliby;

  for( int chipid = -1; chipid < 41; ++chipid ){
    const std::string query = usr::fstr( "%d/%s", chipid, calibz );
    auto rawvec             = dataj.get_child( usr::pt::ptree::path_type( query, '/' ) );
    std::vector<double> vec;

    for( auto x : rawvec ){
      vec.push_back( x.second.get_value<double>() );
    }

    calibx[ chipid + 1 ] = vec[0];
    caliby[ chipid + 1 ] = vec[1];
  }

  TH1D base( "base", "", 10, 25,   45 );
  TH1D xdiff( "xdiff", "", 20, 39.5, 40.5 );
  TH1D ydiff( "ydiff", "", 20, 29.5, 30.5 );
  std::vector<double> xvec;
  std::vector<double> yvec;


  for( int i = 0; i < 42; ++i ){
    if( i % 6 != 5 ){
      const double x = calibx[i+1] - calibx[i];
      const double y = caliby[i+1] - caliby[i];
      xdiff.Fill( std::sqrt( x*x+y*y ) );
      xvec.push_back( std::sqrt( x*x+y*y ) );
    }
    if( i / 6 < 6 ){
      const double x = calibx[i+6] - calibx[i];
      const double y = caliby[i+6] - caliby[i];
      ydiff.Fill( std::sqrt( x*x+y*y ) );
      yvec.push_back( std::sqrt( x*x+y*y ) );
    }
  }

  {
    usr::plt::Simple1DCanvas c;
    c.PlotHist( base );

    c.PlotHist( xdiff,
      usr::plt::EntryText(
        usr::fstr( "X Direction (STD:%.2lfmm, %.1lf%%)",
          usr::StdDev( xvec ),
          100.0*usr::StdDev( xvec )/usr::Mean( xvec ) ) ),
      usr::plt::PlotType( usr::plt::scatter ),
      usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
      usr::plt::MarkerColor( usr::plt::col::blue ),
      usr::plt::MarkerSize( 0.5 ),
      usr::plt::LineColor( usr::plt::col::blue ) );

    c.PlotHist( ydiff,
      usr::plt::EntryText(
        usr::fstr( "Y Direction (STD:%.2lfmm, %.1lf%%)",
          usr::StdDev( yvec ),
          100.0* usr::StdDev( yvec )/usr::Mean( yvec ) ) ),
      usr::plt::PlotType( usr::plt::scatter ),
      usr::plt::MarkerStyle( usr::plt::sty::mkrsquare ),
      usr::plt::MarkerColor( usr::plt::col::red ),
      usr::plt::MarkerSize( 0.5 ),
      usr::plt::LineColor( usr::plt::col::red ) );

    c.Pad().Xaxis().SetTitle( "Distance from neighbor [mm]" );
    c.Pad().Yaxis().SetTitle( "Test points" );

    c.DrawCMSLabel( "", "Visual Alignment" );
    c.DrawLuminosity( "Mock Tile" );

    c.SaveAsPDF( "test.pdf" );
  }

}
