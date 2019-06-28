#include "SiPMCalib/InvSqCalc/interface/LEDFormat.hpp"

#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"

#include "TF2.h"
#include "TGraphErrors.h"
#include "TH2D.h"

#include <algorithm>
#include <fstream>
#include <sstream>

LEDManager::LEDManager( const usr::fs::path file )
{
  std::ifstream fin;
  std::string line;
  fin.open( file, std::ifstream::in );

  double x, y, z, lumi, lumierr;

  while( std::getline( fin, line ) ){
    std::istringstream linestream( line );
    linestream >> x >> y >> z >> lumi >> lumierr;

    _pointlist.push_back( LEDPoint() );
    _pointlist.back().x       = x;
    _pointlist.back().y       = y;
    _pointlist.back().z       = z;
    _pointlist.back().lumi    = lumi;
    _pointlist.back().lumierr = lumierr;
  }

  std::stable_sort( _pointlist.begin(), _pointlist.end(),
    []( const LEDPoint& x, const LEDPoint& y )->bool{
    return x.z < y.z;
  } );
}

void
LEDManager::FindCenter( const double z, double& x, double& y ) const
{
  TH2D* graph = MakeHScanGraph( z );
  TF2 func( "f", "[3]/((x-[0])*(x-[0]) + (y-[1])*(y-[1]) + [2]*[2])  + [4]" );

  auto fit = graph->Fit( &func, "Q N 0 EX0 S", "" );

  x = func.GetParameter( 0 );
  y = func.GetParameter( 1 );

  delete graph;
}

TH2D*
LEDManager::MakeHScanGraph( const double z ) const
{
  std::vector<double> xlist;
  std::vector<double> ylist;
  std::vector<double> zero;
  std::vector<double> lumi;
  std::vector<double> lumierr;

  for( const auto p : _pointlist ){
    if( p.z == z ){
      zero.push_back( 0 );
      xlist.push_back( p.x );
      ylist.push_back( p.y );
      lumi.push_back( p.lumi );
      lumierr.push_back( p.lumierr );
    }
  }

  const double xdiff = fabs( xlist.at( 0 ) - xlist.at( 1 ) );
  const double xmax  = *std::max_element( xlist.begin(), xlist.end() );
  const double xmin  = *std::min_element( xlist.begin(), xlist.end() );
  const double ymin  = *std::min_element( ylist.begin(), ylist.end() );
  const double ymax  = *std::max_element( ylist.begin(), ylist.end() );

  TH2D* ans = new TH2D( ( "hist"+usr::RandomString( 6 ) ).c_str(), "",
    ( xmax-xmin )/xdiff + 1, xmin-0.5*xdiff, xmax+0.5*xdiff,
    ( ymax-ymin )/xdiff + 1, ymin-0.5*xdiff, ymax+0.5*xdiff
     );

  for( unsigned i = 0; i < xlist.size(); ++i ){
    const int binidx = ans->FindBin( xlist.at( i ), ylist.at( i ) );
    ans->SetBinContent( binidx, lumi.at( i ) );
    ans->SetBinError( binidx, lumierr.at( i ) );
  }

  ans->SetStats( 0 );

  return ans;
}

TGraph*
LEDManager::MakeZScanGraph( const double x, const double y ) const
{
  std::vector<double> z;
  std::vector<double> zero;
  std::vector<double> lumi;
  std::vector<double> lumierr;

  for( const auto p : _pointlist ){
    if( p.x == x && p.y == y ){
      z.push_back( p.z );
      zero.push_back( 0 );
      lumi.push_back( p.lumi );
      lumierr.push_back( p.lumierr );
    }
  }

  return new TGraphErrors( z.size(),
    z.data(), lumi.data(),
    zero.data(), lumierr.data() );
}

TGraph*
LEDManager::MakeXScanGraph( const double y, const double z ) const
{
  std::vector<double> x;
  std::vector<double> zero;
  std::vector<double> lumi;
  std::vector<double> lumierr;

  for( const auto p : _pointlist ){
    if( p.y == y && p.z == z ){
      x.push_back( p.x );
      zero.push_back( 0 );
      lumi.push_back( p.lumi );
      lumierr.push_back( p.lumierr );
    }
  }

  return new TGraphErrors( x.size(),
    x.data(), lumi.data(),
    zero.data(), lumierr.data() );
}

double
LEDManager::Xmin() const
{
  return std::min_element( _pointlist.begin(), _pointlist.end(),
    []( const LEDPoint& x, const LEDPoint& y ){
    return x.x < y.x;
  }
    )->x;
}

double
LEDManager::Xmax() const
{
  return std::max_element( _pointlist.begin(), _pointlist.end(),
    []( const LEDPoint& x, const LEDPoint& y ){
    return x.x < y.x;
  }
    )->x;
}

double
LEDManager::Ymin() const
{
  return std::min_element( _pointlist.begin(), _pointlist.end(),
    []( const LEDPoint& x, const LEDPoint& y ){
    return x.y < y.y;
  }
    )->y;
}

double
LEDManager::Ymax() const
{
  return std::max_element( _pointlist.begin(), _pointlist.end(),
    []( const LEDPoint& x, const LEDPoint& y ){
    return x.y < y.y;
  }
    )->y;
}

double
LEDManager::Zmin() const
{
  return std::min_element( _pointlist.begin(), _pointlist.end(),
    []( const LEDPoint& x, const LEDPoint& y ){
    return x.z < y.z;
  }
    )->z;
}

double
LEDManager::Zmax() const
{
  return std::max_element( _pointlist.begin(), _pointlist.end(),
    []( const LEDPoint& a, const LEDPoint& b ){
    return a.z < b.z;
  }
    )->z;
}


double
LEDManager::LumiMin() const
{
  return std::min_element( _pointlist.begin(), _pointlist.end(),
    []( const LEDPoint& x, const LEDPoint& y ){
    return x.lumi < y.lumi;
  }
    )->lumi;
}

double
LEDManager::LumiMax() const
{
  return std::max_element( _pointlist.begin(), _pointlist.end(),
    []( const LEDPoint& x, const LEDPoint& y ){
    return x.lumi < y.lumi;
  }
    )->lumi;
}
