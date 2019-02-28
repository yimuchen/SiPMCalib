#include "SiPMCalib/InvSqCalc/interface/MCFormat.hpp"

#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"
#include "UserUtils/Common/interface/SystemUtils/Time.hpp"
#include "UserUtils/MathUtils/interface/Measurement/CommonDistro.hpp"

#include "TGraphAsymmErrors.h"

#include <boost/format.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>

Model::Model( const double a, const double b, const double c, const double d ) :
  r( a ),
  l( b ),
  o( c ),
  z( d ),
  total( 0 ),
  passed( 0 )
{}

Model::~Model()
{}

double
Model::MaxTheta() const
{
  return atan( 1.5 * ( r + sqrt( l*l/2 ) + fabs( o ) ) / z );
}

double
Model::Sangle() const
{
  return 2* M_PI * ( 1-cos( MaxTheta() ) );
}

void
Model::Run( const unsigned n )
{
  std::uniform_real_distribution<double> pdf( 0, 1 );
  std::mt19937_64 gen( usr::CurrentTimeInNanSec() );

  for( unsigned i = 0; i < n; ++i ){
    const double genr = pdf( gen );
    const double gent = 2 * M_PI * pdf( gen );
    const double x0   = o + r * sqrt( genr ) * sin( gent );
    const double y0   = r * sqrt( genr ) * cos( gent );

    const double genz = 1 - ( 1-cos( MaxTheta() ) ) * pdf( gen );
    const double genp = 2 * M_PI * pdf( gen );
    const double genq = sqrt( 1-genz*genz );

    const double x1 = x0 + ( z/genz ) * genq * sin( genp );
    const double y1 = y0 + ( z/genz ) * genq * cos( genp );

    total++;
    if( fabs( x1 ) < l/2 && fabs( y1 ) < l/2 ){
      passed++;
    }
  }
}

usr::Measurement
Model::Lumi() const
{
  return M_PI * r * r * Sangle() * usr::Efficiency::Bayesian( passed, total );
}

/******************************/

MCManager::MCManager()
{}

MCManager::MCManager( const usr::fs::path& file )
{
  std::ifstream fin;
  std::string line;
  fin.open( file, std::ifstream::in );

  double r, l, o, z;
  unsigned total, passed;

  while( std::getline( fin, line ) ){
    std::istringstream linestream( line );
    linestream >> r >> l >> o >> z >> total >> passed;

    _modellist.push_back( Model( r, l, o, z ) );
    _modellist.back().total  = total;
    _modellist.back().passed = passed;
  }
}

MCManager::~MCManager()
{}


Model&
MCManager::GetModel(
  const double r, const double l, const double o, const double z )
{
  for( auto& m : _modellist ){
    if( m.r == r && m.l == l && m.o == o && m.z == z ){
      return m;
    }
  }

  _modellist.push_back( Model( r, l, o, z ) );

  while( _modellist.back().passed < 100 ){
    _modellist.back().Run( 100 );
  }

  return _modellist.back();
}

void
MCManager::SaveToTXT( const usr::fs::path& file )
{
  std::ofstream fout;
  fout.open( file );

  for( const auto& mod : _modellist ){
    fout << boost::format( "%8.1lf %8.1lf %8.1lf %8.1lf %8u %8u" )
      % mod.r % mod.l % mod.o % mod.z % mod.total % mod.passed
         <<  std::endl;
  }

  fout.close();
}


TGraphAsymmErrors
MCManager::MakeZScanGraph( const double r, const double l, const double o ) const
{
  std::vector<double> z;
  std::vector<double> zerr;
  std::vector<double> lumi;
  std::vector<double> lumihi;
  std::vector<double> lumilow;

  for( const auto& m : _modellist ){
    if( m.r == r && m.l == l && m.o == o && m.z >= 20 && m.z <= 450){
      const auto lumiM = m.Lumi();
      z.push_back( m.z );
      zerr.push_back( 0 );
      lumi.push_back( lumiM.CentralValue() );
      lumihi.push_back( lumiM.AbsUpperError() );
      lumilow.push_back( lumiM.AbsLowerError() );
    }
  }

  TGraphAsymmErrors ans(z.size());

  ans.SetName( usr::RandomString(6).c_str() );

  for( unsigned i = 0 ; i < z.size(); ++i ){
    ans.SetPoint(  i, z.at(i), lumi.at(i) );
    ans.SetPointEXhigh( i , 0 );
    ans.SetPointEXlow( i , 0 );
    ans.SetPointEYhigh( i, lumihi.at(i) );
    ans.SetPointEYlow( i, lumilow.at(i) );
  }

  return ans;
}

TGraphAsymmErrors
MCManager::MakeHScanGraph( const double r, const double l, const double z ) const
{
  std::vector<double> o;
  std::vector<double> oerr;
  std::vector<double> lumi;
  std::vector<double> lumihi;
  std::vector<double> lumilow;

  for( const auto& m : _modellist ){
    if( m.r == r && m.l == l && m.z == z ){
      const auto lumiM = m.Lumi();
      o.push_back( m.o );
      oerr.push_back( 0 );
      lumi.push_back( lumiM.CentralValue() );
      lumihi.push_back( lumiM.AbsUpperError() );
      lumilow.push_back( lumiM.AbsLowerError() );
    }
  }

  return TGraphAsymmErrors(
    o.size(),
    o.data(), lumi.data(),
    oerr.data(), oerr.data(),
    lumilow.data(), lumihi.data() );
}
