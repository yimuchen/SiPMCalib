#include "SiPMCalib/SiPMCalc/interface/SiPMFormat.hpp"

#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/Maths.hpp"
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"
#else
#include "UserUtils/Common/Maths.hpp"
#include "UserUtils/Common/STLUtils/StringUtils.hpp"
#include "UserUtils/PlotUtils/Simple1DCanvas.hpp"
#endif

#include <algorithm>
#include <boost/format.hpp>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TSpectrum.h"

SiPMFormat::SiPMFormat( const std::string& file,
                        const double       binwidth,
                        const unsigned     start,
                        const unsigned     end,
                        const std::string& baseline )
{
  std::string line;
  std::ifstream fin( file, std::ios::in );


  unsigned linecount = 0;

  // Getting first line
  std::getline( fin, line );
  std::istringstream linestream( line );
  linestream >> _timeint >> _ncapture
  >> _presample >> _postsample
  >> _convfactor;

  MakeBaseLine( baseline );

  // Getting all other lines
  while( std::getline( fin, line ) ){
    ++linecount;
    if( linecount % 1000 == 0 ){
      std::cout << "\rLine " <<  linecount << "..." << std::flush;
    }
    double area = 0;

    for( unsigned i = start; i < end && i < _presample + _postsample; ++i ){
      const int8_t v1 = line[2*i] >= 'a' ? 10 + line[2*i] - 'a' :
                        line[2*i] >= 'A' ? 10 + line[2*i] - 'A' :
                        line[2*i] - '0';
      const int8_t v2 = line[2*i+1] >= 'a' ? 10 + line[2*i+1] - 'a' :
                        line[2*i+1] >= 'A' ? 10 + line[2*i+1] - 'A' :
                        line[2*i+1] - '0';
      const int8_t v = v1 << 4 | v2;
      area += ( (double)v - _baseline[i] ) * _timeint * ( -1 );
    }

    _arealist.push_back( area );
  }

  std::sort( _arealist.begin(), _arealist.end() );
  std::cout << "Done reading!" << std::endl;

  _binwidth = binwidth;

  const double xmin  = usr::RoundDown( _arealist.front(), _binwidth );
  const double xmax  = usr::RoundUp(   _arealist.back(),  _binwidth );
  const double nbins = ( xmax - xmin ) / _binwidth;

  _x.reset( new RooRealVar(
    usr::RandomString( 6 ).c_str(), "Readout", xmin, xmax, "ADC #times ns" ) );
  x().setBins( nbins );
  _data.reset( new RooDataHist(
    usr::RandomString( 6 ).c_str(), "data", RooArgList( x() ) ) );

  for( const auto a : _arealist ){
    x() = a;
    data().add( RooArgList( x() ) );
  }
}

void
SiPMFormat::TruncateDataSet( const double maxarea )
{
  const double xmax  = usr::RoundUp( maxarea,  _binwidth );
  const double nbins = (xmax - x().getMin() ) /_binwidth ;
  x().setMax( xmax );
  x().setBins( nbins );

  _data.reset( new RooDataHist(
    usr::RandomString( 6 ).c_str(), "data", RooArgList( x() ) ) );

  for( const auto a : _arealist ){
    if( a < maxarea ){
      x() = a;
      data().add( RooArgList( x() ) );
    }
  }
}

void
SiPMFormat::MakeBaseLine( const std::string& file )
{
  _baseline.clear();
  _baseline_err.clear();
  _baseline.resize( _presample+_postsample, 0 );
  _baseline_err.resize( _presample + _postsample, 0 );

  if( file == "" ){ return; }

  std::string line;
  std::ifstream fin( file, std::ios::in );

  unsigned t, n, pre, post;
  double conv;

  // Getting first line
  std::getline( fin, line );
  std::istringstream linestream( line );
  linestream >> t >> n >>  pre >> post >> conv;
  if( t    != TimeInterval() ||
      pre  != PreSamples()   ||
      post != PostSamples()  ||
      256 * conv != ADCConversion() ){
    std::cout << "Bad format!" << std::endl;
    std::cout << t << " " << pre << " " << post << " " << conv << std::endl;
    std::cout << TimeInterval() << " " << PreSamples() << " " << PostSamples() << " " << ADCConversion() <<  std::endl;
    return;
  }

  unsigned linecount = 0;

  // Getting all other lines
  while( std::getline( fin, line ) ){
    ++linecount;

    for( unsigned i = 0; i < _presample + _postsample; ++i ){
      const int8_t v1 = line[2*i] >= 'a' ? 10 + line[2*i] - 'a' :
                        line[2*i] >= 'A' ? 10 + line[2*i] - 'A' :
                        line[2*i] - '0';
      const int8_t v2 = line[2*i+1] >= 'a' ? 10 + line[2*i+1] - 'a' :
                        line[2*i+1] >= 'A' ? 10 + line[2*i+1] - 'A' :
                        line[2*i+1] - '0';
      const int8_t v = v1 << 4 | v2;
      _baseline[i]     += (double)v;
      _baseline_err[i] += (double)v * (double)v;
    }
  }

  for( unsigned i = 0; i < _baseline.size(); ++i ){
    _baseline[i]    /= linecount;
    _baseline_err[i] = TMath::Sqrt( _baseline_err[i]/linecount
      - _baseline[i]*_baseline[i] );
    std::cout << i << " "
              << _baseline[i] << " "
              << _baseline_err[i] << std::endl;
  }
}

SiPMFormat::~SiPMFormat(){}
