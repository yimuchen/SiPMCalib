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
                        const std::string& baseline)
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
      area += ((double)v - _baseline[i]) * _timeint * ( -1 );
    }

    _arealist.push_back( area );
  }

  std::sort( _arealist.begin(), _arealist.end() );
  std::cout << "Done reading!" << std::endl;

  _binwidth = binwidth;
}

void
SiPMFormat::MakeDataSet( const double maxarea )
{
  const double xmin = usr::RoundDown( _arealist.front(), _binwidth );
  const double xmax = maxarea == 0 ?
                      usr::RoundUp( _arealist.back(), _binwidth )  :
                      usr::RoundUp( maxarea,          _binwidth );
  const double nbins = ( xmax - xmin ) / _binwidth;

  _x.reset( new RooRealVar(
    usr::RandomString( 6 ).c_str(), "Readout", xmin, xmax, "ADC #times ns" ) );
  x().setBins( nbins );
  _data.reset( new RooDataHist(
    usr::RandomString( 6 ).c_str(), "data", RooArgList( x() ) ) );

  for( const auto a : _arealist ){
    x() = a;
    if( a < xmax ){
      data().add( RooArgList( x() ) );
    }
  }
}

void SiPMFormat::MakeBaseLine( const std::string& file )
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
  linestream >> t >> n >>  pre >> post >> conv ;
  if( t    != TimeInterval() ||
      pre  != PreSamples()   ||
      post != PostSamples()  ||
      256 * conv != ADCConversion()  ){
        std::cout << "Bad format!" << std::endl;
        std::cout << t << " " << pre << " " << post << " " << conv << std::endl;
        std::cout << TimeInterval() << " " << PreSamples() << " " << PostSamples() << " " << ADCConversion() <<  std::endl;
        return ;
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
      _baseline[i] += (double)v;
      _baseline_err[i] += (double)v * (double) v;
    }
  }

  for( unsigned i = 0 ; i < _baseline.size() ; ++i ){
    _baseline[i]/= linecount;
    _baseline_err[i] = TMath::Sqrt(_baseline_err[i]/linecount
                                  - _baseline[i]*_baseline[i]);
    std::cout << i << " "
              << _baseline[i] << " "
              << _baseline_err[i] << std::endl;
  }
}

SiPMFormat::~SiPMFormat(){}

void
SiPMFormat::RunDarkEstimate( const std::string& savefit )
{
  const double xmin    = usr::RoundDown( _arealist.front(), _binwidth );
  const double xmax    = usr::RoundUp( _arealist.back(), _binwidth );
  const unsigned nbins = ( xmax - xmin ) / _binwidth;

  TH1D hist( usr::RandomString( 6 ).c_str(),
             usr::RandomString( 6 ).c_str(),
             nbins, xmin, xmax );
  TSpectrum spec( 5 );// 5 peaks should be plenty

  // Filling histograms
  for( const auto a : _arealist ){
    hist.Fill( a );
  }

  spec.Search( &hist, 1, "nobackground" );

  // Fitting first peak ;
  TF1 g0( usr::RandomString( 6 ).c_str(), "gaus",
          spec.GetPositionX()[0]  - 2*_binwidth,
          spec.GetPositionX()[0]  + 2*_binwidth );
  g0.SetParameter( 0, hist.Integral() );
  g0.SetParameter( 1, spec.GetPositionX()[0] );
  g0.SetParameter( 2, 2*_binwidth );
  hist.Fit( &g0, "QN0LR" );

  estped = g0.GetParameter( 1 );
  ests0  = g0.GetParameter( 2 );

  // Creating second spectrum for analysis
  const double xmin2  = usr::RoundUp( estped + 3*ests0, _binwidth );
  const double nbins2 = ( xmax - xmin2 ) / _binwidth;
  TH1D hist2( usr::RandomString( 6 ).c_str(), usr::RandomString( 6 ).c_str(),
              nbins2, xmin2, xmax );

  for( const auto a : _arealist ){
    if( a > xmin2 ){hist2.Fill( a );}
  }

  TSpectrum spec2( 5 );
  spec2.Search( &hist2, 2, "nobackground" );

  // Fitting using only the next peak
  TF1 g1( usr::RandomString( 6 ).c_str(), "gaus",
          spec2.GetPositionX()[0]  - 2*_binwidth,
          spec2.GetPositionX()[0]  + 2*_binwidth );

  hist2.Fit( &g1, "QLR" );

  ests1 = sqrt( g1.GetParameter( 2 )*g1.GetParameter( 2 ) - ests0*ests0 );
  // ests1 = 0;
  estgain   = g1.GetParameter( 1 ) - g0.GetParameter( 1 );
  estdcfrac = g1.GetParameter( 0 ) / g0.GetParameter( 0 );

  if( savefit == "" ){
    return;
  }

  {
    usr::plt::Simple1DCanvas c;
    c.PlotHist( hist, usr::plt::EntryText( "All data") );
    c.PlotHist( hist2, usr::plt::EntryText("Truncated data") );
    auto& graph0 = c.PlotFunc( g0, usr::plt::EntryText("Pedestal peak fit") );
    auto& graph1 = c.PlotFunc( g1, usr::plt::EntryText("1Geiger peak fit") );

    hist.SetLineColor( kBlack );
    hist2.SetLineColor( kBlue );
    graph0.SetLineColor( kRed );
    graph1.SetLineColor( kGreen );

    c.SaveAsPDF( savefit + "_peakfind.pdf" );
  }


}

void
SiPMFormat::RunLumiEstimate( const std::string& plotfit )
{
  const double xmin    = usr::RoundDown( _arealist.front(), _binwidth );
  const double xmax    = usr::RoundUp( _arealist.back(), _binwidth );
  const unsigned nbins = ( xmax - xmin ) / _binwidth;

  TH1D hist( usr::RandomString( 6 ).c_str(),
             usr::RandomString( 6 ).c_str(),
             nbins, xmin, xmax );
  TSpectrum spec( 20 );// 20 peaks should be plenty

  // Filling histograms
  for( const auto a : _arealist ){
    hist.Fill( a );
  }

  spec.Search( &hist, 1, "nobackground" );

  std::vector<double> peaklist;
  std::vector<TF1> fitlist;
  std::map<double, double> peakunc;
  std::map<double, double> widthmap;
  std::map<double, double> widthunc;
  std::map<double, double> heightmap;

  const int    maxbin  = hist.FindBin( spec.GetPositionX()[0] );
  const double histmax = hist.GetBinContent( maxbin );

  for( int i = 0; i < spec.GetNPeaks(); ++i ){
    const double x = spec.GetPositionX()[i];
    const int bin  = hist.FindBin( spec.GetPositionX()[i] );
    if( hist.GetBinContent( bin ) < ( histmax / 20. ) ){ continue; }

    fitlist.push_back( TF1( usr::RandomString( 6 ).c_str(), "gaus",
      x - 2*_binwidth, x + 2*_binwidth ) );
    auto& f = fitlist.back();
    f.SetParameter( 1, x );
    f.SetParameter( 2, 2*_binwidth );
    hist.Fit( &f, "QN0LR" );

    if( std::fabs(x-f.GetParameter(1)) > 3*_binwidth ){ continue; }

    peaklist.push_back( f.GetParameter( 1 ) );
    peakunc[   f.GetParameter( 1 ) ] = f.GetParError( 1 );
    widthmap[  f.GetParameter( 1 ) ] = f.GetParameter( 2 );
    widthunc[  f.GetParameter( 1 ) ] = f.GetParError( 2 );
    heightmap[ f.GetParameter( 1 ) ] = f.GetParameter( 0 );
  }

  const double primpeak = peaklist.front();
  std::sort( peaklist.begin(), peaklist.end() );
  const int primm = std::find( peaklist.begin(), peaklist.end(), primpeak )
                         - peaklist.begin();

  TGraphErrors gaingraph( peaklist.size() );
  TGraphErrors s1graph( peaklist.size() );
  TGraphErrors poigraph( peaklist.size() );

  for( unsigned i = 0; i < peaklist.size(); ++i ){
    const double peak = peaklist.at( i );
    gaingraph.SetPoint( i, i, peak - peaklist.at( 0 ) );
    gaingraph.SetPointError( i, 0, peakunc.at( peak ) );

    s1graph.SetPoint( i, i, widthmap.at( peak ) );
    s1graph.SetPointError( i, 0, widthunc.at( peak ) );

    poigraph.SetPoint( i, i, heightmap.at( peak ) );
    poigraph.SetPointError( i, 0, sqrt( heightmap.at( peak ) ) );
  }

  // Fitting function the peak positions
  auto fgain = []( const double* xx, const double* par ){
                 const double x = xx[0];
                 const double a = par[0];
                 return x*a;
               };

  auto fs1 = []( const double* xx, const double* par ){
               const double x  = xx[0];
               const double s0 = par[0];
               const double s1 = par[1];
               return sqrt( s0*s0 + x * s1 *s1 );
             };

  auto fpoi = []( const double* xx, const double* par ){
                const double x  = xx[0];
                const double N  = par[0];
                const double mu = par[1];
                const double l  = par[2];
                const double y  = ( mu + x * l );
                double prod     = 1;

                for( int i = 1; i <= x; ++i ){
                  prod *= y;
                  prod /= (double)( i );
                }

                prod *= mu / y;
                return N * prod * TMath::Exp( -y );
              };

  TF1 fg( usr::RandomString( 6 ).c_str(), fgain, 0, peaklist.size(), 1 );
  TF1 fs( usr::RandomString( 6 ).c_str(), fs1, 0, peaklist.size(), 2 );
  TF1 fp( usr::RandomString( 6 ).c_str(), fpoi, 0, peaklist.size(), 3 );

  // Initializing fitting requirements
  fg.SetParameter( 0, gaingraph.GetY()[1] );
  fs.SetParameter( 0, s1graph.GetY()[0] );
  fs.SetParameter( 1, 0 );
  fp.SetParameter( 0, hist.Integral() );
  fp.SetParameter( 1, std::max(0.99*primm,0.5) ); // Cannot be zero... breaks fit
  fp.SetParameter( 2, 0.1 );

  // Fitting graphs
  gaingraph.Fit( &fg, "QRN0 EX0" );
  s1graph.Fit( &fs, "QRN0 EX0" );
  poigraph.Fit( &fp, "QRN0 EX0" );

  estped    = peaklist.at( 0 );
  estgain   = fg.GetParameter( 0 );
  ests0     = fs.GetParameter( 0 );
  ests1     = fs.GetParameter( 1 );
  estmean   = fp.GetParameter( 1 );
  estlambda = fp.GetParameter( 2 );

  if( plotfit == "" ){
    return;
  }
  {
    usr::plt::Simple1DCanvas c;

    c.PlotHist( hist, usr::plt::EntryText( "SiPM readout" ) );

    for( auto& f : fitlist ){
      auto& g = &f == &( fitlist.front() ) ?
                c.PlotFunc( f,
        usr::plt::PlotType( usr::plt::simplefunc ),
        usr::plt::EntryText( "Local Gaussian fit" ) ) :
                c.PlotFunc( f, usr::plt::PlotType( usr::plt::simplefunc ) );
      g.SetLineColor( kRed );
    }

    for( int i = 0; i < spec.GetNPeaks(); ++i ){
      auto& line = c.Pad().DrawVLine(
          spec.GetPositionX()[i],
          usr::plt::LineColor( usr::plt::col::darkgray ),
          usr::plt::LineStyle( usr::plt::sty::lindotted ) ) ;
      if( i == 0 ){
        c.Pad().AddLegendEntry( line, "Peak search results", "L" );
      }
    }

    c.DrawCMSLabel( "Peak finding", "Spectral Fit" );
    c.SetHistAxisTitles( "Readout", "ADC #times ns" );
    c.SaveAsPDF( plotfit + "_peakfind.pdf" );
  }
  {
    usr::plt::Simple1DCanvas c;
    auto& g = c.PlotGraph( gaingraph,
      usr::plt::PlotType( usr::plt::scatter ),
      usr::plt::TrackY( usr::plt::TrackY::both ) );
    auto& f = c.PlotFunc( fg, usr::plt::PlotType( usr::plt::simplefunc ) );

    g.SetMarkerStyle( usr::plt::sty::mkrcircle );
    g.SetMarkerSize( 1 );

    f.SetLineColor( kRed );

    c.Xaxis().SetTitle( "Discharge peak" );
    c.Yaxis().SetTitle( "Fitted peak position [ADC #times ns]" );
    c.DrawCMSLabel( "Spectral Fit" );
    c.Pad().WriteLine( ( boost::format( "Est. gain=%.0lf#pm%.1f" )
                         % fg.GetParameter( 0 )
                         % fg.GetParError( 0 ) ).str() );

    c.SaveAsPDF( plotfit + "_gainfit.pdf" );
  }
  {
    usr::plt::Simple1DCanvas c;
    auto& g = c.PlotGraph( s1graph,
      usr::plt::PlotType( usr::plt::scatter ),
      usr::plt::TrackY( usr::plt::TrackY::both ) );
    auto& f = c.PlotFunc( fs, usr::plt::PlotType( usr::plt::simplefunc ) );

    g.SetMarkerStyle( usr::plt::sty::mkrcircle );
    g.SetMarkerSize( 1 );

    f.SetLineColor( kRed );

    c.Xaxis().SetTitle( "Discharge peak" );
    c.Yaxis().SetTitle( "Fitted peak width [ADC #times ns]" );
    c.DrawCMSLabel( "Spectral Fit" );
    c.Pad().WriteLine( ( boost::format( "Est. s_{0}=%.0lf#pm%.1f" )
                         % fs.GetParameter( 0 )
                         % fs.GetParError( 0 ) ).str() );
    c.Pad().WriteLine( ( boost::format( "Est. s_{1}=%.0lf#pm%.1f" )
                         % fs.GetParameter( 1 )
                         % fs.GetParError( 1 ) ).str() );

    c.SaveAsPDF( plotfit + "_widthfit.pdf" );
  }
  {
    usr::plt::Simple1DCanvas c;

    TH1D hist( usr::RandomString( 5 ).c_str(), "",
               peaklist.size(), -0.5, peaklist.size()-0.5 );

    for( unsigned i = 0; i < peaklist.size(); ++i ){
      hist.Fill( i, fp.Eval( i ) );
    }

    auto& g = c.PlotGraph( poigraph,
      usr::plt::PlotType( usr::plt::scatter ),
      usr::plt::TrackY( usr::plt::TrackY::both ) );
    auto& f = c.PlotHist( hist );

    g.SetMarkerStyle( usr::plt::sty::mkrcircle );
    g.SetMarkerSize( 1 );

    f.SetLineColor( kRed );

    c.Xaxis().SetTitle( "Discharge peak" );
    c.Yaxis().SetTitle( "Fitted peak height [Events]" );
    c.DrawCMSLabel( "Spectral Fit" );
    c.Pad().WriteLine( ( boost::format( "Est. mean=%.2lf#pm%.3f" )
                         % fp.GetParameter( 1 )
                         % fp.GetParError( 1 ) ).str() );
    c.Pad().WriteLine( ( boost::format( "Est. #lambda=%.3lf#pm%.4f" )
                         % fp.GetParameter( 2 )
                         % fp.GetParError( 2 ) ).str() );

    c.SaveAsPDF( plotfit + "_poissonfit.pdf" );
  }
}
