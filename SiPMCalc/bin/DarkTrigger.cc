#include "SiPMCalib/Common/interface/MakeRooData.hpp"
#include "SiPMCalib/Common/interface/WaveFormat.hpp"
#include "SiPMCalib/SiPMCalc/interface/CrossTalkPdf.hpp"
#include "SiPMCalib/SiPMCalc/interface/SiPMDarkPdf.hpp"

#include "UserUtils/Common/interface/Maths.hpp"
#include "UserUtils/Common/interface/STLUtils/OStreamUtils.hpp"
#include "UserUtils/MathUtils/interface/Measurement.hpp"
#include "UserUtils/PlotUtils/interface/Ratio1DCanvas.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

#include <algorithm>
#include <boost/format.hpp>
#include <fstream>
#include <iostream>
#include <map>
#include <set>

#include "RooAddPdf.h"
#include "RooChi2Var.h"
#include "RooDataHist.h"
#include "RooExponential.h"
#include "RooRealVar.h"
#include "TF1.h"
#include "TProfile.h"

usr::Measurement CalcCrossTalk( const std::string&,
                                const std::string& );
usr::Measurement CalcDecayTime( const std::string&,
                                const std::string& );
std::vector<usr::Measurement> CalcAP( const std::string&, const std::string& );


int
main( int argc, char* argv[] )
{
  const usr::Measurement crosstalk = CalcCrossTalk( argv[1], argv[2] );
  const usr::Measurement decaytime = CalcDecayTime( argv[1], argv[2] );
  const auto&            vap       = CalcAP(        argv[1], argv[2] );
  const usr::Measurement tap       = vap[0];
  const usr::Measurement tdc       = vap[1];
  const usr::Measurement approb    = vap[2];

  std::cout << "\n\n\n\n"
            << usr::separator( '-', 80 )
            << std::endl;
  boost::format printfmt( "%30s | %20s" );
  usr::fout( "%30s | %20s\n",
             "Cross talk prob. [%]"
             ,
             usr::fmt::decimal( 100.0 * crosstalk, 3 ) );
  usr::fout( "%30s | %20s\n"
             ,
             "Discharge time [ns]"
             ,
             usr::fmt::decimal( decaytime, 3 ) );
  usr::fout( "%30s | %20s\n"
             ,
             "Dark discharge rate [Hz]"
             ,
             usr::fmt::decimal( 1. / ( tdc * 1e-9 ), 3 ) );
  usr::fout( "%30s | %20s\n"
             ,
             "After pulsing time [ns]"
             ,
             usr::fmt::decimal( tap, 3 ) );
  usr::fout( "%30s | %20s\n"
             ,
             "After pulsing prob. [%]"
             ,
             usr::fmt::decimal( 100.0 * approb, 3 ) );
  ;

  return 0;
}


std::vector<usr::Measurement>
CalcAP( const std::string& input,
        const std::string& output )
{
  WaveFormat   wformat( input );
  const double tmin = wformat.Time();
  const double tmax = wformat.Time() * wformat.NSamples();

  RooRealVar x( "x", "Second peak time delay", tmin, tmax, "ns" );
  x.setBins( wformat.NSamples() / 4  );
  RooDataHist data( "data", "", RooArgSet( x ) );

  for( unsigned i = 0; i < wformat.NWaveforms(); ++i  ){
    int16_t        cv0 = 0, cv1 = 0, cv2 = 0;// value for caching peak position
    const unsigned nopeak    = wformat.NSamples()+1;
    unsigned       firstpeak = nopeak;

    for( unsigned j = 0; j < wformat.NSamples(); ++j ){
      cv2 = wformat.WaveformRaw( i )[ j ];

      if( firstpeak == nopeak && cv1 <= cv0 && cv1 <= cv2 && cv1 < -1 ){
        firstpeak = i;
      }
      if( firstpeak != nopeak
          && firstpeak < i
          && cv1 < cv0-2 && cv1 < cv2
          && cv1 < -1// Adding a minimal threshold
          ){
        x = ( i-firstpeak ) * wformat.Time();
        data.add( RooArgSet( x ) );
        break;
      }
      cv0 = cv1;
      cv1 = cv2;
    }
  }

  RooRealVar ndc( "ndc", "", 1, data.sumEntries() * 2 );
  RooRealVar nap( "nap", "", 1, data.sumEntries() * 2 );
  RooRealVar dcf( "dcf", "", -1e-2, -1e-30 );
  RooRealVar apf( "apf", "", -100, -1e-30 );

  RooExponential pdc( "pdc", "", x, dcf );
  RooExponential pap( "pap", "", x, apf );

  RooAddPdf p( "p", "",
               RooArgList( pdc, pap ),
               RooArgList( ndc, nap ) );

  dcf = -1e-6;
  apf = -1e-2;

// Fitting two distinct ranges
  pdc.fitTo( data, RooFit::Range( 100, tmax ) );
  pap.fitTo( data, RooFit::Range( 30, 50 ) );

  p.fitTo( data, RooFit::Range( 20, tmax ) );

  const usr::Measurement ap_c( -apf.getVal(), apf.getError(), apf.getError() );
  const usr::Measurement dc_c( -dcf.getVal(), dcf.getError(), dcf.getError() );
  const usr::Measurement tap = 1. / ( ap_c );
  const usr::Measurement tdc = 1. / ( dc_c );

  const usr::Measurement Nap( nap.getVal(), nap.getError(), nap.getError() );
  const usr::Measurement Ndc( ndc.getVal(), ndc.getError(), ndc.getError() );

  const double           tr      = 10;  // Is this number valid??
  const usr::Measurement Nap_int = Nap * ( tap ) * ( 1. / ( tap+tr ) );
  const usr::Measurement Ndc_int = Ndc * ( tdc ) * ( 1. / ( tdc+tr ) );
  const usr::Measurement prob    = Nap_int * ( 1. / ( Nap_int+Ndc_int ) );

  std::cout << "T_{ap}=" << usr::fmt::decimal( tap ) << " [ns]" << std::endl;
  std::cout << "T_{dc}=" << usr::fmt::decimal( tdc ) << " [ns]" << std::endl;
  std::cout << "N_{ap}=" << usr::fmt::decimal( Nap ) << " [events]" <<
    std::endl;
  std::cout << "N_{ap,int}=" << usr::fmt::decimal( Nap_int )
            << " [events]" << std::endl;
  std::cout << "N_{dc}=" << usr::fmt::decimal( Ndc ) << " [events]" <<
    std::endl;
  std::cout << "N_{dc,int}=" << usr::fmt::decimal( Ndc_int )
            << " [events]" << std::endl;
  std::cout << "P_ap  =" << usr::fmt::decimal( prob ) << std::endl;

  usr::plt::Ratio1DCanvas c( x );

  auto& dg = c.PlotData( data,
                         usr::plt::EntryText( "SiPM readout" ) );
  auto& dcg = c.PlotPdf( p,
                         RooFit::Components( pdc ),
                         usr::plt::EntryText( "Dark current" ) );
  auto& apg = c.PlotPdf( p,
                         RooFit::Components( pap ),
                         usr::plt::EntryText( "After pulse" ) );
  auto& pg = c.PlotPdf( p,
                        usr::plt::EntryText( "Fit" ) );

  dg.SetMarkerSize( 0.2 );
  dcg.SetLineColor( kGreen );
  apg.SetLineColor( kRed );
  pg.SetLineColor( kBlue );

  c.PlotScale( pg, pg );
  c.PlotScale( dg,
               pg,
               usr::plt::PlotType( usr::plt::scatter ) );

  c.DrawLuminosity( "Dark current trigger" );
  c.DrawCMSLabel( "", "Noise Parameters" );
  c.TopPad()
  .WriteLine( usr::fstr( "#tau_{AP} = %.1lf_{#pm%.2lf} [ns]",
                         tap.CentralValue(),
                         tap.AbsAvgError() ) )
  .WriteLine( usr::fstr( "#tau_{DC} = %.1lf_{#pm%.2lf} [#mu s]",
                         tdc.CentralValue() / 1000.,
                         tdc.AbsAvgError() / 1000. ) )
  .WriteLine( usr::fstr( "P_{AP} = %.1lf_{#pm%.2lf}[%%]",
                         100. * prob.CentralValue(),
                         100. * prob.AbsAvgError() ) );

  c.TopPad().SetLogy( 1 );
  c.TopPad().SetLogx( 1 );
  c.TopPad().SetYaxisMax( c.TopPad().GetYaxisMax() * 300 );
  c.BottomPad().SetLogx( 1 );
  c.BottomPad().Yaxis().SetTitle( "Data/Fit" );

  c.SaveAsPDF( output+"_ap_timeplot.pdf" );

  return {tap, tdc, prob};
}


usr::Measurement
CalcDecayTime( const std::string& input,
               const std::string& output )
{
  static const unsigned start = 3;
  static const unsigned end   = 50;

  WaveFormat wformat( input );

  RooRealVar  x( "x", "x", -10000, 20000 );
  RooRealVar  ped( "ped", "ped", -200, 200 );
  RooRealVar  gain( "gain", "gain", 1, 500 );
  RooRealVar  s0( "s0", "s0", 0.01, 100 );
  RooRealVar  s1( "s1", "s1", 0.01, 50 );
  RooRealVar  dcfrac( "dcfrac", "dcfrac", 0, 0.2 );
  RooRealVar  epsilon( "epslion", "epsilon", 1e-5, 1e-1 );
  SiPMDarkPdf pdf( "dark", "dark",
                   x, ped, gain, s0, s1, dcfrac, epsilon );

  // Making the data sets from the waveform data
  const auto list = wformat.SumList( start, end );
  SetRange( x, 2, -1, list );
  std::unique_ptr<RooDataHist> data( MakeData( x, list, -1 ) );
  pdf.RunEstimate( *data );

  const unsigned nbins = start+end;
  const double   xmin  = 0;
  const double   xmax  = wformat.Time() * nbins;

  TProfile p( usr::RandomString( 6 ).c_str(), "",
              nbins, xmin, xmax );
  p.SetErrorOption( "s" );
  p.SetStats( 0 );

  for( unsigned i = 0; i < wformat.NWaveforms(); ++i ){
    // Finding the peak
    int16_t      cv0 = 0, cv1 = 0, cv2 = 0;
    const double nopeak    = nbins+1;
    unsigned     localpeak = nopeak;

    // Finding the local peak position
    for( unsigned j = start; j < end; ++j ){
      const int16_t v = wformat.WaveformRaw( i )[j];
      if( localpeak == nopeak ){
        cv2 = v;
      }
      if( localpeak == nopeak && cv1 <= cv0 && cv1 <= cv2 && cv1 < -1 ){
        localpeak = i;
      }
      cv0 = cv1;
      cv1 = cv2;
    }

    const double area = wformat.WaveformSum( i, start, end );

    if( ped.getVal()-0.2 * gain.getVal()  < area
        && area < ped.getVal()+0.2 * gain.getVal() ){
      for( unsigned j = localpeak; j < localpeak+nbins; ++j ){
        const double x = wformat.Time() * (double)( j-localpeak );
        const double y = wformat.Waveform( i )[j];
        p.Fill( x, y );
      }
    }
  }

  TF1 f( usr::RandomString( 6 ).c_str(), "[0]*TMath::Exp(-x/[1]) + [2]",
         xmin, xmax * 0.8 );
  f.SetParameter( 0, p.GetBinContent( 0 ) );
  f.SetParameter( 1, 1 );
  p.Fit( &f, "QN0R" );// Chi squared fit

  usr::plt::Ratio1DCanvas c;

  auto& g = c.PlotFunc( f,
                        usr::plt::TrackY( usr::plt::tracky::both ),
                        usr::plt::EntryText( "Fit exp." ),
                        usr::plt::LineColor( usr::plt::col::red ) );
  auto& gp = c.PlotHist( p,
                         usr::plt::PlotType( usr::plt::scatter ),
                         usr::plt::EntryText( "1 Geiger Readout" ),
                         usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
                         usr::plt::MarkerSize( 0.2 ),
                         usr::plt::LineColor( usr::plt::col::darkgray ) );

  c.PlotScale( gp,
               g,
               usr::plt::PlotType( usr::plt::scatter ) );

  c.TopPad().DrawLuminosity( "Dark current trigger" );
  c.TopPad().DrawCMSLabel( "", "Noise Parameter" );
  c.TopPad().WriteLine( usr::fstr( "#tau = %.1lf_{#pm%.2lf}",
                                   f.GetParameter( 1 ),
                                   f.GetParError( 1 ) ) )
  .WriteLine( usr::fstr( "#chi^{2}/%d = %.2lf",
                         nbins-3,
                         f.GetChisquare() / ( nbins-3 ) ) );
  c.TopPad().Yaxis().SetTitle( "ADC value" );
  c.BottomPad().Xaxis().SetTitle( "Time [ns]" );
  c.BottomPad().Yaxis().SetTitle( "Data/Fit" );

  c.SaveAsPDF( output+"_decay_tracefit.pdf" );

  std::cout << f.GetParameter( 0 ) << " " << f.GetParError( 0 ) << std::endl
            << f.GetParameter( 1 ) << " " << f.GetParError( 1 ) << std::endl
            << f.GetParameter( 2 ) << " " << f.GetParError( 2 ) << std::endl;

  return usr::Measurement(
    f.GetParameter( 1 ),
    f.GetParError( 1 ),
    f.GetParError( 1 ) );
  return usr::Measurement( 0, 0, 0 );
}


usr::Measurement
CalcCrossTalk( const std::string& input,
               const std::string& output )
{
  static const unsigned start = 3;
  static const unsigned end   = 6;

  WaveFormat wformat( input );

  RooRealVar  x( "x", "x", -10000, 200000 );
  RooRealVar  ped( "ped", "ped", -200, 200 );
  RooRealVar  gain( "gain", "gain", 1, 100 );
  RooRealVar  s0( "s0", "s0", 0.01, 15 );
  RooRealVar  s1( "s1", "s1", 0.01, 5 );
  RooRealVar  dcfrac( "dcfrac", "dcfrac", 0, 0.2 );
  RooRealVar  epsilon( "epslion", "epsilon", 1e-5, 1e-1 );
  SiPMDarkPdf pdf( "dark", "dark",
                   x, ped, gain, s0, s1, dcfrac, epsilon );

  const auto list = wformat.SumList( start, end );
  SetRange( x, 2, -1, list );
  std::unique_ptr<RooDataHist> data( MakeData( x, list, -1 ) );

  pdf.RunEstimate( *data );

  std::vector<double> uniquearea = usr::RemoveDuplicate( list );
  TGraph              threshold( uniquearea.size()+1  );
  // Adding first point for aesthetic reasons.
  threshold.SetPoint( 0, list.front() / 2, list.size()+1 );

  for( unsigned i = 0, j = 0; i < uniquearea.size() && j < list.size(); ++i ){
    threshold.SetPoint( i+1, list.at( j ), list.size()-j );

    while( j < list.size() && list.at( j ) == uniquearea.at( i ) ){
      ++j;
    }
  }

  const unsigned           fitwidth = 7;
  std::map<double, double> derivmap;
  std::vector<double>      derivminraw;

  for( unsigned i = 0; i < uniquearea.size()-fitwidth; ++i ){
    const double xmin = uniquearea.at( i );
    const double xmax = uniquearea.at( i+fitwidth-1 );
    const double xcen = ( xmin+xmax ) / 2;

    TF1 f( usr::RandomString( 6 ).c_str(), "[0]*x + [1]", xmin, xmax );
    threshold.Fit( &f, "QN0R EX0" );

    derivmap[xcen] = f.GetParameter( 0 );
  }

  derivminraw.push_back( usr::GetMinimum( list ) * 0.75 );

  for( auto mapit = ++derivmap.begin(); mapit != derivmap.end(); ++mapit ){
    auto prev = mapit;
    --prev;
    auto next = mapit;
    ++next;
    if( next == derivmap.end() ){ break; }

    if( fabs( mapit->second ) < fabs( prev->second ) &&
        fabs( mapit->second ) < fabs( next->second ) ){
      derivminraw.push_back( mapit->first );
    }
  }

  std::vector<double> derivmin;

  // Removing neighbors
  for( unsigned i = 0; i < derivminraw.size();  ){
    std::vector<double> group;
    group.push_back( derivminraw.at( i ) );

    for( unsigned j = i+1; j < derivminraw.size(); ++j ){
      if( derivminraw.at( j )-group.front() < gain.getVal() / 2 ){
        group.push_back( derivminraw.at( j ) );
        if( j == derivminraw.size()-1 ){
          i = j;
          break;
        }
      } else {
        i = j;
        break;
      }
    }

    double mindiv  = 100000000000;
    double mindivx = 0;

    for( const double x : group ){
      if( derivmap[x] < mindiv ){
        mindiv  = derivmap[x];
        mindivx = x;
      }
    }

    derivmin.push_back( mindivx );
    if( i == derivminraw.size()-1 ){break;}
  }

  usr::plt::Simple1DCanvas c;

  c.PlotGraph( threshold,
               usr::plt::PlotType( usr::plt::simplefunc ) );

  c.Pad().SetTextAlign( usr::plt::font::bottom_right );

  for( unsigned i = 0; i < derivmin.size() && i < 3; ++i ){
    const double x = derivmin.at( i );
    c.Pad().DrawHLine(
      threshold.Eval( x ),
      usr::plt::LineColor( usr::plt::col::red ),
      usr::plt::LineStyle( usr::plt::sty::lindotted ) );
    c.Pad().WriteAtData(
      usr::GetMaximum( list ),
      threshold.Eval( x ) * 1.1,
      usr::fstr( "%.1f Threshold", i+0.5 ) );
  }

  const usr::Measurement p0515 = derivmin.size() >= 2 ? usr::Efficiency::Minos(
    threshold.Eval( derivmin.at( 1 ) ),
    threshold.Eval( derivmin.at( 0 ) ) ) : usr::Measurement( 0, 0, 0 );

  c.DrawLuminosity( "Dark current trigger" );
  c.DrawCMSLabel( "", "Noise Parameter" );
  c.Pad()
  .WriteLine( usr::fstr( "Int. Window: %d[ns]"
                         ,
                         ( end-start ) * wformat.Time() ) )
  .WriteLine( usr::fstr( "Cross-talk_{1.5/0.5} = %.2lf_{#pm%.3lf}%%"
                         ,
                         100 * p0515.CentralValue()
                         ,
                         100 * p0515.AbsAvgError() ) );

  c.Xaxis().SetTitle( "Area threshold [ADC #times ns]" );
  c.Yaxis().SetTitle( "Remaining events" );
  c.SetLogy( 1 );

  c.SaveAsPDF( output+"_xtalk_threshold.pdf" );

  return p0515;
  return usr::Measurement( 0, 0, 0 );
}
