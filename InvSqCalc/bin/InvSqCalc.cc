#include "UserUtils/Common/interface/SystemUtils/Time.hpp"
#include "UserUtils/MathUtils/interface/Measurement.hpp"
#include "UserUtils/PlotUtils/interface/Constants.hpp"
#include "UserUtils/PlotUtils/interface/Ratio1DCanvas.hpp"

#include <boost/format.hpp>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "TF1.h"
#include "TGraphAsymmErrors.h"

static const double xmax = 1;

TGraph* ScanZ( const double   r,
               const double   offset,
               const unsigned nrun = 3e4 );
TGraph* ScanO( const double   r,
               const double   z,
               const unsigned nrun = 3e5  );
TGraph* MakeFitGraph( TGraph* );

TGraph* MakeAlignGraph( TGraph* );

usr::Measurement CalcLumi(
  const double   r,
  const double   offset,
  const double   z,
  const unsigned nrun
  );

int
main( int argc, char** argv )
{
  const std::vector<int> color = {
    kBlue, kRed, kGreen, kViolet, kOrange + 10
  };

  {// Making hole comparison diagram.
    std::vector<double> rlist = {0.5, 1, 2, 3};
    std::vector<TGraph*> rgraphlist;
    std::vector<TGraph*> fgraphlist;

    for( const auto r : rlist ){
      rgraphlist.push_back( ScanZ( r, 0 ) );
    }

    for( const auto g : rgraphlist ){
      fgraphlist.push_back( MakeFitGraph( g ) );
    }

    usr::plt::Ratio1DCanvas c;
    boost::format entry( "Aperture R = %lgmm" );

    for( unsigned i = 0; i < rlist.size(); ++i  ){
      c.PlotGraph( rgraphlist.at( i ),
        usr::plt::PlotType( usr::plt::scatter ),
        usr::plt::TrackY( usr::plt::TrackY::both ),
        usr::plt::EntryText( ( entry% rlist.at( i ) ).str() )
        );
      c.PlotGraph( fgraphlist.at( i ),
        usr::plt::PlotType( usr::plt::simplefunc ),
        usr::plt::TrackY( usr::plt::TrackY::none )
        );

      rgraphlist.at( i )->SetLineColor( color.at( i ) );
      rgraphlist.at( i )->SetMarkerColor( color.at( i ) );
      rgraphlist.at( i )->SetMarkerStyle( 20 );
      rgraphlist.at( i )->SetMarkerSize( 0.1 );
      fgraphlist.at( i )->SetLineColor( color.at( i ) );

      c.PlotScale(
        rgraphlist.at( i ), fgraphlist.at( i ),
        usr::plt::PlotType( usr::plt::scatter )
        );

      c.BottomPad().Xaxis().SetTitle( "Z distance [mm]" );
      c.BottomPad().Yaxis().SetTitle( "Exp./Fit" );
      c.TopPad().Yaxis().SetTitle( "Expected Luminosity (arbitrary units)" );
    }

    c.TopPad().DrawVLine( 100, kGray, 2 );
    c.TopPad().DrawCMSLabel( usr::plt::cap::sim, "HGCal" );
    c.TopPad().WriteLine( "1mm^{2} SiPM" );
    c.TopPad().WriteLine( "(aligned)" );

    c.TopPad().SetLogx( true );
    c.TopPad().SetLogy( true );
    c.BottomPad().SetLogx( true );
    c.SaveAsPDF( "HoleSize.pdf" );
  }

  for( const auto r : {1, 3, 5} ){// Making Offset comparison diagram
    std::vector<double> olist = {0,1,2,3,4,5};
    std::vector<TGraph*> ographlist;
    TGraph* fgraph;

    for( const auto o : olist ){
      ographlist.push_back( ScanZ( r, o, 3e5 ) );
    }

    fgraph = MakeFitGraph( ographlist.front() );

    usr::plt::Ratio1DCanvas c;
    boost::format entry( "Offset = %lgmm" );

    fgraph->SetLineColor( kBlack );
    c.PlotGraph( fgraph,
      usr::plt::PlotType( usr::plt::simplefunc ),
      usr::plt::EntryText( "Fit" )
      );
    c.PlotScale( fgraph, fgraph, usr::plt::PlotType( usr::plt::simplefunc ) );
    c.BottomPad().Xaxis().SetTitle( "Z distance [mm]" );
    c.BottomPad().Yaxis().SetTitle( "Exp./Fit" );
    c.TopPad().Yaxis().SetTitle( "Expected Luminosity (arbitrary units)" );

    c.TopPad().SetAxisFont();
    c.BottomPad().SetAxisFont();

    for( unsigned i = 0; i < olist.size(); ++i  ){
      c.PlotGraph( ographlist.at( i ),
        usr::plt::PlotType( usr::plt::scatter ),
        usr::plt::TrackY( usr::plt::TrackY::both ),
        usr::plt::EntryText( ( entry% olist.at( i ) ).str() )
        );

      ographlist.at( i )->SetLineColor( color.at( i%color.size() ) );
      ographlist.at( i )->SetMarkerColor( color.at( i%color.size() ) );
      ographlist.at( i )->SetMarkerStyle( 20 );
      ographlist.at( i )->SetMarkerSize( 0.1 );

      c.PlotScale( ographlist.at( i ), fgraph,
        usr::plt::PlotType( usr::plt::scatter )
        );

    }

    c.TopPad().DrawVLine( 100, kGray, 2 );
    c.TopPad().DrawCMSLabel( usr::plt::cap::sim, "HGCal" );
    c.TopPad().WriteLine( "1mm^{2} SiPM" );
    c.TopPad().WriteLine( ( boost::format( "r=%lgmm Aperture" )%r ).str() );

    c.TopPad().SetLogx( true );
    c.TopPad().SetLogy( true );
    c.BottomPad().SetLogx( true );
    c.SaveAsPDF( ( boost::format( "Offset_r%lf.pdf" )%r ).str() );
  }

  for( const auto r : {1, 3, 5} ){// Making plots for horizontal scan
    const std::vector<double> zlist = {5,10,15, 20};
    std::vector<TGraph*> ographlist;
    std::vector<TGraph*> fgraphlist;

    for( const auto z : zlist ){
      ographlist.push_back( ScanO( r, z, 3e6 ) );
      fgraphlist.push_back( MakeAlignGraph( ographlist.back() ) );
    }

    usr::plt::Ratio1DCanvas c;
    boost::format entry( "z=%lg" );

    for( unsigned i = 0; i < zlist.size(); ++i ){
      c.PlotGraph( ographlist.at( i ),
        usr::plt::PlotType( usr::plt::scatter ),
        usr::plt::TrackY( usr::plt::TrackY::both ),
        usr::plt::EntryText( ( entry%zlist.at( i ) ).str() ) );
      c.PlotGraph( fgraphlist.at( i ),
        usr::plt::PlotType( usr::plt::simplefunc )
        );

      ographlist.at( i )->SetLineColor( color.at( i ) );
      ographlist.at( i )->SetMarkerColor( color.at( i ) );
      ographlist.at( i )->SetMarkerStyle( 20 );
      ographlist.at( i )->SetMarkerSize( 0.1 );
      fgraphlist.at( i )->SetLineColor( color.at( i ) );

      c.PlotScale( ographlist.at( i ), fgraphlist.at( i ),
        usr::plt::PlotType( usr::plt::scatter ) );
    }

    c.BottomPad().Xaxis().SetTitle( "Offset [mm]" );
    c.BottomPad().Yaxis().SetTitle( "Exp./Fit" );
    c.TopPad().Yaxis().SetTitle( "Expected Luminosity (arbitrary units)" );
    c.TopPad().DrawCMSLabel( usr::plt::cap::sim, "HGCal" );
    c.TopPad().WriteLine( "1mm^{2} SiPM" );
    c.TopPad().WriteLine( ( boost::format( "R=%lgmm Aperture" )%r ).str() );

    c.TopPad().SetLogy( true );
    c.SaveAsPDF( ( boost::format( "Scan_r%lf.pdf" )%r ).str() );
  }


  return 0;
}


TGraph*
ScanZ( const double r, const double offset, const unsigned nrun )
{
  std::vector<double> vz;
  std::vector<double> vl;
  std::vector<double> ez;
  std::vector<double> elh;
  std::vector<double> ell;

  for( double z = 5; z <= 500; z += 15 ){
    const auto lumi = CalcLumi( r, offset, z, nrun );

    vz.push_back( z );
    vl.push_back( lumi.CentralValue() );
    ez.push_back( 0 );
    ell.push_back( lumi.AbsLowerError() );
    elh.push_back( lumi.AbsUpperError() );
  }

  std::cout << "... Done!" << std::endl;

  return new TGraphAsymmErrors( vz.size(),
    vz.data(), vl.data(),
    ez.data(), ez.data(),
    ell.data(), elh.data() );
}

TGraph*
ScanO( const double r, const double z, const unsigned nrun )
{
  std::vector<double> vo;
  std::vector<double> vl;
  std::vector<double> ez;
  std::vector<double> elh;
  std::vector<double> ell;

  for( double o = -10; o <= 10; o += 0.5  ){
    const auto lumi = CalcLumi( r, o, z, nrun );
    vo.push_back( o );
    vl.push_back( lumi.CentralValue() );
    ez.push_back( 0 );
    ell.push_back( lumi.AbsLowerError() );
    elh.push_back( lumi.AbsUpperError() );
  }

  std::cout << "... Done!" << std::endl;

  return new TGraphAsymmErrors(
    vo.size(),
    vo.data(), vl.data(),
    ez.data(), ez.data(),
    ell.data(), elh.data()
    );
}


TGraph*
MakeFitGraph( TGraph* r )
{
  TF1 f( "f", "[0]/(x*x)", 10, 500 );
  r->Fit( &f, "Q EX0 N", "", 100, 500 );
  std::vector<double> vx;
  std::vector<double> vy;

  for( double x = 5; x <= 500; x++ ){
    vx.push_back( x );
    vy.push_back( f.Eval( x ) );
  }

  return new TGraph( vx.size(), vx.data(), vy.data() );
}

TGraph*
MakeAlignGraph( TGraph* g )
{
  TF1 m( "g", "[0] + [1]*x*x", -5, +5 );
  g->Fit( &m, "Q EX0 N", "", -2, 2 );
  std::vector<double> vx;
  std::vector<double> vy;

  for( double x = -3; x < 3; x += 0.1 ){
    vx.push_back( x );
    vy.push_back( m.Eval( x ) );
  }

  return new TGraph( vx.size(), vx.data(), vy.data() );
}


usr::Measurement
CalcLumi(
  const double   r,
  const double   offset,
  const double   z,
  const unsigned nrun
  )
{
  std::uniform_real_distribution<double> pdf( 0, 1 );
  std::mt19937_64 gen( usr::CurrentTimeInNanSec() );

  const double tmax
    = atan( 1.5 * ( r + sqrt( 2*xmax*xmax ) + fabs(offset) ) / z );
  const double Sangle
    = 2* M_PI * ( 1-cos( tmax ) );
  unsigned passed = 0;

  for( unsigned i = 0; i < nrun; ++i ){
    if( i % 1000 == 0 ){
      std::cout << boost::format( "\rr=%lf o=%lf z=%lf %d/%d" )
        % r % offset % z % i %nrun << std::flush;
    }

    const double genr = pdf( gen );
    const double gent = 2 * M_PI * pdf( gen );
    const double x0   = offset + r * sqrt( genr ) * sin( gent );
    const double y0   = r * sqrt( genr ) * cos( gent );

    const double genz = cos( tmax ) + ( 1-cos( tmax ) ) * pdf( gen );
    const double genp = 2 * M_PI * pdf( gen );
    const double genq = sqrt( 1-genz*genz );

    const double x1 = x0 + ( z/genz ) *genq * sin( genp );
    const double y1 = y0 + ( z/genz ) *genq * cos( genp );

    if( fabs( x1 ) < xmax/2 && fabs( y1 ) < xmax/2 ){
      passed += 1;
    }
  }

  return M_PI * r * r * Sangle * usr::Efficiency::Bayesian( passed, nrun );

}
