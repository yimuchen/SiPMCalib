#include "SiPMCalib/SiPMCalc/interface/SiPMDarkPdf.hpp"
#include "UserUtils/PlotUtils/interface/Ratio1DCanvas.hpp"

#include "RooChi2Var.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooRealVar.h"

#include <algorithm>
#include <boost/format.hpp>
#include <fstream>
#include <iostream>

int
main( int argc, char const* argv[] )
{
  std::vector<double> vecdata;
  std::ifstream input( argv[1], std::ios::in );

  char buffer[65536];
  double t, peak, tarea, xin;
  input.getline( buffer, 65536 );
  input.getline( buffer, 65536 );// Getting dummy lines

  while( input >> t >> peak >> tarea >> xin ){
    std::cout << "\r"
              << boost::format( "%lf\t%lf\t%lf\t%lf" ) % t % peak % tarea % xin
              << std::flush;
    if( peak == 0.200000000 || peak == 0.100000000 ){
      continue;
    }
    vecdata.push_back( xin * -1000 );
  }

  std::sort( vecdata.begin(), vecdata.end() );

  RooRealVar x( "x", "Readout (mV#times ns)",
                vecdata.front() - 50,
                vecdata.back() + 50 );
  x.setBins( 250 );

  RooDataHist data( "data", "data", RooArgList( x ) );
  RooDataSet udata( "udata", "udata", RooArgList( x ) );

  for( const auto p : vecdata ){
    x = p;
    data.add( RooArgSet( x ) );
    udata.add( RooArgSet( x ) );
  }

  RooRealVar ped( "ped", "ped", 0, -100, 100 );
  RooRealVar gain( "gain", "gain", 150, 1, 500 );
  RooRealVar s0( "s0", "s0", 20, 10, 100 );
  RooRealVar s1( "s1", "s1", 10, 0.001, 50 );
  RooRealVar dcfrac( "dcfrac1", "dcfrac1", 0.005, 0, 0.2 );

  SiPMDarkPdf pdf( "dark", "dark", x, ped, gain, s0, s1, dcfrac );

  auto fit = pdf.fitTo( data, RooFit::Range( -100, 400 ), RooFit::Save() );

  unsigned it = 0 ;
  while( fit->status() && it < 10 ){
    delete fit;
    fit = pdf.fitTo( data, RooFit::Save(),
      RooFit::Range( -100, 400 ) );
    it++;
    std::cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << it << std::endl;
  }

  usr::plt::Ratio1DCanvas c( x );
  auto & fitgraph = c.PlotPdf( pdf,
    RooFit::Normalization( data.sumEntries() ),
    usr::plt::EntryText( "Model Fit" ),
    RooFit::Range(-100,600) );
  auto & datgraph = c.PlotData( udata,
    usr::plt::EntryText( "SiPM readout" ) );

  fitgraph.SetLineColor( kBlue );
  fitgraph.SetFillColor( kCyan );
  datgraph.SetMarkerSize( 0.2 );

  auto & r1 = c.PlotScale( datgraph, fitgraph,
    usr::plt::PlotType( usr::plt::scatter ) );

  c.TopPad().SetHistAxisTitles( "Readout", "mV #times ns" );
  c.BottomPad().Yaxis().SetTitle( "Data/fit" );
  c.SetLogy( true );
  c.SaveAsPNG( "darktest.png" );
  c.SetLogy( false );
  c.SaveAsPNG( "darktest_nlog.png" );

  std::cout
    << "ped " << ped.getVal() << "  " << ped.getError() << std::endl
    << "gain " << gain.getVal() << "  " << gain.getError() << std::endl
    << "s0 " << s0.getVal() << "  " << s0.getError() << std::endl
    << "s1 " << s1.getVal() << "  " << s1.getError() << std::endl
    << "dcfrac " << dcfrac.getVal() << "  " << dcfrac.getError() << std::endl
    << std::endl;

  return 0;
}
