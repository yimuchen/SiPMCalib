#include "SiPMCalib/SiPMCalc/interface/SiPMPdf.hpp"
#include "UserUtils/PlotUtils/interface/Ratio1DCanvas.hpp"

#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooChi2Var.h"

#include <boost/format.hpp>
#include <fstream>
#include <iostream>
#include <algorithm>


/// Getting better results:
// * Pedestal subtraction.
// * Hardware issues (double check readout requirements)

int
main( int argc, char const* argv[] )
{
  std::vector<double> vecdata;
  std::ifstream input( argv[1], std::ios::in );

  char buffer[65536];
  double t, peak, tarea, xin ;
  input.getline( buffer, 65536 );
  input.getline( buffer, 65536 );// Getting dummy lines

  while( input >> t >> peak >> tarea >> xin ){
    std::cout << "\r"
              << boost::format("%lf\t%lf\t%lf\t%lf") %t%peak%tarea%xin
              << std::flush;
    if( peak == 0.200000000 || peak == 0.100000000 ) { continue; }
    vecdata.push_back( xin*-1000 );
  }
  std::sort( vecdata.begin(), vecdata.end() );



  RooRealVar x( "x", "Readout (mV#times ns)",
      vecdata.front()-50,
      vecdata.back()+50 );
  x.setBins( 500 );
  x.setRange( "full", x.getMin(), x.getMax() );
  x.setRange( "Wide", -50, x.getMax() * 3 / 4 );

  RooDataHist data( "data", "data", RooArgList( x ) );
  RooDataSet  udata("udata","udata", RooArgList( x ) );

  for( const auto p : vecdata ){
    x = p;
    data.add( RooArgSet(x) );
    udata.add( RooArgSet(x) );
  }



  RooRealVar ped( "ped", "ped", 10, -100, 1500 );
  RooRealVar gain( "gain", "gain", 320, 0, 10000 );
  RooRealVar s0( "s0", "s0", 60, 10, 100000 );
  RooRealVar s1( "s1", "s1", 20, 0.001, 100 );
  RooRealVar mean( "mean", "mean", 7.7, 0.0001, 50 );
  RooRealVar lambda( "lambda", "lambda", 0.01, 0, 0.50 );
  RooRealVar acfrac("acfrac", "acfrac", 0.01,0,0.4);
  RooRealVar acshift("acshift","acshift",10,5,50);
  RooRealVar acwidth("acwidth","acwidth",20,0,100);
  RooRealVar alpha("alpha","alpha",0.01,0,0.5);
  RooRealVar beta("beta","beta",100,1,2000);
  SiPMPdf p0( "p0", "p", x, ped, gain, s0, s1, mean, lambda
    , acfrac, acshift, acwidth
    , alpha,beta );
  //  );

  RooRealVar ped1( "ped1", "ped", 10, -50, 1500 );
  RooRealVar gain1( "gain1", "gain", 500, 0, 10000 );
  RooRealVar s01( "s01", "s0", 50, 10, 1000 );
  RooRealVar s11( "s11", "s1", 0, 0, 100 );
  RooRealVar mean1( "mean1", "mean", 0.6, 0.0001, 50 );
  RooRealVar lambda1( "lambda1", "lambda", 0.0, 0, 0.5 );
  SiPMPdf p1( "p1", "p1", x, ped1, gain1, s01, s11, mean1, lambda1 );



  RooFitResult* fit0 = p0.fitTo( data, RooFit::Save(),
    RooFit::Range( "Wide" ) );

  unsigned it = 0 ;
  while( fit0->status() && it < 10 ){
    delete fit0;
    fit0 = p0.fitTo( data, RooFit::Save(),
      RooFit::Range( "Wide" ) );
    it++;
    std::cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << it << std::endl;
  }

  //delete fit0 ;
  //fit0 = p0.fitTo( udata, RooFit::Save(),
  //  RooFit::Range( "Wide" ) );
//
//  while( fit0->status() ){
//    delete fit0;
//    fit0 = p0.fitTo( udata, RooFit::Save(),
//      RooFit::Range( "Wide" ) );
//  }

  // Evaluating second fit range from first fit
  x.setRange("Narrow",
    x.getMin("Wide"),
    x.getMax("Wide")
     );
  ped1 = ped.getVal();
  gain1 = gain.getVal();
  s01 = s0.getVal();
  s11 = s1.getVal();
  mean1 = mean.getVal();
  lambda1 = lambda.getVal();

  RooFitResult* fit1 = p1.fitTo(
    data, RooFit::Save(),
    RooFit::Range( "Narrow" ) );

  it = 0 ;
  while( fit1->status() && it < 10 ){
    delete fit1;
    fit1 = p1.fitTo(
      data, RooFit::Save(),
      RooFit::Range( "Narrow" ) );
    ++it;
  }

  RooChi2Var gof("gof","gof",p0,data,RooFit::Range("Wide"));
  RooChi2Var gof1("gof1","gof1",p1,data,RooFit::Range("Narrow"));
  const int dof = (x.getMax("Wide") - x.getMin("Wide"))
                   / x.getBinning().averageBinWidth()
                   - 6;
  const int dof1 = (x.getMax("Narrow") - x.getMin("Narrow"))
                   / x.getBinning().averageBinWidth()
                   - 6;
  const double gofv = gof.getVal() / dof;
  const double gofv1 = gof1.getVal() / dof1 ;

  usr::plt::Ratio1DCanvas c( x );
  auto& fit1graph = c.PlotPdf( p1,
    RooFit::Normalization( data.sumEntries() ),
    RooFit::NormRange( "full" ),
    RooFit::Range("Narrow"),
    RooFit::VisualizeError( *fit1, 1 ),
    RooFit::Precision(-1),
    usr::plt::EntryText(
      (boost::format("Zero corr.(#chi^{2}_{#nu}=%.2d) ")
        %(gofv1)).str() ) );
  auto& fit0graph = c.PlotPdf( p0,
    RooFit::Normalization( data.sumEntries() ),
    RooFit::NormRange( "full" ),
    RooFit::Range("Wide"),
    RooFit::VisualizeError( *fit0, 1 ),
    RooFit::Precision(-1),
    usr::plt::EntryText(
      (boost::format("Poi x Gaus (#chi^{2}_{#nu}=%.2d) ")
        %(gofv)).str() ) );
  auto& datgraph = c.PlotData( data,
    usr::plt::EntryText(         "SiPM readout" ) );

  fit1graph.SetLineColor( kRed );
  fit1graph.SetFillColor( kPink -9 );
  fit0graph.SetLineColor( kBlue );
  fit0graph.SetFillColor( kCyan );
  datgraph.SetMarkerSize( 0.2 );

  auto& r1 = c.PlotScale( datgraph, fit1graph,
    usr::plt::PlotType( usr::plt::scatter ) );
  auto& r0 = c.PlotScale( datgraph, fit0graph,
    usr::plt::PlotType( usr::plt::scatter ) );

  r1.SetLineColor( kRed );
  r1.SetMarkerColor( kRed );
  r0.SetLineColor( kBlue );
  r0.SetMarkerColor( kBlue );


  c.TopPad().SetHistAxisTitles( "Readout", "mV #times ns" );
  c.BottomPad().Yaxis().SetTitle( "Data/fit" );
  c.SetLogy( true );
  c.SaveAsPNG( "test.png" );
  c.SetLogy( false );
  c.SaveAsPNG( "test_nlog.png" );


  std::cout
    << "ped " << ped.getVal() <<  "  " << ped.getError() << std::endl
    << "gain " << gain.getVal() <<  "  " << gain.getError() << std::endl
    << "s0 " << s0.getVal() <<  "  " << s0.getError() << std::endl
    << "s1 " << s1.getVal() <<  "  " << s1.getError() << std::endl
    << "mean " << mean.getVal() <<  "  " << mean.getError() << std::endl
    << "lambda " << lambda.getVal() <<  "  " << lambda.getError() << std::endl
    << "acfrac " << acfrac.getVal() <<  "  " << acfrac.getError() << std::endl
    << "acshift " << acshift.getVal() <<  "  " << acshift.getError() << std::endl
    << "acwidth " << acwidth.getVal() <<  "  " << acwidth.getError() << std::endl
    << "alpha " << alpha.getVal() <<  "  " << alpha.getError() << std::endl
    << "beta " << beta.getVal() <<  "  " << beta.getError() << std::endl
    << gof.getVal() << "/" << dof << " = " << gof.getVal() / dof << std::endl
  ;
  std::cout << std::endl;
  std::cout
    << "ped1 " << ped1.getVal() <<  "  " << ped1.getError() << std::endl
    << "gain1 " << gain1.getVal() <<  "  " << gain1.getError() << std::endl
    << "s01 " << s01.getVal() <<  "  " << s01.getError() << std::endl
    << "s11 " << s11.getVal() <<  "  " << s11.getError() << std::endl
    << "mean1 " << mean1.getVal() <<  "  " << mean1.getError() << std::endl
    << "lambda1 " << lambda1.getVal() <<  "  " << lambda1.getError() << std::endl
    << gof1.getVal() << "/" << dof1 << " = " << gof1.getVal() / dof1 << std::endl
  ;
  std::cout << std::endl;

  return 0;
}
