#include "SiPMCalib/SiPMCalc/interface/SiPMPdf.hpp"
#include "SiPMCalib/SiPMCalc/interface/SiPMSumGaus.hpp"

#include "UserUtils/PlotUtils/interface/Ratio1DCanvas.hpp"

#include "RooDataSet.h"
#include "RooRealVar.h"



int
main()
{
  RooRealVar x( "x", "Readout (mV#times ns)", 300, 1500 );

  RooRealVar ped( "ped", "ped", 350, 100, 500 );
  RooRealVar gain( "gain", "gain", 200, 100, 500 );
  RooRealVar s0( "s0", "s0", 50, 15, 100 );
  RooRealVar s1( "s1", "s1", 10, 0.05, 100 );
  RooRealVar mean( "mean", "mean", 1.5, 0.01, 50 );
  RooRealVar lambda( "lambda", "lambda", 0.25, 0, 0.2 );
  RooRealVar alpha( "alpha", "alpha", 0.1, 0, 1 );
  RooRealVar beta( "beta", "beta", 50, 20, 1000 );
  SiPMPdf p0( "p0", "p0", x, ped, gain, s0, s1, mean, lambda, alpha, beta );

  RooRealVar ped1( "ped", "ped", 350, 100, 500 );
  RooRealVar gain1( "gain", "gain", 200, 100, 500 );
  RooRealVar s01( "s0", "s0", 50, 15, 100 );
  RooRealVar s11( "s1", "s1", 10, 0.05, 100 );
  RooRealVar mean1( "mean", "mean", 1.5, 0.01, 50 );
  RooRealVar lambda1( "lambda", "lambda", 0.25, 0, 0.2 );
  SiPMPdf p1( "p1", "p1", x, ped1, gain1, s01, s11, mean1, lambda1 );

  RooRealVar ped2( "ped", "ped", 350, 100, 500 );
  RooRealVar gain2( "gain", "gain", 200, 100, 500 );
  RooRealVar s02( "s0", "s0", 50, 15, 100 );
  RooRealVar s12( "s1", "s1", 10, 0.05, 100 );
  RooRealVar mean2( "mean", "mean", 1.5, 0.01, 50 );
  SiPMPdf p2( "p2", "p2", x, ped2, gain2, s02, s12, mean2 );

  RooDataSet* mock = p0.generate( RooArgSet( x ), 3e4 );

  p1.fitTo( *mock );
  p2.fitTo( *mock );

  usr::plt::Ratio1DCanvas c( x );
  auto& g2 = c.PlotPdf(
    p2, RooFit::Normalization( 3e4 ),
    usr::plt::EntryText( "Poisson" ) );
  auto& g1 = c.PlotPdf(
    p1, RooFit::Normalization( 3e4 ),
    usr::plt::EntryText( "Generalized Poisson" ) );
  auto& g0 = c.PlotPdf(
    p0, RooFit::Normalization( 3e4 ),
    usr::plt::EntryText( "Gen. Poisson w/ AP" ) );
  auto& gdat = c.PlotData( mock, usr::plt::EntryText("Fake data") );

  gdat.SetMarkerSize(0.08);
  g1.SetLineColor(kRed);
  g2.SetLineColor(kGreen);

  c.PlotScale( g0, g0, usr::plt::TrackY(usr::plt::TrackY::none),
             usr::plt::PlotType(usr::plt::fittedfunc) );
  c.PlotScale( g2, g0, usr::plt::TrackY(usr::plt::TrackY::none),
             usr::plt::PlotType(usr::plt::fittedfunc) );
  c.PlotScale( g1, g0, usr::plt::TrackY(usr::plt::TrackY::none),
             usr::plt::PlotType(usr::plt::fittedfunc) );
  c.PlotScale( gdat, g0, usr::plt::TrackY(usr::plt::TrackY::none),
             usr::plt::PlotType(usr::plt::scatter) );

  c.TopPad().SetHistAxisTitles("Pulse Area", "mV #times ns");
  c.BottomPad().Yaxis().SetTitle("/Gen.with.AP");
  c.SetLogy(true);

  c.SaveAsPNG( "ModelCompare.png" );

  delete mock;
}
