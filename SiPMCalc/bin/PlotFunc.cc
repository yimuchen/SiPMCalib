#include "SiPMCalib/SiPMCalc/interface/SiPMPdf.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

int
main( int argc, char const* argv[] )
{

  RooRealVar x( "x", "Readout (mV#times ns)", -300, 1700 );
  RooRealVar ped( "ped", "ped", 0, -50, 1500 );
  RooRealVar gain( "gain", "gain", 300, 0, 10000 );
  RooRealVar s0( "s0", "s0", 20, 10, 100000 );
  RooRealVar s1( "s1", "s1",  5, 0.001, 100 );
  RooRealVar mean( "mean", "mean", 1.2, 0.0001, 50 );
  RooRealVar lambda( "lambda", "lambda", 0.10, 0, 0.50 );
  RooRealVar acfrac("acfrac","acfrac",0,0,0.1);
  RooRealVar acshift("acshift","acshift",0,0,0.1);
  RooRealVar acwidth("acwidth","acwidth",1,0.5,2);
  RooRealVar alpha("alpha","alpha",0.1,0,1);
  RooRealVar beta("beta","beta",30,5,1000);
  SiPMPdf p0( "p0", "p0",
    x, ped, gain, s0, s1, mean, lambda,
    acfrac, acshift, acwidth,
    alpha,beta );


  usr::plt::Simple1DCanvas c(x);

  auto& g0 = c.PlotPdf( p0, usr::plt::TrackY(usr::plt::TrackY::both) );

  beta = 50;

  auto& g1 = c.PlotPdf( p0, usr::plt::TrackY(usr::plt::TrackY::both) );

  g0.SetLineColor(kBlue);
  g1.SetLineColor(kGreen);

  c.SetLogy( true );
  c.Pad().SetYaxisMin(1e-12);

  c.SaveAsPNG("funcshape.png");

  return 0;
}
