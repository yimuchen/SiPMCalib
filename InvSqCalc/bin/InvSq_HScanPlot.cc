#include "SiPMCalib/InvSqCalc/interface/LEDFormat.hpp"
#include "SiPMCalib/InvSqCalc/interface/MCFormat.hpp"

#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/PlotUtils/interface/Flat2DCanvas.hpp"


#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF2.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TFitResult.h"

#include "boost/format.hpp"

double ExpFunc( const double* xy, const double * param )
{
  const double x = xy[0];
  const double y = xy[1];
  const double x0 = param[0];
  const double y0 = param[1];
  const double z  = param[2];
  const double N  = param[3];
  const double P  = param[4];

  const double D2 = (x-x0)*(x-x0) + (y-y0)*(y-y0) + z*z;
  // return N/D2 + P;
  return (N*z)/(D2 * sqrt(D2)) + P;
}

int main( int argc, char** argv )
{
  usr::po::options_description desc("Options for plot making");
  desc.add_options()
    ("data,d", usr::po::value<std::string>(), "Data file")
    ("output,o", usr::po::value<std::string>(), "Output file")
  ;

  usr::ArgumentExtender arg;
  arg.AddOptions( desc );
  arg.ParseOptions( argc, argv );

  LEDManager data( arg.Arg("data") );
  const double z = data._pointlist.front().z ;

  TH2D* hist = data.MakeHScanGraph( z );
  TF2* func1 = new TF2("func1", ExpFunc, data.Xmin(), data.Xmax(), data.Ymin(), data.Ymax(), 5 );

  func1->SetParameters(
    (data.Xmin() + data.Xmax())/2,
    (data.Ymin() + data.Ymax())/2,
    z ,
    ( data.LumiMax() - data.LumiMin() )* (z*z) ,
    data.LumiMin()
  );

  auto fit = hist->Fit( func1,"EX0 N 0 S");

  TGraphErrors cen;
  cen.SetPoint( 0, func1->GetParameter(0), func1->GetParameter(1) );
  cen.SetPointError( 0, func1->GetParError(0), func1->GetParError(1) );

  usr::plt::Flat2DCanvas c;


  c.PlotHist( hist,
    usr::plt::Plot2DF(usr::plt::heatcont),
    usr::plt::EntryText( "Data (est. contour)" ) );

  auto& funcg = c.PlotFunc( func1,
    usr::plt::Plot2DF(usr::plt::cont),
    usr::plt::EntryText( "Fitted Contour") );

  c.Plot1DGraph( cen,
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::EntryText( "Fitted (x_{0},y_{0})") );

  c.DrawCMSLabel("Preliminary","HGCal");
  c.DrawLuminosity("LED setup");
  c.Pad().SetTextCursor( 0.015, 0.5, usr::plt::font::top_left )
    .WriteLine( "#frac{N z_{0}}{ #sqrt{(x-x_{0})^{2} + (y-y_{0})^{2} + z_{0}^{2} }^{3}} + P" )
    .WriteLine( "" )
    .WriteLine( "" )
    .WriteLine(
      (boost::format( "x_{0} = %.1lf_{#pm%.2lf} [mm]")
      % func1->GetParameter(0)
      % func1->GetParError(0)).str() )
    .WriteLine(
      (boost::format( "y_{0} = %.1lf_{#pm%.2lf} [mm]")
      % func1->GetParameter(1)
      % func1->GetParError(1)).str() )
    .WriteLine(
      (boost::format( "z_{0} = %.1lf_{%.2lf} [mm]")
      % func1->GetParameter(2)
      % func1->GetParError(2)
      ).str()  )
    .WriteLine( (boost::format("Gantry Z = %.1lf [mm]") % z).str() )
    .WriteLine( (boost::format("#chi^{2}/DoF = %.3lf")% (fit->Chi2()/fit->Ndf())).str())
    ;

  hist->SetLineColorAlpha( kBlue, 1 );
  funcg.SetLineColorAlpha( kGreen, 1 );
  funcg.SetLineWidth(1);

  cen.SetLineColor(kRed);
  cen.SetMarkerColor(kRed);
  cen.SetMarkerStyle(20);
  cen.SetMarkerSize(0.5);

  c.Xaxis().SetTitle( "Gantry x [mm]");
  c.Yaxis().SetTitle( "Gantry y [mm]");
  c.Zaxis().SetTitle( "Luminosity [A.U.]");

  c.SaveAsPDF( arg.Arg("output") );

  return 0;
}