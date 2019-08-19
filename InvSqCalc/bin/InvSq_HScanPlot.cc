#include "SiPMCalib/InvSqCalc/interface/LEDFormat.hpp"
#include "SiPMCalib/InvSqCalc/interface/MCFormat.hpp"

#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"
#include "UserUtils/PlotUtils/interface/Flat2DCanvas.hpp"


#include "TCanvas.h"
#include "TColor.h"
#include "TF2.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"

#include "boost/format.hpp"

double
ExpFunc( const double* xy, const double* param )
{
  const double x  = xy[0];
  const double y  = xy[1];
  const double x0 = param[0];
  const double y0 = param[1];
  const double z  = param[2];
  const double N  = param[3];
  const double P  = param[4];

  const double D2 = ( x-x0 )*( x-x0 ) + ( y-y0 )*( y-y0 ) + z*z;
  // return N/D2 + P;
  return ( N*z )/( D2 * sqrt( D2 ) ) + P;
}

int
main( int argc, char** argv )
{
  usr::po::options_description desc( "Options for plot making" );
  desc.add_options()
    ( "data,d", usr::po::value<std::string>(), "Data file" )
    ( "output,o", usr::po::value<std::string>(), "Output file" )
    ( "type,t", usr::po::value<std::string>()->default_value( "static" ), "Type of LED configuration" )
  ;

  usr::ArgumentExtender arg;
  arg.AddOptions( desc );
  arg.ParseOptions( argc, argv );

  LEDManager data( arg.Arg( "data" ) );
  const double z = data._pointlist.front().z;

  TH2D* hist = data.MakeHScanGraph( z );
  TF2* func1 = new TF2( "func1", ExpFunc,
    data.Xmin(), data.Xmax(),
    data.Ymin(), data.Ymax(), 5 );

  func1->SetParameters(
    hist->GetMean( 1 ),
    hist->GetMean( 2 ),
    z,
    ( data.LumiMax() - data.LumiMin() )* ( z*z ),
    data.LumiMin()
    );

  std::cout << hist->GetMean( 1 ) << std::endl;
  std::cout << hist->GetMean( 2 ) << std::endl;


  auto fit = hist->Fit( func1, "EX0 N 0 S" );

  const double xfit    = func1->GetParameter( 0 );
  const double yfit    = func1->GetParameter( 1 );
  const double zfit    = func1->GetParameter( 2 );
  const double Nfit    = func1->GetParameter( 3 );
  const double Pfit    = func1->GetParameter( 4 );
  const double sipmlen = 1.3;

  TGraphErrors cen;
  cen.SetPoint( 0, xfit, yfit );
  cen.SetPointError( 0, func1->GetParError( 0 ), func1->GetParError( 1 ) );

  TGraph sipm( 5 );
  sipm.SetPoint( 0, xfit + sipmlen/2, yfit+sipmlen/2 );
  sipm.SetPoint( 1, xfit - sipmlen/2, yfit+sipmlen/2 );
  sipm.SetPoint( 2, xfit - sipmlen/2, yfit-sipmlen/2 );
  sipm.SetPoint( 3, xfit + sipmlen/2, yfit-sipmlen/2 );
  sipm.SetPoint( 4, xfit + sipmlen/2, yfit+sipmlen/2 );

  // gStyle->SetPalette(kBird);

  usr::plt::Flat2DCanvas c;

  c.PlotHist( hist,
    usr::plt::Plot2DF( usr::plt::heat ) );

  auto& histc = c.PlotHist( (TH2D*)hist->Clone(),
    usr::plt::Plot2DF( usr::plt::cont ),
    usr::plt::EntryText( "Data (est. contour)" ),
    usr::plt::LineColor( usr::plt::col::blue ) );

  auto& funcg = c.PlotFunc( func1,
    usr::plt::Plot2DF( usr::plt::cont ),
    usr::plt::EntryText( "Fitted Contour" ),
    usr::plt::LineColor( usr::plt::col::green ) );

  c.Plot1DGraph( sipm,
    usr::plt::PlotType( usr::plt::simplefunc ),
    usr::plt::EntryText( "SiPM (Expected)" ),
    usr::plt::LineColor( usr::plt::col::red ),
    usr::plt::LineStyle( usr::plt::sty::lindashed ) );

  c.Plot1DGraph( cen,
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::EntryText( "Fitted (x_{0},y_{0})" ),
    usr::plt::MarkerColor( usr::plt::col::red ),
    usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
    usr::plt::MarkerSize( 0.5 ),
    usr::plt::LineColor( usr::plt::col::red ) );

  // Styling of contour levels
  std::vector<double> contlevel = {
    Nfit / ( zfit*zfit ) + Pfit,
    Nfit / ( zfit*zfit ) * 99./100. + Pfit,
    Nfit / ( zfit*zfit ) * 90./100. + Pfit,
    Nfit / ( zfit*zfit ) * 75./100. + Pfit,
    Nfit / ( zfit*zfit ) * 50./100. + Pfit
  };
  std::sort( contlevel.begin(), contlevel.end() );

  funcg.GetHistogram()->SetContour( contlevel.size(), contlevel.data() );
  histc.SetContour( contlevel.size(), contlevel.data() );



  c.DrawCMSLabel( "Preliminary", "HGCal" );
  c.DrawLuminosity( usr::fstr("LED (%s)", arg.Arg("type") ) );
  c.Pad().SetTextCursor( 0.015, 0.6, usr::plt::font::top_left )
  .WriteLine( "Contours: #frac{99}{100}, #frac{9}{10}, #frac{3}{4}, #frac{1}{2}lumi" )
  .WriteLine( "#frac{N z_{0}}{ #sqrt{(x-x_{0})^{2} + (y-y_{0})^{2} + z_{0}^{2} }^{3}} + P" )
  .WriteLine( usr::fstr( "x_{0} = %.1lf_{#pm%.2lf} [mm]",
    func1->GetParameter( 0 ),
    func1->GetParError( 0 ) ) )
  .WriteLine( usr::fstr( "y_{0} = %.1lf_{#pm%.2lf} [mm]",
    func1->GetParameter( 1 ),
    func1->GetParError( 1 ) ) )
  .WriteLine( usr::fstr( "z_{0} = %.1lf_{#pm%.2lf} [mm]",
    func1->GetParameter( 2 ),
    func1->GetParError( 2 ) ) )
  .WriteLine( usr::fstr( "Gantry Z = %.1lf [mm]", z ) )
  .WriteLine( usr::fstr( "#chi^{2}/DoF = %.3lf", fit->Chi2()/fit->Ndf() ) );

  c.Xaxis().SetTitle( "Gantry x [mm]" );
  c.Yaxis().SetTitle( "Gantry y [mm]" );
  c.Zaxis().SetTitle( "Luminosity [pA]" );

  cen.GetHistogram()->SetMaximum( 1000 );
  cen.GetHistogram()->SetMinimum( 1000 );

  c.SaveAsPDF( arg.Arg( "output" ) );
  c.SaveToROOT( "test.root", "c" );

  return 0;
}
