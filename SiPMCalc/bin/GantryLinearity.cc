#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/PlotUtils/interface/PlotCommon.hpp"
#include "UserUtils/PlotUtils/interface/Ratio1DCanvas.hpp"

#include "SiPMCalib/InvSqCalc/interface/LEDFormat.hpp"

#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"

double LinearModel( const double*, const double* );
double LOModel( const double*, const double* );
double NLOModel( const double*, const double* );

double
ZProf( const double* vz, const double* param )
{
  const double z  = vz[0];
  const double z0 = param[0];
  const double o  = param[1];
  const double N  = param[2];
  const double P  = param[3];

  const double D2 = o*o + ( z-z0 )*( z-z0 );
  // return N/D2 + P;
  return ( N*( z-z0 ) )/( D2 * sqrt( D2 ) ) + P;
}

int
main( int argc, char* argv[] )
{
  usr::po::options_description desc( "Options for plot making" );
  desc.add_options()
    ( "data,d", usr::po::value<std::string>(), "SiPM Data file" )
    ( "reference,r", usr::po::value<std::string>(), "Reference LED file" )
    ( "fitz", usr::po::value<double>(), "Low light fit z value" )
    ( "fitn", usr::po::value<double>(), "Low light fit number of photons" )
    ( "npix", usr::po::value<double>(), "Number of pixels on SiPM board" )
    ( "output,o", usr::po::value<std::string>(), "Output file" )
  ;

  usr::ArgumentExtender arg;
  arg.AddOptions( desc );
  arg.ParseOptions( argc, argv );

  const double fitz = arg.Arg<double>( "fitz" );
  const double fitn = arg.Arg<double>( "fitn" );
  const double npix = arg.Arg<double>( "npix" );
  LEDManager reference( arg.Arg( "reference" ) );
  LEDManager data( arg.Arg( "data" ) );
  const std::string output = arg.Arg("output");

  TGraph* referencez = reference.MakeZScanGraph(
    reference._pointlist.front().x,
    reference._pointlist.front().y );
  TGraph* dataz = data.MakeZScanGraph(
    data._pointlist.front().x,
    data._pointlist.front().y );

  // Getting fit for plotting
  TF1 func = TF1( "func", ZProf,
    usr::plt::GetXmin( referencez ), usr::plt::GetXmax( referencez ),
    4 );
  referencez->Fit( &func, "EX0 M E N 0" );
  std::cout << func.GetParameter( 0 ) << " "
           << func.GetParameter( 1 ) << " "
           << func.GetParameter( 2 ) << " "
           << func.GetParameter( 3 ) << std::endl;



  TGraphErrors lumigraph( dataz->GetN() );

  for( int i = 0; i < dataz->GetN(); ++i ){
    const double z           = dataz->GetX()[i];
    const double rawlumi     = func.Eval( z ) - func.GetParameter(3);
    const double reflumi     = func.Eval( fitz ) - func.GetParameter(3);
    const double lumi        = (rawlumi / reflumi) * (fitn / npix) ;
    const double readout     = dataz->GetY()[i] /100;
    const double readout_unc = dataz->GetErrorY( i ) /100 ;


    lumigraph.SetPoint( i, lumi, readout );
    lumigraph.SetPointError( i, 0, readout_unc );
  }

  const double xmin = usr::plt::GetXmin( lumigraph );
  const double xmax = usr::plt::GetXmax( lumigraph );
  const double ymax = usr::plt::GetYmax( lumigraph );

  // Fitting leading order model
  TF1 expf( "ExpF", LOModel,  xmin, 1.3, 4 );
  expf.SetParameter( 0, 1 );
  expf.FixParameter( 1, 1 );
  expf.SetParameter( 2, 0 );
  expf.FixParameter( 3, npix );
  lumigraph.Fit( &expf, "QN0 EX0 R" );
  lumigraph.Fit( &expf, "QN0 EX0 R" );
  lumigraph.Fit( &expf, "QN0 EX0 R" );
  auto fitexp = lumigraph.Fit( &expf, "QN0 EX0 R S" );

  std::cout << "==== LEADING ORDER ====================" << std::endl;
  std::cout << "PDE :"
            << expf.GetParameter( 1 )  << " +- " << expf.GetParError( 1 )
            << std::endl;
  std::cout << "Gain: "
            << expf.GetParameter( 0 ) << " +- " << expf.GetParError( 0 )
            << std::endl;

  TF1 nlof( "NLOF", NLOModel, xmin, xmax, 6 );
  nlof.SetParameter( 0, expf.GetParameter( 0 ) );
  nlof.FixParameter( 1, expf.GetParameter( 1 ) );
  nlof.SetParameter( 2, 0 );
  nlof.SetParameter( 3, 17 );
  nlof.SetParameter( 4, expf.GetParameter( 2 ) );
  nlof.FixParameter( 5, npix );
  lumigraph.Fit( &nlof, "QN0 EX0 R" );
  lumigraph.Fit( &nlof, "QN0 EX0 R" );
  lumigraph.Fit( &nlof, "QN0 EX0 R" );
  auto fitnlo = lumigraph.Fit( &nlof, "QN0 EX0 R S" );

  std::cout << "==== NEXT TO LEADING ORDER ====================" << std::endl;
  std::cout << "PDE :"
            << nlof.GetParameter( 1 )  << " +- " << nlof.GetParError( 1 )
            << std::endl;
  std::cout << "Gain: "
            << nlof.GetParameter( 0 ) << " +- " << nlof.GetParError( 0 )
            << std::endl;
  std::cout << "alpha :"
            << nlof.GetParameter( 2 )  << " +- " << nlof.GetParError( 2 )
            << std::endl;
  std::cout << "beta: "
            << nlof.GetParameter( 3 ) << " +- " << nlof.GetParError( 3 )
            << std::endl;

  usr::plt::Ratio1DCanvas c;

  c.PlotGraph( lumigraph,
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::TrackY( usr::plt::tracky::both )
    );

  auto& fexpg = c.PlotFunc( expf,
    usr::plt::PlotType( usr::plt::fittedfunc ),
    usr::plt::VisualizeError( fitexp ),
    //usr::plt::EntryText(
    //  usr::fstr( " LO_{GOF=%.1lf/%d}", fitexp->Chi2(), fitexp->Ndf() ) ),
    RooFit::Precision( 1e-4 ),
    usr::plt::LineColor( usr::plt::col::red ),
    usr::plt::FillColor( usr::plt::col::pink, 0.3 )
    );

  auto& fnlog = c.PlotFunc( nlof,
    usr::plt::PlotType( usr::plt::fittedfunc ),
    usr::plt::VisualizeError( fitnlo ),
    usr::plt::EntryText(
      usr::fstr( "NLO_{GOF=%.1lf/%d}", fitnlo->Chi2(), fitnlo->Ndf() ) ),
    RooFit::Precision( 1e-4 ),
    usr::plt::LineColor( usr::plt::col::green ),
    usr::plt::FillColor( usr::plt::col::green, 0.3 )
    );


  c.PlotGraph( lumigraph,
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::EntryText( "SiPM Readout" ),
    usr::plt::MarkerColor( usr::plt::col::black ),
    usr::plt::MarkerStyle( usr::plt::sty::mkrcircle  ),
    usr::plt::MarkerSize( 0.1 ),
    usr::plt::LineColor( usr::plt::col::gray, 0.5 ),
    usr::plt::TrackY( usr::plt::tracky::both )
    );
  c.TopPad().SetDataMax( ymax * 31 );// Expanding y range


  c.PlotScale( fnlog,     fnlog, usr::plt::PlotType( usr::plt::fittedfunc ) );
  c.PlotScale( fexpg,     fnlog, usr::plt::PlotType( usr::plt::simplefunc ) );
  c.PlotScale( lumigraph, fnlog, usr::plt::PlotType( usr::plt::scatter ) );

  c.DrawLuminosity( "Pulsed LED" );
  c.DrawCMSLabel( "", "Linearity Test" );
  c.TopPad()
  .WriteLine( usr::fstr("N_{pixel} = %d", (int)npix) )
  .WriteLine( "V_{bias} = 54V" )
  .WriteLine( usr::fstr( "P_{sec} = %.1lf#pm%.1lf%%",
    nlof.GetParameter( 2 )*100,
    nlof.GetParError( 2 )*100 ) )
  .WriteLine( usr::fstr( "#tau_{SiPM}/#tau_{pulse} = %.1lf#pm%.1lf",
    nlof.GetParameter( 3 ),
    nlof.GetParError( 3 ) ) );

  c.BottomPad().Xaxis().SetTitle( "#epsilon N_{#gamma} / N_{pix}" );
  c.TopPad().Yaxis().SetTitle( "Measured Intensity (A.U.)" );
  c.BottomPad().Yaxis().SetTitle( "Data/Fit" );
  c.BottomPad().SetLogx( 1 );
  c.TopPad().SetLogx( 1 );
  c.TopPad().SetLogy( 1 );

  c.SaveAsPDF( output );
  std::cout << "Saved to " << output << std::endl;

  return 0;
}
