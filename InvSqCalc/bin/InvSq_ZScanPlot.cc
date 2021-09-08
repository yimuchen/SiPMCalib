#include "SiPMCalib/Common/interface/StdFormat.hpp"

#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/Common/interface/STLUtils/OStreamUtils.hpp"
#include "UserUtils/PlotUtils/interface/Ratio1DCanvas.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

#include "TFitResult.h"
#include "TGraphErrors.h"

double
ZProf( const double* vz, const double* param )
{
  const double z  = vz[0];
  const double z0 = param[0];
  const double o  = param[1];
  const double N  = param[2];
  const double P  = param[3];
  const double N0 = param[4];

  const double D2 = o*o + ( z-z0 )*( z-z0 );
  // return N/D2 + P;
  return ( N*N0*( z-z0 ) )/( D2 * sqrt( D2 ) ) + P;
}

TGraph*     MakeZScanGraph( const StdFormat&, double uncscale );
TFitResult* FitAndShiftData( TGraph* data, TF1* func,
                             double& ped, double& zoffset );

int
main( int argc, char** argv )
{
  usr::po::options_description desc( "Options for plot making" );
  desc.add_options()
    ( "data,d", usr::po::value<std::string>(), "Data file" )
    ( "output,o", usr::po::value<std::string>(), "Output plot file" )
    ( "norm,n", usr::po::defvalue<double>( 1.0 ),
    "Estimated norm in the inverse square formula, this helps with the "
    "convergence of fitter" )
    ( "uncscale,u", usr::po::defvalue<double>( 1 ),
    "Additional uncertainty scaling factor" )
  ;

  usr::ArgumentExtender arg;
  arg.AddOptions( desc );
  arg.ParseOptions( argc, argv );

  double pedestal = 0;// For storing the original fit results
  double zoffset  = 0;  // For storing the original fit results

  StdFormat data( arg.Arg( "data" ) );
  TGraph* dataz     = MakeZScanGraph( data, arg.Arg<double>( "uncscale" ) );
  const double xmin = usr::plt::GetXmin( dataz );
  const double xmax = usr::plt::GetXmax( dataz );
  TF1* func         = new TF1( "func", ZProf, xmin, xmax, 5 );
  func->FixParameter( 4, arg.Arg<double>( "norm" ) );
  TFitResult* fit = FitAndShiftData( dataz, func, pedestal, zoffset );

  // Begin plotting
  usr::plt::Ratio1DCanvas c;
  auto& fitg = c.PlotFunc( func,
    usr::plt::PlotType( usr::plt::fittedfunc ),
    usr::plt::VisualizeError( fit ),
    usr::plt::EntryText( "Fitted Data" ),
    usr::plt::TrackY( usr::plt::tracky::max ),
    usr::plt::LineColor( usr::plt::col::blue ),
    usr::plt::FillColor( usr::plt::col::cyan ) );
  auto& datag = c.PlotGraph( dataz,
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::EntryText( "Readout" ),
    usr::plt::MarkerColor( usr::plt::col::black ),
    usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
    usr::plt::MarkerSize( 0.3 ),
    usr::plt::LineColor( usr::plt::col::gray, 0.5 ) );
  c.PlotScale( fitg, fitg,
    usr::plt::PlotType( usr::plt::fittedfunc ),
    usr::plt::FillColor( usr::plt::col::cyan ) );
  c.PlotScale( datag, fitg,
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::MarkerSize( 0.3 ) );

  c.TopPad().Yaxis().SetTitle( "Luminosity - Ped. [mV-ns]" );
  c.BottomPad().Xaxis().SetTitle( "Gantry z + z_{0} [mm]" );
  c.BottomPad().Yaxis().SetTitle( "Data/Fit" );

  c.TopPad().SetLogy( kTRUE );
  c.TopPad().SetLogx( kTRUE );
  c.TopPad().SetYaxisMin( usr::plt::GetYmin( dataz ) / 30  );
  c.BottomPad().SetLogx( kTRUE );

  c.DrawLuminosity( "LED setup" );
  c.DrawCMSLabel( "Preliminary", "HGCal" );
  c.TopPad()
  .SetTextCursor( c.TopPad().InnerTextLeft(), 0.45 )
  .WriteLine( "#frac{L_{0}(z+z_{0})}{(x_{0}^{2} + (z+z_{0})^{2})^{3/2}} + Ped" )
  // .WriteLine( "" )
  .WriteLine( ( boost::format( "Fit Ped = %.2lf_{#pm%.3lf}" )
                % pedestal % func->GetParError( 3 ) ).str() )
  .WriteLine( ( boost::format( "Fit z_{0} = %.2lf_{#pm%.3lf}" )
                % -zoffset % func->GetParError( 0 ) ).str() )
  .WriteLine( ( boost::format( "Fit x_{0} = %.2lf_{#pm%.3lf}" )
                % func->GetParameter( 1 ) % func->GetParError( 1 ) ).str() )
  .WriteLine( ( boost::format( "#chi^{2}/D.o.F = %.2lf" )
                % ( fit->Chi2()/fit->Ndf() ) ).str() );

  c.SaveAsPDF( arg.Arg( "output" ) );

  return 0;
}


TGraph*
MakeZScanGraph( const StdFormat& sformat, const double uncscale )
{
  std::vector<double> z       = sformat.Z();
  std::vector<double> lumi    = sformat.DataCol( 0 );
  std::vector<double> lumierr = sformat.DataCol( 1 );
  std::vector<double> zero( z.size(), 0.0 );
  std::vector<double> lerr_scaled( z.size(), 0.0 );
  lerr_scaled.resize( z.size() );

  std::transform( lumierr.begin(), lumierr.end(), lerr_scaled.begin(),
    [uncscale](double x ){
    return x * uncscale;
  } );

  return new TGraphErrors( z.size(),
    z.data(), lumi.data(),
    zero.data(), lerr_scaled.data() );
}


TFitResult*
FitAndShiftData( TGraph* data, TF1* func, double& ped, double& zoffset )
{
  usr::log::PrintLog( usr::log::INFO,
    "Getting the estimate values for fitting" );
  const double lumimax = usr::plt::GetYmax( data );
  const double lumimin = usr::plt::GetYmin( data );
  const double zmin    = usr::plt::GetXmin( data );

  func->SetParameter( 0, 0 );
  func->SetParameter( 1, 0 );
  func->SetParameter( 2, ( lumimax*zmin*zmin )/func->GetParameter( 4 ) );
  func->SetParameter( 3, lumimin );

  usr::log::PrintLog( usr::log::INFO, "Running the fit" );
  data->Fit( func, "Q EX0 M E N 0" ).Get();

  // Saving the original fit results for the pedestal and zoffset
  ped     = func->GetParameter( 3 );
  zoffset = func->GetParameter( 0 );

  usr::log::PrintLog( usr::log::INFO,
    "Shifting the data for better plotting" );

  for( int j = 0; j < data->GetN(); ++j ){
    data->GetX()[j] -= zoffset;
    data->GetY()[j] -= ped;
  }

  usr::log::PrintLog( usr::log::INFO,
    "Shifting the function for better plotting" );
  func->SetParameter( 0, 0.0 );
  func->SetRange( usr::plt::GetXmin( data ), usr::plt::GetXmax( data ) );

  usr::log::PrintLog( usr::log::INFO,
    "Refitting to get a correct result container" );
  TFitResult* ans = data->Fit( func, "Q EX0 M E N 0 S" ).Get();

  usr::log::PrintLog( usr::log::INFO,
    "Complete fitting routine" );

  return ans;
}
