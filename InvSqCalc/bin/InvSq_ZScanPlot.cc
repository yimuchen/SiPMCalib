/**
 * @file InvSq_ZScanPlot.cc
 * @brief Command `SiPM_InvSq_ZScanPlot`: Running the inverse square fit a data
 * file collected with z scan.
 *
 * @details
 * This function expects the a standard data format, collection via a zscan
 * command with a singular power setting. The data columns are expected to
 * represent the readout value and the uncertainty in the readout value (poisson
 * uncertainty already reduced), this command then performs a inverse square fit
 * to the data, and generates a plot displaying the full output. A summary of the
 * fit results in text format is also displayed on the terminal for rapid
 * debugging.
 *
 * ### List of command specific options:
 *
 * - `--data, -d`: The path to the input data file.
 * - `--output, -o`: The path to store the output plot file (PDF).
 * - `--uncscale,　-u`:　Scale factor to apply to the uncertainty column, this is
 *   in place in case you which to run the fit without the Poisson uncertainty
 *   subtraction, or which to modify the uncertainty readout in any way.
 *
 * ### Additional details:
 *
 * Notice that the plot generated doesn't place the data at the original
 * coordinate position and the origian readout luminosity values. As for the
 * perfect fit with an inverse square, the reader will typically expected that
 * the fit result in the log-log plot is a straight line, which the human eye is
 * very sensitive at picking up deviations from. For this reason, we shift the
 * data to remove the effects of pedestal and z directional offset so that the
 * data presented be a straight line. The removal of these effects are still
 * explicitly spelt out in the axis labels.
 */
#include "SiPMCalib/Common/interface/StdFormat.hpp"

#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/Common/interface/STLUtils/OStreamUtils.hpp"
#include "UserUtils/PlotUtils/interface/Ratio1DCanvas.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

#include "TFitResult.h"
#include "TGraphErrors.h"

/**
 * @brief The inverse square law formula with just the z coordinate as the
 * variable.
 *
 * Notice the value of N0 is used to scale the fit function by some constant
 * amount. As this parameter has perfect correlation with N, this parameter
 * should only ever be set as a constant.
 */
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

TGraphErrors MakeZScanGraph( const StdFormat&, double uncscale );

TFitResult   FitAndShiftData( TGraph& data, TF1& func,
                              double& ped, double& zoffset );

/**
 * @brief Main routine for plotting inverse z scan. Go to the source file for
 * detailed documentation.
 */
int
main( int argc, char** argv )
{
  usr::po::options_description desc( "Options for plot making" );
  desc.add_options()
    ( "data,d", usr::po::reqvalue<std::string>(), "Data file" )
    ( "output,o", usr::po::reqvalue<std::string>(), "Output plot file" )
    ( "uncscale,u", usr::po::defvalue<double>( 1 ), "Additional uncertainty scaling factor" )
  ;

  usr::ArgumentExtender arg;
  arg.AddVerboseOpt();
  arg.AddOptions( desc );
  arg.ParseOptions( argc, argv );

  double pedestal = 0;// For storing the original fit results
  double zoffset  = 0;  // For storing the original fit results

  StdFormat data( arg.Arg( "data" ) );
  TGraph dataz      = MakeZScanGraph( data, arg.Arg<double>( "uncscale" ) );
  const double xmin = usr::plt::GetXmin( dataz );
  const double xmax = usr::plt::GetXmax( dataz );

  // Making the data
  TF1 func( "func", ZProf, xmin, xmax, 5 );
  const TFitResult fit = FitAndShiftData( dataz, func, pedestal, zoffset );
  // Begin plotting
  usr::plt::Ratio1DCanvas c;
  auto& fitg = c.PlotFunc( func,
    usr::plt::PlotType( usr::plt::fittedfunc ),
    usr::plt::EntryText( "Fitted Data" ),
    usr::plt::VisualizeError( fit ),
    usr::plt::FillColor( usr::plt::col::cyan ),
    usr::plt::TrackY( usr::plt::tracky::max ),
    usr::plt::LineColor( usr::plt::col::blue )
    );
  auto& datag = c.PlotGraph( dataz,
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::EntryText( "Readout" ),
    usr::plt::MarkerColor( usr::plt::col::black ),
    usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
    usr::plt::MarkerSize( 0.2 ),
    usr::plt::LineColor( usr::plt::col::gray, 1.0 ) );
  c.PlotScale( fitg, fitg,
    usr::plt::PlotType( usr::plt::fittedfunc ),
    usr::plt::FillColor( usr::plt::col::cyan ) );
  c.PlotScale( datag, fitg,
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::MarkerSize( 0.2 ) );


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
  .WriteLine( usr::fstr( "Fit Ped = %.2lf_{#pm%.3lf}",  pedestal, func.GetParError( 3 ) ) )
  .WriteLine( usr::fstr( "Fit z_{0} = %.2lf_{#pm%.3lf}", -zoffset, func.GetParError( 0 ) ) )
  .WriteLine( usr::fstr( "Fit x_{0} = %.2lf_{#pm%.3lf}", func.GetParameter( 1 ), func.GetParError( 1 ) ) )
  .WriteLine( usr::fstr( "#chi^{2}/D.o.F = %.2lf/%d", fit.Chi2(), fit.Ndf() ) );

  c.SaveAsPDF( arg.Arg( "output" ) );

  return 0;
}

/**
 * @brief Morphing the standard format into a ROOT TGraph container.
 */
TGraphErrors
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

  return TGraphErrors( z.size(),
    z.data(), lumi.data(),
    zero.data(), lerr_scaled.data() );
}

/**
 * @brief Fitting the data to the inverse square function and shift the data
 * points in the data graph to remove pedestal and vertical offset effect.
 *
 * The fit will be performed a first time to determined the pedestal and vertical
 * off set. Then, the data will be shifted according to the fitted pedestal and
 * vertical positions. Here we also scale the luminosity normalization such that
 * the fit uncertainty on the luminosity should be the same as the vertical
 * offset uncertainty. The pedestal and vertical offset will be stored in the
 * latter input inputs.
 *
 * Once the data shifting is complete, we then perform the fit again such that
 * the function parameters properly match the data.
 */
TFitResult
FitAndShiftData( TGraph& data, TF1& func, double& ped, double& zoffset )
{
  usr::log::PrintLog( usr::log::INFO,
    "Getting the estimate values for fitting" );
  const double lumimax = usr::plt::GetYmax( data );
  const double lumimin = usr::plt::GetYmin( data );
  const double zmin    = usr::plt::GetXmin( data );

  func.SetParameter( 0, 5 );
  func.SetParameter( 1, 0 );
  func.SetParameter( 2, ( lumimax*( zmin+5 )*( zmin+5 ) ) );
  func.SetParameter( 3, lumimin );
  func.FixParameter( 4, 1.0 );

  usr::log::PrintLog( usr::log::INFO, "Running the fit" );
  data.Fit( &func, "Q EX0 M E N 0" ).Get();

  // Saving the original fit results for the pedestal and zoffset
  ped     = func.GetParameter( 3 );
  zoffset = func.GetParameter( 0 );
  const double scale = func.GetParError( 0 ) / func.GetParError( 2 );
  func.FixParameter( 4, 1/ scale  );
  func.SetParameter( 2, func.GetParameter( 2 ) *  scale );

  usr::log::PrintLog( usr::log::INFO,
    "Shifting the data for better plotting" );

  for( int j = 0; j < data.GetN(); ++j ){
    data.GetX()[j] -= zoffset;
    data.GetY()[j] -= ped;
  }

  usr::log::PrintLog( usr::log::INFO,
    "Shifting the function for better plotting" );
  func.SetParameter( 0, 0.0 );
  func.SetRange( usr::plt::GetXmin( data ), usr::plt::GetXmax( data ) );

  usr::log::PrintLog( usr::log::INFO,
    "Refitting to get a correct result container" );
  TFitResult* ans = data.Fit( &func, "Q EX0 M E N 0 S" ).Get();

  usr::log::PrintLog( usr::log::INFO,
    "Complete fitting routine" );

  usr::log::PrintLog( usr::log::INFO,
    usr::fstr( "Fitted z0: %.2lf +- %.2lf", zoffset, func.GetParError( 0 ) ) );
  usr::log::PrintLog( usr::log::INFO,
    usr::fstr( "Fitted x0: %.2lf +- %.2lf", func.GetParameter( 1 ), func.GetParError( 1 ) ) );
  usr::log::PrintLog( usr::log::INFO,
    usr::fstr( "Fitted L0: %.2lf +- %.2lf", func.GetParameter( 2 ), func.GetParError( 2 ) ) );
  usr::log::PrintLog( usr::log::INFO,
    usr::fstr( "Fitted ped: %.2lf +- %.2lf", ped, func.GetParError( 3 ) ) );

  ans->Print();

  return *ans;
}
