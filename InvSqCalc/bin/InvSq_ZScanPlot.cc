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

  const double D2 = o*o + ( z-z0 )*( z-z0 );
  // return N/D2 + P;
  return ( N*( z-z0 ) )/( D2 * sqrt( D2 ) ) + P;
}

TGraph* MakeZScanGraph( const StdFormat&, double uncscale );
void    FitAndShiftData( TGraph* data, double& ped, double& zoffset );

int
main( int argc, char** argv )
{
  usr::po::options_description desc( "Options for plot making" );
  desc.add_options()
    ( "data,d", usr::po::value<std::string>(), "Data file" )
    ( "output,o", usr::po::value<std::string>(), "Output plot file" )
    ( "uncscale,u", usr::po::defvalue<double>( 1 ),
    "Additional uncertainty scaling factor" )
  ;

  usr::ArgumentExtender arg;
  arg.AddOptions( desc );
  arg.ParseOptions( argc, argv );

  double pedestal = 0;
  double zoffset  = 0;

  StdFormat data( arg.Arg( "data" ) );
  TGraph* dataz = MakeZScanGraph( data, arg.Arg<double>( "uncscale" ) );
  FitAndShiftData( dataz, pedestal, zoffset );

  TF1 func = TF1( "func", ZProf,
    usr::plt::GetXmin( dataz ), usr::plt::GetXmax( dataz ),
    4 );
  auto fit = dataz->Fit( &func, "EX0 M E N 0 S" );
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
    usr::plt::PlotType( usr::plt::fittedfunc ) );
  c.PlotScale( datag, fitg,
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::MarkerSize( 0.3 ) );

  c.TopPad().Yaxis().SetTitle( "Luminosity - Ped. [pA]" );
  c.BottomPad().Xaxis().SetTitle( "Gantry z + z_{0} [mm]" );
  c.BottomPad().Yaxis().SetTitle( "Data,MC/Fit" );

  c.TopPad().SetLogy( kTRUE );
  c.TopPad().SetLogx( kTRUE );
  c.TopPad().SetYaxisMin( dataz->GetY()[dataz->GetN()-1] / 30  );
  c.BottomPad().SetLogx( kTRUE );

  c.DrawLuminosity( "LED setup" );
  c.DrawCMSLabel( "Preliminary", "HGCal" );
  c.TopPad()
  .SetTextCursor( c.TopPad().InnerTextLeft(), 0.45 )
  .WriteLine( "#frac{L_{0}(z+z_{0})}{(x_{0}^{2} + (z+z_{0})^{2})^{3/2}} + Ped" )
  // .WriteLine( "" )
  .WriteLine( ( boost::format( "Fit Ped = %.2lf_{#pm%.3lf}" )
                % pedestal % func.GetParError( 3 ) ).str() )
  .WriteLine( ( boost::format( "Fit z_{0} = %.2lf_{#pm%.3lf}" )
                % -zoffset % func.GetParError( 0 ) ).str() )
  .WriteLine( ( boost::format( "Fit x_{0} = %.2lf_{#pm%.3lf}" )
                % func.GetParameter( 1 ) % func.GetParError( 0 ) ).str() )
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

  std::cout << lumierr << std::endl;
  std::cout << lerr_scaled << std::endl;

  return new TGraphErrors( z.size(),
    z.data(), lumi.data(),
    zero.data(), lerr_scaled.data() );
}


void
FitAndShiftData( TGraph* data, double& ped, double& zoffset )
{
  TF1 func = TF1( "shiftfunc", ZProf,  0, 500, 4 );
  TGraph datacopy( *data );

  for( unsigned i = 0; i < 100; ++i ){
    const double lumimax = usr::plt::GetYmax( datacopy );
    const double lumimin = usr::plt::GetYmin( datacopy );
    const double zmin    = usr::plt::GetXmin( datacopy );

    func.SetParameters(
      0.0, 0.0,
      lumimax * zmin * zmin - lumimin,
      lumimin
      );

    datacopy.Fit( &func, "Q EX0 N 0" );

    if( fabs( func.GetParameter( 0 ) ) < fabs( func.GetParError( 0 ) ) &&
        fabs( func.GetParameter( 3 ) ) < fabs( func.GetParError( 3 ) ) ){
      break;
    }

    for( int j = 0; j < datacopy.GetN(); ++j ){
      datacopy.GetX()[j] -= func.GetParameter( 0 );
      datacopy.GetY()[j] -= func.GetParameter( 3 );
    }
  }

  zoffset = data->GetX()[0] - datacopy.GetX()[0];
  ped     = data->GetY()[0] - datacopy.GetY()[0];

  for( int i = 0; i < data->GetN(); ++i ){
    data->GetX()[i] = datacopy.GetX()[i];
    data->GetY()[i] = datacopy.GetY()[i];
  }
}
