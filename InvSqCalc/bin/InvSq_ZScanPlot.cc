#include "SiPMCalib/InvSqCalc/interface/LEDFormat.hpp"
#include "SiPMCalib/InvSqCalc/interface/MCFormat.hpp"

#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/PlotUtils/interface/Ratio1DCanvas.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

#include <boost/format.hpp>

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

void ShiftData( TGraph* data, double& ped, double& zoffset );

void ScaleMC( TGraphAsymmErrors& mc, const TGraph* data );

int
main( int argc, char** argv )
{
  usr::po::options_description desc( "Options for plot making" );
  desc.add_options()
    ( "data,d", usr::po::value<std::string>(), "Data file" )
    ( "mc,m", usr::po::value<std::string>(), "MC file" )
    ( "output,o", usr::po::value<std::string>(), "Output file" )
  ;

  usr::ArgumentExtender arg;
  arg.AddOptions( desc );
  arg.ParseOptions( argc, argv );

  LEDManager data( arg.Arg( "data" ) );
  MCManager mc( arg.Arg( "mc" ) );

  TGraph* dataz = data.MakeZScanGraph(
    data._pointlist.front().x,
    data._pointlist.front().y );

  TGraphAsymmErrors mc0 = mc.MakeZScanGraph( 1, 1.5, 0 );
  TGraphAsymmErrors mc5 = mc.MakeZScanGraph( 1, 1.5, 0.5 );

  double pedestal = 0;
  double zoffset  = 0;

  // Fixing data with pedestal and offset
  ShiftData( dataz, pedestal, zoffset );

  // Scaling MC to data
  ScaleMC( mc0, dataz );
  ScaleMC( mc5, dataz );


  // Getting fit for plotting
  TF1 func = TF1( "func", ZProf,
    usr::plt::GetXmin( dataz ), usr::plt::GetXmax( dataz ),
    4 );
  auto fit = dataz->Fit( &func, "EX0 N 0 S" );


  // Begin plotting
  usr::plt::Ratio1DCanvas c;

  auto& fitg = c.PlotFunc( func,
    usr::plt::PlotType( usr::plt::fittedfunc ),
    usr::plt::ShowFitErr( fit ),
    usr::plt::TrackY( usr::plt::TrackY::max ),
    usr::plt::EntryText( "Fitted Data" ) );
  auto& datag = c.PlotGraph( dataz,
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::TrackY( usr::plt::TrackY::max ),
    usr::plt::EntryText( "LED Readout" ) );
  auto& fitg0 = c.PlotGraph( mc5,
    usr::plt::PlotType( usr::plt::simplefunc ),
    usr::plt::EntryText( "MC (x_{0} = 5mm)" ) );
  auto& fitg5 = c.PlotGraph( mc0,
    usr::plt::PlotType( usr::plt::simplefunc ),
    usr::plt::EntryText( "MC (x_{0} = 0mm)" ) );

  datag.SetLineColor( kGray );
  datag.SetMarkerColor( kBlack );
  datag.SetMarkerStyle( 20 );
  datag.SetMarkerSize( 0.1 );
  fitg.SetLineColor( kBlue );
  fitg.SetFillColorAlpha( kCyan, 0.3 );
  fitg.SetFillStyle( usr::plt::sty::fillsolid );
  fitg0.SetLineColor( kRed );
  fitg5.SetLineColor( kGreen );

  auto& rfit = c.PlotPull( datag, datag,
    usr::plt::PlotType( usr::plt::fittedfunc ) );

  auto& rfit0 = c.PlotPull( fitg0, datag,
    usr::plt::PlotType( usr::plt::simplefunc ) );
  auto& rfit5 = c.PlotPull( fitg5, datag,
    usr::plt::PlotType( usr::plt::simplefunc ) );

  auto& rdata = c.PlotPull( fitg, datag,
    usr::plt::PlotType( usr::plt::fittedfunc ) );

  rfit.SetFillStyle( usr::plt::sty::fillsolid );
  rfit.SetFillColor( 19 );

  rdata.SetFillColorAlpha( kCyan, 0.3 );

  c.TopPad().Yaxis().SetTitle( "Luminosity - Ped. [A.U.]" );
  c.BottomPad().Xaxis().SetTitle( "Gantry z + z_{0} [mm]" );
  c.BottomPad().Yaxis().SetTitle( "Pull (v.s. Data)" );

  c.TopPad().SetLogy( kTRUE );
  c.TopPad().SetLogx( kTRUE );
  c.TopPad().SetYaxisMin( dataz->GetY()[dataz->GetN()-1] / 10  );
  c.BottomPad().SetLogx( kTRUE );

  c.DrawLuminosity( "LED setup" );
  c.DrawCMSLabel( "Preliminary", "HGCal" );
  c.TopPad()
  .WriteLine( "" ).WriteLine( "" ).WriteLine( "" )
  .WriteLine( "" ).WriteLine( "" )
  .WriteLine( "#frac{L_{0}(z+z_{0})}{(x_{0}^{2} + (z+z_{0})^{2})^{3/2}} + Ped" )
  .WriteLine( "" ).WriteLine( "" ).WriteLine( "" )
  .WriteLine( ( boost::format( "Fit Ped = %.2lf_{#pm%.3lf}" )
                % pedestal % func.GetParError( 3 ) ).str() )
  .WriteLine( ( boost::format( "Fit z_{0} = %.2lf_{#pm%.3lf}" )
                % -zoffset % func.GetParError( 0 ) ).str() )
  .WriteLine( ( boost::format( "Fit x_{0} = %.2lf_{#pm%.3lf}" )
                % func.GetParameter( 1 ) % func.GetParError( 0 ) ).str() );

  c.SaveAsPDF( arg.Arg( "output" ) );

  return 0;
}


void
ShiftData( TGraph* data, double& ped, double& zoffset )
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
    std::cout << func.GetParameter( 0 ) << " " << func.GetParError( 0 ) << " "
              << func.GetParameter( 3 ) << " " << func.GetParError( 3 )<< " "
              << std::endl;

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
  std::cout << zoffset << " " << ped << std::endl;

  for( int i = 0; i < data->GetN(); ++i ){
    data->GetX()[i] = datacopy.GetX()[i];
    data->GetY()[i] = datacopy.GetY()[i];
  }
}


void
ScaleMC( TGraphAsymmErrors& mc, const TGraph* data )
{
  // MC scaling
  const double r = data->Eval( 100 ) / mc.Eval( 100 );

  for( int i = 0; i < mc.GetN(); ++i ){
    mc.GetY()[i] *= r;
    mc.SetPointEYhigh( i,
      mc.GetErrorYhigh( i ) * r
      );
    mc.SetPointEYlow( i,
      mc.GetErrorYlow( i ) * r
      );
  }
}
