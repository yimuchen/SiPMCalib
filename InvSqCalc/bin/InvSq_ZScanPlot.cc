#include "SiPMCalib/InvSqCalc/interface/LEDFormat.hpp"
#include "SiPMCalib/InvSqCalc/interface/MCFormat.hpp"

#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/PlotUtils/interface/Ratio1DCanvas.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

#include "TFitResult.h"

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

void ScaleMC( TGraphAsymmErrors& mc, const TGraph* data, const TF1& datafit );

int
main( int argc, char** argv )
{
  usr::po::options_description desc( "Options for plot making" );
  desc.add_options()
    ( "data,d", usr::po::value<std::string>(), "Data file" )
  //  ( "mc,m", usr::po::value<std::string>(), "MC file" )
    ( "output,o", usr::po::value<std::string>(), "Output file" )
  ;

  usr::ArgumentExtender arg;
  arg.AddOptions( desc );
  arg.ParseOptions( argc, argv );

  LEDManager data( arg.Arg( "data" ) );
  // MCManager mc( arg.Arg( "mc" ) );

  TGraph* dataz = data.MakeZScanGraph(
    data._pointlist.front().x,
    data._pointlist.front().y );

  // TGraphAsymmErrors mc0 = mc.MakeZScanGraph( 1, 1.5, 0 );
  // TGraphAsymmErrors mc5 = mc.MakeZScanGraph( 1, 1.5, 0.5 );

  double pedestal = 0;
  double zoffset  = 0;

  // Fixing data with pedestal and offset
  ShiftData( dataz, pedestal, zoffset );


  // Getting fit for plotting
  TF1 func = TF1( "func", ZProf,
    usr::plt::GetXmin( dataz ), usr::plt::GetXmax( dataz ),
    4 );
  auto fit = dataz->Fit( &func, "EX0 M E N 0 S" );

  // Scaling MC to data
  // ScaleMC( mc0, dataz, func );
  // ScaleMC( mc5, dataz, func );

  // Begin plotting
  usr::plt::Ratio1DCanvas c;

  // auto& mcg = c.PlotGraph( mc0,
  //   usr::plt::PlotType( usr::plt::simplefunc ),
  //   usr::plt::EntryText( "Simulation" ),
  //   usr::plt::TrackY( usr::plt::tracky::max ),
  //    );
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
    usr::plt::MarkerSize( 0.1 ),
    usr::plt::LineColor( usr::plt::col::gray, 0.5 ) );

  c.PlotScale( fitg, fitg,
    usr::plt::PlotType( usr::plt::fittedfunc ) );

  c.PlotScale( datag, fitg,
    usr::plt::PlotType( usr::plt::scatter ) );

  c.TopPad().Yaxis().SetTitle( "Luminosity - Ped. [pA]" );
  c.BottomPad().Xaxis().SetTitle( "Gantry z + z_{0} [mm]" );
  c.BottomPad().Yaxis().SetTitle( "Data,MC/Fit" );

  c.TopPad().SetLogy( kTRUE );
  c.TopPad().SetLogx( kTRUE );
  c.TopPad().SetYaxisMin( dataz->GetY()[dataz->GetN()-1] / 10  );
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

  c.SaveAsPDF( arg.MakePDFFile( arg.Arg( "output" ) ) );

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
ScaleMC( TGraphAsymmErrors& mc, const TGraph* data,  const TF1& func )
{
  // MC scaling
  const double r = ( data->Eval( 100 + func.GetParameter( 0 ) )
                     - func.GetParameter( 3 ) ) / mc.Eval( 100 );

  std::vector<double> z;
  std::vector<double> lumi;
  std::vector<double> lumiup;
  std::vector<double> lumidown;

  for( int i = 0; i < mc.GetN(); ++i ){
    if( mc.GetX()[i] < usr::plt::GetXmin( data )
        || mc.GetX()[i] > usr::plt::GetXmax( data ) ){
      continue;
    }

    z.push_back( mc.GetX()[i] + func.GetParameter( 0 ) );
    lumi.push_back( mc.GetY()[i] * r + func.GetParameter( 3 ) );
    lumiup.push_back( mc.GetErrorYhigh( i ) * r  );
    lumidown.push_back( mc.GetErrorYlow( i ) * r  );
  }

  mc.Set( z.size() );

  for( int i = 0; i < mc.GetN(); ++i ){
    mc.GetX()[i] = z.at( i );
    mc.GetY()[i] = lumi.at( i );
    mc.SetPointEYhigh( i, lumiup.at( i ) );
    mc.SetPointEYlow( i, lumidown.at( i ) );
  }
}
