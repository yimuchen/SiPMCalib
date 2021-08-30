#include "SiPMCalib/InvSqCalc/interface/LEDFormat.hpp"

#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/PlotUtils/interface/Ratio1DCanvas.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

#include "TGraphErrors.h"
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

void FitShiftData( TGraph* data, double& zoffset );
void ShiftData( TGraph* data, const double zoffset  );

int
main( int argc, char* argv[] )
{
  double zoffset;
  TGraph* g4p3 = LEDManager( "LED_4p3V.txt" ).MakeZScanGraph( 10.5, 222.1 );
  TGraph* g3p2 = LEDManager( "LED_3p2V.txt" ).MakeZScanGraph( 10.5, 222.1 );
  TGraph* g2p7 = LEDManager( "LED_2p7V_scaled.txt" ).MakeZScanGraph( 10.5, 222.1 );

  FitShiftData( g4p3, zoffset );
  ShiftData( g3p2, zoffset  );
  ShiftData( g2p7, zoffset  );

  // Getting fit for plotting
  TF1 func = TF1( "func", ZProf,
    usr::plt::GetXmin( g4p3 ),
    usr::plt::GetXmax( g4p3 ),  4 );
  auto fit = g4p3->Fit( &func, "EX0 M E N 0 S" );

  usr::plt::Ratio1DCanvas c;

  c.PlotGraph( g2p7,
    usr::plt::EntryText( "LED@2.7V" ),
    usr::plt::TrackY( usr::plt::tracky::both ),
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::LineColor( usr::plt::col::darkgreen ),
    usr::plt::MarkerColor( usr::plt::col::darkgreen ),
    usr::plt::MarkerStyle( usr::plt::sty::mkrdiamond ),
    usr::plt::MarkerSize( 0.4 )
    );

  c.PlotGraph( g3p2,
    usr::plt::EntryText( "LED@3.2V" ),
    usr::plt::TrackY( usr::plt::tracky::both ),
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::LineColor( usr::plt::col::darkred ),
    usr::plt::MarkerColor( usr::plt::col::darkred ),
    usr::plt::MarkerStyle( usr::plt::sty::mkrsquare ),
    usr::plt::MarkerSize( 0.4 )
    );

  auto& gf = c.PlotFunc( func,
    usr::plt::PlotType( usr::plt::fittedfunc ),
    usr::plt::EntryText( "Fit (4.3V)" ),
    usr::plt::LineColor( usr::plt::col::blue ),
    usr::plt::FillColor( usr::plt::col::skyblue )
    );

  c.PlotGraph( g4p3,
    usr::plt::EntryText( "LED@4.3V" ),
    usr::plt::TrackY( usr::plt::tracky::both ),
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::LineColor( usr::plt::col::black ),
    usr::plt::MarkerColor( usr::plt::col::black ),
    usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
    usr::plt::MarkerSize( 0.4 )
    );

  c.PlotScale( g4p3 , &gf, usr::plt::PlotType( usr::plt::scatter ) );
  c.PlotScale( g3p2 , &gf, usr::plt::PlotType( usr::plt::scatter ) );
  c.PlotScale( g2p7 , &gf, usr::plt::PlotType( usr::plt::scatter ) );

  c.TopPad().SetLogx( kTRUE );
  c.TopPad().SetLogy( kTRUE );

  c.TopPad().Yaxis().SetTitle( "Luminosity - Ped [pA]" );
  c.BottomPad().Xaxis().SetTitle( "Gantry z - Vert. offset [mm]" );

  c.TopPad().DrawLuminosity( "Static LED");
  c.TopPad().DrawCMSLabel( "Luminosity Control", "HGCal" );

  c.SaveAsPDF( "WideZProfile.pdf" );

  return 0;

}

void
FitShiftData( TGraph* data, double& zoffset )
{
  TF1 func = TF1( "shiftfunc", ZProf,  0, 500, 4 );
  TGraph datacopy( *data );
  double ped;

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
ShiftData( TGraph* data, const double zoffset )
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
    if( fabs( func.GetParameter( 3 ) ) < fabs( func.GetParError( 3 ) ) ){
      break;
    }

    for( int j = 0; j < datacopy.GetN(); ++j ){
      datacopy.GetY()[j] -= func.GetParameter( 3 );
    }
  }

  for( int i = 0; i < data->GetN(); ++i ){
    data->GetX()[i] = datacopy.GetX()[i] - zoffset;
    data->GetY()[i] = datacopy.GetY()[i];
  }
}
