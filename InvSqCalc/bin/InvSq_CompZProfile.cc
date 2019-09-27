#include "SiPMCalib/InvSqCalc/interface/LEDFormat.hpp"
#include "SiPMCalib/InvSqCalc/interface/MCFormat.hpp"

#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/PlotUtils/interface/Ratio1DCanvas.hpp"
// #include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

#include <boost/format.hpp>
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

void ShiftData( TGraph* data, double& ped, double& zoffset );

int
main( int argc, char* argv[] )
{
  usr::po::options_description desc( "Options for plot making" );
  desc.add_options()
    ( "data,d", usr::po::value<std::vector<std::string> >()->multitoken(), "Data files " )
    ( "output,o", usr::po::value<std::string>(), "Output file " )
  ;

  usr::ArgumentExtender arg;
  arg.AddOptions( desc );
  arg.ParseOptions( argc, argv );

  std::vector<LEDManager> datalist;

  // Making a list of data
  for( const auto& datafile : arg.ArgList<std::string>( "data" ) ){
    datalist.push_back( LEDManager( datafile ) );
  }


  std::vector<TGraph*> graphlist;
  const double xcentral = datalist.front()._pointlist.front().x;
  const double ycentral = datalist.front()._pointlist.front().y;
  double ped;
  double off;

  for( unsigned i = 0; i < datalist.size(); ++i ){
    const auto& data   = datalist.at( i );
    const double xscan = data._pointlist.front().x;
    const double yscan = data._pointlist.front().y;

    graphlist.push_back( data.MakeZScanGraph( xscan, yscan ) );

    if( i == 0 ){
      ShiftData( graphlist.back(), ped, off );
    } else {
      for( int j = 0; j < graphlist.back()->GetN(); ++j ){
        graphlist.back()->GetY()[j] -= ped;
        //graphlist.back()->GetY()[j] *= std::pow(0.96,i);
        ((TGraphErrors*)(graphlist.back()))->SetPointError(j,0,0);
        graphlist.back()->GetX()[j] -= off;
        // graphlist.back()->GetX()[j] *= std::pow(0.96,i);
      }
    }

  }
  // Getting fit for plotting
  TF1 func = TF1( "func", ZProf,
    usr::plt::GetXmin( graphlist.front() ),
    usr::plt::GetXmax( graphlist.front() ),  4 );
  auto fit   = graphlist.front()->Fit( &func, "EX0 M E N 0 S" );

  usr::plt::Ratio1DCanvas c;

  const std::vector<int> colorlist  = {
     usr::plt::col::blue,
     usr::plt::col::red,
     usr::plt::col::green };
  const std::vector<int> markerlist = {
    usr::plt::sty::mkrcircle,
    usr::plt::sty::mkrdiamond,
    usr::plt::sty::mkrsquare
  };

  auto& fitg = c.PlotFunc( func,
    usr::plt::PlotType( usr::plt::fittedfunc ),
    usr::plt::VisualizeError( fit ),
    usr::plt::EntryText( "Fit (offset=0mm)" ),
    usr::plt::FillColor( kGray, 0.5 ),
    usr::plt::LineColor( kBlack )
     );

  for( int i = graphlist.size() -1; i >= 0; --i ){
    const auto& data   = datalist.at( i );
    const double xscan = data._pointlist.front().x;
    const double yscan = data._pointlist.front().y;
    const double dist = sqrt( ( xscan-xcentral ) * ( xscan-xcentral )
      + ( yscan-ycentral ) * ( yscan-ycentral ) );
    c.PlotGraph( graphlist.at(i),
      usr::plt::PlotType( usr::plt::scatter ),
      usr::plt::TrackY( usr::plt::TrackY( usr::plt::TrackY::both ) ),
      usr::plt::EntryText( usr::fstr( "offset=%.1lf [mm]", dist ) ),
      usr::plt::MarkerColor( colorlist.at(i) ),
      usr::plt::MarkerStyle( markerlist.at(i) ),
      usr::plt::MarkerSize( 0.15 ),
      usr::plt::LineColor( colorlist.at(i), 0.3 )
      );
  }

  // Drawing bottom pad
  c.PlotScale( fitg, fitg, usr::plt::PlotType( usr::plt::fittedfunc ) );

  for( const auto& g : graphlist ){
    c.PlotScale( g, fitg, usr::plt::PlotType( usr::plt::scatter ) );
  }


  c.DrawCMSLabel( "Preliminary", "HGCal" );
  c.DrawLuminosity( "LED Setup" );

  c.TopPad().SetDataMin( 0.5 );
  c.TopPad().SetLogx( 1 );
  c.TopPad().SetLogy( 1 );
  c.BottomPad().SetLogx( 1 );

  c.TopPad().Yaxis().SetTitle( "Luminosity - Ped [pA]" );
  c.BottomPad().Xaxis().SetTitle( "Gantry z - Vert. offset [mm]" );
  c.BottomPad().Yaxis().SetTitle( "Data/Fit" );

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
