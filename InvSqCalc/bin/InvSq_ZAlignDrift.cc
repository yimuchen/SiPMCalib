#include "SiPMCalib/InvSqCalc/interface/LEDFormat.hpp"
#include "SiPMCalib/InvSqCalc/interface/MCFormat.hpp"

#include "UserUtils/Common/interface/STLUtils/OStreamUtils.hpp"
#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

#include "TF2.h"
#include "TGraphErrors.h"

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

  std::vector<double> zlist;
  std::vector<double> zero;
  std::vector<double> xdrift;
  std::vector<double> xdrifterr;
  std::vector<double> ydrift;
  std::vector<double> ydrifterr;

  for( const auto& data : datalist ){
    const double z = data._pointlist.front().z;
    auto* lumiscan = data.MakeHScanGraph( z );
    TF2 func( "func", ExpFunc,
              data.Xmin(), data.Xmax(), data.Ymin(), data.Ymax(), 5 );
    func.SetParameters(
      ( data.Xmin() + data.Xmax() )/2,
      ( data.Ymin() + data.Ymax() )/2,
      z,
      ( data.LumiMax() - data.LumiMin() )* ( z*z ),
      data.LumiMin()
      );

    lumiscan->Fit( &func, "EX0 N 0" );

    zlist.push_back( z );
    zero.push_back( 0 );
    xdrift.push_back( func.GetParameter( 0 ) );
    xdrifterr.push_back( func.GetParError( 0 ) );
    ydrift.push_back( func.GetParameter( 1 ) );
    ydrifterr.push_back( func.GetParError( 1 ) );

    usr::fout( "%lf %lf %lf %lf %lf\n"
      , zlist.back()
      , xdrift.back()
      , xdrifterr.back()
      , ydrift.back()
      , ydrifterr.back() );

    delete lumiscan;
  }

  // Shifting to zero
  const double xcen = xdrift.front();
  const double ycen = ydrift.front();

  for( unsigned i = 0; i < xdrift.size(); ++i ){
    xdrift[i] -= xcen;
    ydrift[i] -= ycen;
    std::cout << ydrift.at(i) << std::endl;
  }

  TGraphErrors xdriftgraph( zlist.size(),
                            zlist.data(), xdrift.data(), zero.data(), xdrifterr.data() );
  TGraphErrors ydriftgraph( zlist.size(),
                            zlist.data(), ydrift.data(), zero.data(), ydrifterr.data() );

  usr::plt::Simple1DCanvas c;

  c.PlotGraph( xdriftgraph,
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::EntryText( "Lumi-aligned X" ),
    usr::plt::TrackY( usr::plt::TrackY::both )
    );
  c.PlotGraph( ydriftgraph,
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::EntryText( "Lumi-aligned Y" ),
    usr::plt::TrackY( usr::plt::TrackY::both )
    );


  xdriftgraph.SetLineColor( kBlue );
  xdriftgraph.SetMarkerColor( kBlue );
  xdriftgraph.SetMarkerStyle( usr::plt::sty::mkrcircle );
  xdriftgraph.SetMarkerSize( 0.5 );

  ydriftgraph.SetLineColor( kRed );
  ydriftgraph.SetMarkerColor( kRed );
  ydriftgraph.SetMarkerStyle( usr::plt::sty::mkrsquare );
  ydriftgraph.SetMarkerSize( 0.5 );

  c.DrawCMSLabel( "Preliminary", "HGCal" );
  c.DrawLuminosity( "LED Setup" );

  c.Xaxis().SetTitle( "Scan gantry z [mm]" );
  c.Yaxis().SetTitle( "Horizontal drift [mm]" );
  // c.Zaxis().SetTitle("Lumin")

  c.SaveAsPDF( arg.Arg( "output" ) );

  return 0;
}
