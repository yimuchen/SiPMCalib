#include <UserUtils/PlotUtils/interface/Ratio1DCanvas.hpp>
#include <UserUtils/PlotUtils/interface/Simple1DCanvas.hpp>

#include "TF1.h"

#include <fstream>
#include <sstream>
#include <vector>

struct ZprofileGraphs
{
  TGraph* sharpness;
  TGraph* recox;
  TGraph* recoy;
  TGraph* recoarea;
  TGraph* recod;
};

ZprofileGraphs MakeProfileGraph( const std::string& );

int
main( int argc, char const* argv[] )
{
  auto zprof = MakeProfileGraph( argv[1] );

  // TF1 func( "func", "[0]*(x-[1])*(x-[1]) + [2]", 35, 42 );
  // func.SetParameters( usr::plt::GetYmax( zprof ), 40, 0 );
  // zprof->Fit( &func, "\1\2 0 R" );

  {// Sharpness graphs
    usr::plt::Simple1DCanvas c;

    auto& graph = c.PlotGraph( zprof.sharpness,
      usr::plt::PlotType( usr::plt::scatter ),
      usr::plt::EntryText( "Camera readout" ),
      usr::plt::TrackY( usr::plt::TrackY::both ) );

    graph.SetMarkerStyle( 20 );
    graph.SetMarkerSize( 0.2 );
    graph.SetMarkerColor( kRed );

    c.Xaxis().SetTitle( "Grantry z [mm]" );
    c.Yaxis().SetTitle( "Sharpness measure [A.U.]" );

    c.DrawCMSLabel( "Preliminary", "HGCal" );
    c.DrawLuminosity( "Computer vision" );

    c.SaveAsPDF( "zscan_sharpness.pdf" );
  }
  {// Drift plot
    usr::plt::Simple1DCanvas c;
    auto& xgraph = c.PlotGraph( zprof.recox,
      usr::plt::PlotType( usr::plt::scatter ),
      usr::plt::EntryText( "FOV X" ),
      usr::plt::TrackY( usr::plt::TrackY::both ) );
    auto& ygraph = c.PlotGraph( zprof.recoy,
      usr::plt::PlotType( usr::plt::scatter ),
      usr::plt::EntryText( "FOV Y" ),
      usr::plt::TrackY( usr::plt::TrackY::both ) );

    xgraph.SetMarkerStyle( 20 );
    xgraph.SetMarkerSize( 0.2 );
    xgraph.SetMarkerColor( kRed );
    xgraph.SetLineColor( kRed );

    ygraph.SetMarkerStyle( 21 );
    ygraph.SetMarkerSize( 0.2 );
    ygraph.SetMarkerColor( kBlue );
    ygraph.SetLineColor( kBlue );

    // c.DrawHLine(  7, kGray, 2 );
    // c.DrawHLine( -7, kGray, 2 );

    c.Xaxis().SetTitle( "Gantry Z [mm]" );
    c.Yaxis().SetTitle( "Drift from Center [pix]" );

    c.DrawCMSLabel( "Preliminary", "HGCal" );
    c.DrawLuminosity( "Computer Vision" );

    c.SaveAsPDF( "zscan_xydrift.pdf" );
  }
  {// Reconstructe\1\2re\1\2lot
    TF1 func( "func", "[0]/((x-[1])*(x-[1]))",
              usr::plt::GetXmin( zprof.recoarea ),
              usr::plt::GetXmax( zprof.recoarea ) );
    func.SetParameters( usr::plt::GetYmax( zprof.recoarea ), 10 );
    zprof.recoarea->Fit( &func, "W N 0 R" );

    usr::plt::Ratio1DCanvas c;
    auto& graph = c.PlotGraph( zprof.recoarea,
      usr::plt::PlotType( usr::plt::scatter ),
      usr::plt::EntryText( "Camera readout" ),
      usr::plt::TrackY( usr::plt::TrackY::both ) );
    auto& fgraph = c.PlotFunc( func,
      usr::plt::PlotType( usr::plt::simplefunc ),
      usr::plt::EntryText( "Fit (inv. sq.)" ),
      usr::plt::TrackY( usr::plt::TrackY::both ) );

    graph.SetMarkerStyle( 20 );
    graph.SetMarkerSize( 0.2 );
    graph.SetMarkerColor( kRed );
    graph.SetLineColor( kRed );

    fgraph.SetLineColor( kBlue );

    c.PlotScale( graph,  fgraph, usr::plt::PlotType( usr::plt::scatter ) );
    c.PlotScale( fgraph, fgraph );

    c.DrawCMSLabel( "Preliminary", "HGCal" );
    c.DrawLuminosity( "Computer vision" );

    c.TopPad().Yaxis().SetTitle( "Reconstructed area [Num. Pixels]" );
    c.BottomPad().Xaxis().SetTitle( "Gantry Z [mm]" );
    c.BottomPad().Yaxis().SetTitle( "Data/Fit" );
    c.SaveAsPDF( "zscan_area.pdf" );

  }
  {// Max measure plot
    TF1 func( "func2", "[0]/((x-[1]))",
              usr::plt::GetXmin( zprof.recoarea ),
              usr::plt::GetXmax( zprof.recoarea ) );
    func.SetParameters( usr::plt::GetYmax( zprof.recod ), 10 );
    zprof.recod->Fit( &func, "W N 0 R" );

    usr::plt::Ratio1DCanvas c;
    auto& graph = c.PlotGraph( zprof.recod,
      usr::plt::PlotType( usr::plt::scatter ),
      usr::plt::EntryText( "Camera readout" ),
      usr::plt::TrackY( usr::plt::TrackY::both ) );
    auto& fgraph = c.PlotFunc( func,
      usr::plt::PlotType( usr::plt::simplefunc ),
      usr::plt::EntryText( "Fit (N/(z+z_{0}))" ),
      usr::plt::TrackY( usr::plt::TrackY::both ) );

    graph.SetMarkerStyle( 20 );
    graph.SetMarkerSize( 0.2 );
    graph.SetMarkerColor( kRed );
    graph.SetLineColor( kRed );

    fgraph.SetLineColor( kBlue );

    c.PlotScale( graph,  fgraph, usr::plt::PlotType( usr::plt::scatter ) );
    c.PlotScale( fgraph, fgraph );

    c.DrawCMSLabel( "Preliminary", "HGCal" );
    c.DrawLuminosity( "Computer vision" );

    c.TopPad().Yaxis().SetTitle( "Max Measure [Pixels]" );
    c.BottomPad().Xaxis().SetTitle( "Gantry Z [mm]" );
    c.BottomPad().Yaxis().SetTitle( "Data/Fit" );
    c.SaveAsPDF( "zscan_measure.pdf" );
  }

  return 0;
}


ZprofileGraphs
MakeProfileGraph( const std::string& file )
{
  std::ifstream fin;
  std::string line;
  fin.open( file, std::ifstream::in );

  std::vector<double> zlist;
  std::vector<double> recozlist;
  std::vector<double> slist;
  std::vector<double> recoxlist;
  std::vector<double> recoylist;
  std::vector<double> recoalist;
  std::vector<double> recodlist;

  double x, y, z, sharp, recox, recoy, recoa, recod;

  while( std::getline( fin, line ) ){
    std::istringstream linestream( line );
    linestream >> x >> y >> z >> sharp >> recox >> recoy >> recoa >> recod;

    zlist.push_back( z );
    slist.push_back( sharp );
    if( recox > 0 ){
      recozlist.push_back( z );
      recoxlist.push_back( recox - 640 );
      recoylist.push_back( recoy - 512 );
      recoalist.push_back( recoa );
      recodlist.push_back( recod );
    }
  }

  return ZprofileGraphs {
    new TGraph( zlist.size(), zlist.data(), slist.data() ),
    new TGraph( recozlist.size(), recozlist.data(), recoxlist.data() ),
    new TGraph( recozlist.size(), recozlist.data(), recoylist.data() ),
    new TGraph( recozlist.size(), recozlist.data(), recoalist.data() ),
    new TGraph( recozlist.size(), recozlist.data(), recodlist.data() )
  };
}
