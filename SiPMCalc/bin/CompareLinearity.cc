#include "SiPMCalib/SiPMCalc/interface/NonLinearModel.hpp"

#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/Common/interface/STLUtils/OStreamUtils.hpp"
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"
#include "UserUtils/MathUtils/interface/Measurement/Measurement.hpp"
#include "UserUtils/PlotUtils/interface/Ratio1DCanvas.hpp"

#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TStyle.h"

#include <boost/format.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

struct LumiPoint
{
  unsigned Aidx;
  unsigned Bidx;
  double   lumival;
  double   lumiunc;
};


TGraph* MakeLinearityGraph( const std::string&,
                            std::string&,
                            std::string&,
                            std::string& );

int
main( int argc, char* argv[] )
{
  usr::po::options_description desc(
    "Options for fitting nonlinear readout from Filter wheel data" );
  desc.add_options()
    ( "data", usr::po::multivalue<std::string>(),
    "Input data .txt files" )
    ( "output", usr::po::value<std::string>(),
    "Primary output file name" )
  ;

  usr::ArgumentExtender arg;
  arg.AddOptions( desc );
  arg.ParseOptions( argc, argv );

  const std::vector<std::string> infile_list = arg.ArgList<std::string>( "data" );
  const std::string outputfile               = arg.Arg<std::string>( "output" );

  std::vector<TGraph*> graphlist;
  std::vector<std::string> entrylist;
  std::string model = "";
  std::string bias  = "";
  double xmin       = 100;
  double xmax       = 0;

  for( const auto infile : infile_list ){
    std::string entry;
    TGraph* graph = MakeLinearityGraph( infile, model, bias, entry );
    graphlist.push_back( graph );
    entrylist.push_back( entry );
    xmin = std::min( usr::plt::GetXmin( graph ), xmin );
    xmax = std::max( usr::plt::GetXmax( graph ), xmax );
  }

  TGraph dummy( 2 );
  dummy.SetPoint( 0, xmin, 1 );
  dummy.SetPoint( 1, xmax, 1 );


  const std::vector<int> colorlist = {
    usr::plt::col::black,
    usr::plt::col::darkblue,
    usr::plt::col::darkred,
    usr::plt::col::darkgreen
  };


  usr::plt::Ratio1DCanvas c;
  unsigned i = 0;

  c.PlotGraph( dummy,
    usr::plt::PlotType( usr::plt::simplefunc ),
    usr::plt::LineColor( usr::plt::col::white, 0 ) );

  for( const auto g : graphlist ){
    c.PlotGraph( g,
      usr::plt::EntryText( entrylist.at( i ) ),
      usr::plt::PlotType( usr::plt::scatter ),
      usr::plt::TrackY( usr::plt::tracky::both ),
      usr::plt::LineColor( colorlist[i] ),
      usr::plt::MarkerColor( colorlist[i] ),
      usr::plt::MarkerSize( 0.2 ),
      usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
      usr::plt::FillStyle( usr::plt::col::lightgray )
      );
    ++i;
  }

  c.BottomPad().PlotGraph( dummy );
  c.PlotScale( graphlist[0], graphlist[0],
    usr::plt::PlotType( usr::plt::fittedfunc ),
    usr::plt::FillStyle( usr::plt::sty::fillsolid ),
    usr::plt::FillColor( usr::plt::col::lightgray )
    );

  for( unsigned i = 1; i < graphlist.size(); ++i  ){
    c.PlotScale( graphlist[i], graphlist[0],
      usr::plt::PlotType( usr::plt::scatter ) );
  }

  c.TopPad().DrawVLine( usr::plt::GetXmax( graphlist[0] ),
    usr::plt::LineColor( usr::plt::col::gray ),
    usr::plt::LineStyle( usr::plt::sty::lindensedot )
    );
  c.BottomPad().DrawVLine( usr::plt::GetXmax( graphlist[0] ),
    usr::plt::LineColor( usr::plt::col::gray ),
    usr::plt::LineStyle( usr::plt::sty::lindensedot )
    );
  c.BottomPad().DrawHLine( 1.1,
    usr::plt::LineColor( usr::plt::col::gray ),
    usr::plt::LineStyle( usr::plt::sty::lindashed ) );
  c.BottomPad().DrawHLine( 0.9,
    usr::plt::LineColor( usr::plt::col::gray ),
    usr::plt::LineStyle( usr::plt::sty::lindashed ) );


  c.TopPad().SetLogx( 1 );
  c.BottomPad().SetLogx( 1 );
  c.TopPad().SetLogy( 1 );
  c.BottomPad().Xaxis().SetTitle( "#bar{N}(p.e.) / N_{pix}" );
  c.BottomPad().Yaxis().SetTitle( "Readout / Laser" );
  c.TopPad().Yaxis().SetTitle( "Readout [V-#mu s]" );
  c.TopPad().FinalizeLegend( usr::plt::align::bottom_right );

  c.DrawLuminosity( "Filter Wheel Setup" );
  c.DrawCMSLabel( "", "Linearity Test" );
  c.TopPad().WriteLine( model );
  c.TopPad().WriteLine( bias );

  c.SaveAsPDF( outputfile );

  return 0;

}

TGraph*
MakeLinearityGraph( const std::string& input,
                    std::string&       model,
                    std::string&       bias,
                    std::string&       entry_text )
{
  std::vector<LumiPoint> raw;
  std::string line;
  double N;
  double BV;
  unsigned centralA;
  unsigned centralB;
  double en0;

  std::ifstream inputfile( input );
  std::stringstream ss;
  std::getline( inputfile, entry_text );
  std::getline( inputfile, line );
  ss.str( line );
  ss >> model >> N >> BV;
  std::getline( inputfile, line );
  ss.clear();
  ss.str( line );
  ss >> centralA >> centralB >> en0;
  en0 /= A( centralA ).CentralValue() * B( centralB ).CentralValue();


  while( std::getline( inputfile, line ) ){
    unsigned a, b;
    double lumi, lumiunc;
    std::stringstream ss( line );
    ss >> a >> b >> lumi >> lumiunc;
    raw.push_back( {a, b, lumi, lumiunc} );
  }

  std::sort( raw.begin(), raw.end(),
    []( const LumiPoint& l, const LumiPoint& r ) -> bool {
    return ( A( l.Aidx ).CentralValue()*B( l.Bidx ).CentralValue() )
    < ( A( r.Aidx ).CentralValue()*B( r.Bidx ).CentralValue() );
  } );


  TGraphErrors* graph = new TGraphErrors( raw.size() );

  for( unsigned i = 0; i < raw.size(); ++i ){
    const unsigned idx1 = raw[i].Aidx;
    const unsigned idx2 = raw[i].Bidx;
    const double xval   = ( A( idx1 ).CentralValue()*B( idx2 ).CentralValue() );
    const double xunc   = ( A( idx1 )*B( idx2 ) ).AbsAvgError();
    // Scaling of luminosity is required to keep fit in a reasonable value
    const double lumi = raw[i].lumival / 1e3;
    const double lunc = raw[i].lumiunc / 1e3;
    graph->SetPoint( i, en0*xval /N, lumi );
    graph->SetPointError( i, en0*xunc/ N, lunc );
  }

  bias       = usr::fstr( "V_{bias} = %.1lf", BV );
  model      = usr::fstr( "%s(%d pix)", model, N );
  entry_text = usr::fstr( "%s", entry_text );

  return graph;
}
