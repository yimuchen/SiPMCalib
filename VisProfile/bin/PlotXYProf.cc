#include <UserUtils/PlotUtils/interface/Flat2DCanvas.hpp>

#include "TF2.h"
#include "TFitResult.h"
#include "TGraph2DErrors.h"

#include <boost/format.hpp>
#include <fstream>
#include <sstream>
#include <vector>

TH2D* MakeVisXGraph( const std::string& );
TH2D* MakeVisYGraph( const std::string& );

double
Func( const double* xy, const double* param )
{
  const double x = xy[0];
  const double y = xy[1];
  const double a = param[0];
  const double b = param[1];
  const double c = param[2];

  return a*x + b*y + c;
}

int
main( int argc, char const* argv[] )
{
  {
    auto xprof = MakeVisXGraph( argv[1] );

    TF2 func( "func", Func,
              xprof->GetXaxis()->GetXmin(), xprof->GetXaxis()->GetXmax(),
              xprof->GetYaxis()->GetXmin(), xprof->GetYaxis()->GetXmax(), 3 );
    auto fit = xprof->Fit( &func, "EX0 N 0 S" );

    usr::plt::Flat2DCanvas c;


    auto& xgraph = c.PlotHist( xprof,
      usr::plt::Plot2DF( usr::plt::heat ),
      usr::plt::EntryText( "Camera readout" ) );
    auto& fgraph = c.PlotFunc( func,
      usr::plt::Plot2DF( usr::plt::cont ),
      usr::plt::EntryText( "Linear fit" ) );

    xgraph.SetLineColor( kBlack );
    fgraph.SetLineColor( kRed );

    c.DrawCMSLabel( "Preliminary", "HGCal" );
    c.DrawLuminosity( "Open light" );
    c.Pad().SetTextCursor( 0.015, 0.6, usr::plt::font::top_left )
    .WriteLine( "Cam X = a x + b y + C" )
    .WriteLine(
      ( boost::format( "a = %.1lf_{#pm%.2lf} [pix/mm]" )
        % func.GetParameter( 0 )
        % func.GetParError( 0 ) ).str() )
    .WriteLine(
      ( boost::format( "b = %.1lf_{#pm%.2lf} [pix/mm]" )
        % func.GetParameter( 1 )
        % func.GetParError( 1 ) ).str() )
    .WriteLine(
      ( boost::format( "C = %.1lf_{#pm%.2lf} [pix]" )
        % func.GetParameter( 2 )
        % func.GetParError( 2 )
      ).str() )
    .WriteLine( ( boost::format( "#chi^{2}/DoF = %.3lf" )% ( fit->Chi2()/fit->Ndf() ) ).str() );


    c.Xaxis().SetTitle( "Gantry X [mm]" );
    c.Yaxis().SetTitle( "Gantry Y [mm]" );
    c.Zaxis().SetTitle( "Camera chip X position [pixels]" );

    c.SaveAsPDF( "xprofile.pdf" );
  }
  {
    auto xprof = MakeVisYGraph( argv[1] );

    TF2 func( "func", Func,
              xprof->GetXaxis()->GetXmin(), xprof->GetXaxis()->GetXmax(),
              xprof->GetYaxis()->GetXmin(), xprof->GetYaxis()->GetXmax(), 3 );
    auto fit = xprof->Fit( &func, "EX0 N 0 S" );

    usr::plt::Flat2DCanvas c;


    auto& xgraph = c.PlotHist( xprof,
      usr::plt::Plot2DF( usr::plt::heat ),
      usr::plt::EntryText( "Camera readout" ) );
    auto& fgraph = c.PlotFunc( func,
      usr::plt::Plot2DF( usr::plt::cont ),
      usr::plt::EntryText( "Linear fit" ) );

    xgraph.SetLineColor( kBlack );
    fgraph.SetLineColor( kRed );

    c.DrawCMSLabel( "Preliminary", "HGCal" );
    c.DrawLuminosity( "Open light" );
    c.Pad().SetTextCursor( 0.015, 0.6, usr::plt::font::top_left )
    .WriteLine( "Cam Y = a x + b y + C" )
    .WriteLine(
      ( boost::format( "a = %.1lf_{#pm%.2lf} [pix/mm]" )
        % func.GetParameter( 0 )
        % func.GetParError( 0 ) ).str() )
    .WriteLine(
      ( boost::format( "b = %.1lf_{#pm%.2lf} [pix/mm]" )
        % func.GetParameter( 1 )
        % func.GetParError( 1 ) ).str() )
    .WriteLine(
      ( boost::format( "C = %.1lf_{#pm%.2lf} [pix]" )
        % func.GetParameter( 2 )
        % func.GetParError( 2 )
      ).str() )
    .WriteLine( ( boost::format( "#chi^{2}/DoF = %.3lf" )% ( fit->Chi2()/fit->Ndf() ) ).str() );


    c.Xaxis().SetTitle( "Gantry X [mm]" );
    c.Yaxis().SetTitle( "Gantry Y [mm]" );
    c.Zaxis().SetTitle( "Camera chip Y position [pixels]" );

    c.SaveAsPDF( "yprofile.pdf" );
  }
  return 0;
}


TH2D*
MakeVisXGraph( const std::string& file )
{
  std::ifstream fin;
  std::string line;
  fin.open( file, std::ifstream::in );

  std::vector<double> xlist;
  std::vector<double> ylist;
  std::vector<double> visxlist;

  double x, y, z, visx, visy;

  while( std::getline( fin, line ) ){
    std::istringstream linestream( line );
    linestream >> x >> y >> z >> visx >> visy;

    if( visx < 0 || visy < 0 ){ continue; }

    xlist.push_back( x );
    ylist.push_back( y );
    visxlist.push_back( visx );
  }

  const double diff = std::max(
    fabs( xlist.at( 0 )- xlist.at( 1 ) ),
    fabs( ylist.at( 0 )- ylist.at( 1 ) ) );

  const double xmax = *std::max_element( xlist.begin(), xlist.end() );
  const double xmin = *std::min_element( xlist.begin(), xlist.end() );
  const double ymin = *std::min_element( ylist.begin(), ylist.end() );
  const double ymax = *std::max_element( ylist.begin(), ylist.end() );

  TH2D* ans = new TH2D( ( "hist"+usr::RandomString( 6 ) ).c_str(), "",
    ( xmax-xmin )/diff + 1, xmin-0.5*diff, xmax+0.5*diff,
    ( ymax-ymin )/diff + 1, ymin-0.5*diff, ymax+0.5*diff
     );

  for( unsigned i = 0; i < xlist.size(); ++i ){
    const int binidx = ans->FindBin( xlist.at( i ), ylist.at( i ) );
    ans->SetBinContent( binidx, visxlist.at( i ) );
    ans->SetBinError( binidx, 0.5 );
  }

  ans->SetStats( 0 );

  return ans;
}

TH2D*
MakeVisYGraph( const std::string& file )
{
  std::ifstream fin;
  std::string line;
  fin.open( file, std::ifstream::in );

  std::vector<double> xlist;
  std::vector<double> ylist;
  std::vector<double> visylist;

  double x, y, z, visx, visy;

  while( std::getline( fin, line ) ){
    std::istringstream linestream( line );
    linestream >> x >> y >> z >> visx >> visy;

    if( visx < 0 || visy < 0 ){ continue; }

    xlist.push_back( x );
    ylist.push_back( y );
    visylist.push_back( visy );
  }

  const double diff = std::max(
    fabs( xlist.at( 0 )- xlist.at( 1 ) ),
    fabs( ylist.at( 0 )- ylist.at( 1 ) ) );

  const double xmax = *std::max_element( xlist.begin(), xlist.end() );
  const double xmin = *std::min_element( xlist.begin(), xlist.end() );
  const double ymin = *std::min_element( ylist.begin(), ylist.end() );
  const double ymax = *std::max_element( ylist.begin(), ylist.end() );

  TH2D* ans = new TH2D( ( "hist"+usr::RandomString( 6 ) ).c_str(), "",
    ( xmax-xmin )/diff + 1, xmin-0.5*diff, xmax+0.5*diff,
    ( ymax-ymin )/diff + 1, ymin-0.5*diff, ymax+0.5*diff
     );

  for( unsigned i = 0; i < xlist.size(); ++i ){
    const int binidx = ans->FindBin( xlist.at( i ), ylist.at( i ) );
    ans->SetBinContent( binidx, visylist.at( i ) );
    ans->SetBinError( binidx, 0.5 );
  }

  ans->SetStats( 0 );

  return ans;
}

