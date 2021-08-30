#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

#include "TSpectrum.h"

#include <vector>

int
main( int argc, char* argv[] )
{
  usr::po::options_description desc( "Plotting 1D stability profiles\n"
                                     "Leave one of biasrange, temppulser, or tempsipm empty for it to be used"
                                     " as the x axis" );
  desc.add_options()
    ( "data", usr::po::value<std::string>(),
    "Input data .txt file" )
    ( "outpre", usr::po::defvalue<std::string>( "vardist" ),
    "Prefix of output file name" )
  ;

  usr::ArgumentExtender args;
  args.AddOptions( desc );
  args.ParseOptions( argc, argv );
  args.SetFilePrefix( args.Arg( "outpre" ) );

  const std::string inputfile = args.Arg( "data" );


  std::vector<double> bias_list;
  std::vector<double> led_list;
  std::vector<double> sipm_list;

  // Reading in files
  std::string line;
  std::ifstream fin( inputfile, std::ios::in );

  while( std::getline( fin, line ) ){
    std::istringstream linestream( line );
    double time, readout, readouterr;
    double chipid, x, y, z;
    double bias, ledtemp, sipmtemp;
    linestream >> time >> chipid >> x >> y >> z
    >> bias >> ledtemp >> sipmtemp
    >> readout >> readouterr
    ;

    bias_list.push_back( bias );
    led_list.push_back( ledtemp );
    sipm_list.push_back( sipmtemp );
  }

  const double b_min = *std::min_element( bias_list.begin(), bias_list.end() );
  const double b_max = *std::max_element( bias_list.begin(), bias_list.end() );
  const double s_min = *std::min_element( sipm_list.begin(), sipm_list.end() );
  const double s_max = *std::max_element( sipm_list.begin(), sipm_list.end() );
  const double l_min = *std::min_element( led_list.begin(), led_list.end() );
  const double l_max = *std::max_element( led_list.begin(), led_list.end() );

  TH1D b_hist( usr::RandomString( 6 ).c_str(), "", 100, b_min*0.9, b_max*1.1 );
  TH1D s_hist( usr::RandomString( 6 ).c_str(), "", 100, s_min*0.9, s_max*1.1 );
  TH1D l_hist( usr::RandomString( 6 ).c_str(), "", 100, l_min*0.9, l_max*1.1 );

  for( unsigned i = 0; i < bias_list.size(); ++i ){
    b_hist.Fill( bias_list.at( i ) );
    s_hist.Fill( sipm_list.at( i ) );
    l_hist.Fill( led_list.at( i ) );
  }


  {
    usr::plt::Simple1DCanvas c;
    c.PlotHist( b_hist,
      usr::plt::PlotType( usr::plt::hist ),
      usr::plt::LineColor( usr::plt::col::blue ) );

    TSpectrum s( 20 );
    unsigned npeak = s.Search( &b_hist, 2, "nobackground" );

    for( unsigned i = 0; i < npeak; ++i ){
      c.DrawVLine( s.GetPositionX()[i],
        usr::plt::LineColor( usr::plt::col::red ) );
    }

    c.Xaxis().SetTitle( "LED Bias [mV]" );
    c.Yaxis().SetTitle( "Data points" );
    c.SaveAsPDF( args.MakePDFFile( "Bias" ) );
  }
  {
    usr::plt::Simple1DCanvas c;
    c.PlotHist( s_hist,
      usr::plt::PlotType( usr::plt::hist ),
      usr::plt::LineColor( usr::plt::col::blue ) );
    TSpectrum s( 20 );
    unsigned npeak = s.Search( &s_hist, 2, "nobackground" );

    for( unsigned i = 0; i < npeak; ++i ){
      c.DrawVLine( s.GetPositionX()[i],
        usr::plt::LineColor( usr::plt::col::red ) );
    }

    c.Xaxis().SetTitle( "SiPM Temperature [^{#circ}C]" );
    c.Yaxis().SetTitle( "Data points" );
    c.SaveAsPDF( args.MakePDFFile( "SiPM" ) );
  }
  {
    usr::plt::Simple1DCanvas c;
    c.PlotHist( l_hist,
      usr::plt::PlotType( usr::plt::hist ),
      usr::plt::LineColor( usr::plt::col::blue ) );
    TSpectrum s( 20 );
    unsigned npeak = s.Search( &l_hist, 2, "nobackground" );

    for( unsigned i = 0; i < npeak; ++i ){
      c.DrawVLine( s.GetPositionX()[i],
        usr::plt::LineColor( usr::plt::col::red ) );
    }

    c.Xaxis().SetTitle( "LED Temperature [^{#circ}C]" );
    c.Yaxis().SetTitle( "Data points" );
    c.SaveAsPDF( args.MakePDFFile( "LED" ) );
  }



  return 0;
}
