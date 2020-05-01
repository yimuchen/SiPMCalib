#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/Common/interface/Maths.hpp"
#include "UserUtils/PlotUtils/interface/CommonXCanvas.hpp"
#include "UserUtils/PlotUtils/interface/Flat2DCanvas.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

#include <sstream>
#include <string>
#include <vector>

#include "TGraphErrors.h"
#include "TH2D.h"
#include "TProfile.h"

int main( int argc, char* argv[] )
{
  const std::string inputfile  = argv[1];
  const std::string outputfile = argv[2];

  std::string line;
  std::ifstream fin( inputfile, std::ios::in );

  std::vector<double> zeros;
  std::vector<double> time;
  std::vector<double> readout;
  std::vector<double> read_err;
  std::vector<double> bias;
  std::vector<double> led_temp;
  std::vector<double> sipm_temp;

  double time_in;
  double readout_in;
  double read_err_in;
  double bias_in;
  double led_temp_in;
  double sipm_temp_in;

  while( std::getline( fin, line ) ){
    std::istringstream linestream( line );
    linestream >> time_in >> readout_in >> read_err_in
    >> bias_in >> led_temp_in >> sipm_temp_in;

    zeros.push_back( 0 );
    time.push_back( time_in/60 );
    readout.push_back( -readout_in );
    read_err.push_back( read_err_in );
    bias.push_back( bias_in );
    led_temp.push_back( led_temp_in );
    sipm_temp.push_back( sipm_temp_in );
  }

  TGraphErrors g_readout( time.size()
                        , time.data(), readout.data()
                        , zeros.data(), read_err.data() );
  TGraph g_bias( time.size(), time.data(), bias.data() );
  TGraph g_led( time.size(), time.data(), led_temp.data() );
  TGraph g_sipm( time.size(), time.data(), sipm_temp.data() );

  usr::plt::CommonXCanvas c( 3 );

  c.PlotGraph<0>( g_readout
                , usr::plt::TrackY( usr::plt::TrackY::both )
                , usr::plt::PlotType( usr::plt::fittedfunc )
                , usr::plt::LineColor( usr::plt::col::darkblue )
                , usr::plt::FillColor( usr::plt::col::cyan )
                , usr::plt::EntryText( "SiPM readout" ) );
  c.PlotGraph<1>( g_bias
                , usr::plt::TrackY( usr::plt::TrackY::both )
                , usr::plt::PlotType( usr::plt::simplefunc )
                , usr::plt::LineColor( usr::plt::col::darkgreen )
                , usr::plt::EntryText( "LED Bias" ) );
  c.PlotGraph<2>( g_led
                , usr::plt::TrackY( usr::plt::TrackY::both )
                , usr::plt::PlotType( usr::plt::simplefunc )
                , usr::plt::LineColor( usr::plt::col::red )
                , usr::plt::EntryText( "LED" ) );
  c.PlotGraph<2>( g_sipm
                , usr::plt::TrackY( usr::plt::TrackY::both )
                , usr::plt::PlotType( usr::plt::simplefunc )
                , usr::plt::LineColor( usr::plt::col::black )
                , usr::plt::EntryText( "SiPM" ) );

  c.Pad<0>().Yaxis().SetTitle( "[mV-ns]" );
  c.Pad<1>().Yaxis().SetTitle( "[mV]" );
  c.Pad<2>().Yaxis().SetTitle( "Temp. [^{#circ}C]" );
  c.Pad<2>().Xaxis().SetTitle( "Time [min]" );

  c.Pad<0>().SetDataMax( usr::Mean( readout ) + 2.5*usr::StdDev( readout ) );
  c.Pad<0>().SetDataMin( usr::Mean( readout ) - 2.5*usr::StdDev( readout ) );

  c.Pad<1>().SetYaxisMax( usr::Mean( bias ) + 2*usr::StdDev( bias ) );
  c.Pad<1>().SetYaxisMin( usr::Mean( bias ) - 2*usr::StdDev( bias ) );

  c.Pad<2>().SetYaxisMax( usr::Mean( led_temp ) + 3 * usr::StdDev( led_temp ) );
  c.Pad<2>().SetYaxisMin( usr::Mean( led_temp ) - 3 * usr::StdDev( led_temp ) );


  c.SaveAsPDF( outputfile );

  return 0;
}
