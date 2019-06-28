#include "UserUtils/PlotUtils/interface/Flat2DCanvas.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

#include <boost/format.hpp>
#include <fstream>
#include <vector>

#include "RooDataSet.h"
#include "RooRealVar.h"


int
main( int argc, char* argv[] )
{
  std::ifstream fin;
  std::string line;
  fin.open( argv[1], std::ifstream::in );

  int linecount = 0;
  unsigned t, ncapture, presample, postsample;
  int8_t v;
  double convfactor;
  std::vector<std::vector<int16_t> > vecprofile;

  // Getting first line
  std::getline( fin, line );
  std::istringstream linestream( line );
  linestream >> t >> ncapture >> presample >> postsample >> convfactor;

  // Getting all other lines
  while( std::getline( fin, line ) ){
    ++linecount;
    if( linecount % 1000 == 0 ){
      std::cout << "\rLine " <<  linecount << "..." << std::flush;
    }
    vecprofile.push_back( std::vector<int16_t>( presample + postsample,  0 ) );

    for( unsigned i = 0; i < presample + postsample; ++i ){
      v = 0;
      const int16_t v1 = line[2*i] >= 'a' ? 10 + line[2*i] - 'a' :
                         line[2*i] >= 'A' ? 10 + line[2*i] - 'A' :
                         line[2*i] - '0';
      const int16_t v2 = line[2*i+1] >= 'a' ? 10 + line[2*i+1] - 'a' :
                         line[2*i+1] >= 'A' ? 10 + line[2*i+1] - 'A' :
                         line[2*i+1] - '0';
      v                    = v1 << 4 | v2;
      vecprofile.back()[i] = v;
    }
  }

  std::cout << "Done reading!" << std::endl;

  linecount = 0;

  RooRealVar area( "area", "Pulse area", -400, 400, "ADC-ns" );
  RooDataSet data( "data", "data", RooArgSet( area ) );
  area.setBins( 50 );

  // Filling in the profile
  for( const auto& profile : vecprofile ){
    ++linecount;

    if( linecount % 1000 == 0 ){
      std::cout << "\rFilling profile [" << linecount << "]..." << std::flush;
    }

    double profilearea = 0;

    for( unsigned i = 0; i < profile.size(); ++i ){
      profilearea += ( -(double)profile.at( i ) ) * t ;
    }

    area = profilearea;
    data.add( RooArgSet( area ) );

  }

  {
    usr::plt::Simple1DCanvas c( usr::plt::RangeByVar( area ),
                               usr::plt::Simple1DCanvas::default_width  );
    auto& dgraph = c.PlotData( data );

    dgraph.SetMarkerStyle( usr::plt::sty::mkrcircle );
    dgraph.SetMarkerSize( 0.2 );

    c.DrawLuminosity( "Laser setup" );
    c.DrawCMSLabel( "Preliminary", "HGCal" );
    c.Pad().WriteLine( "Bias voltage: 44V" );
    c.Pad().WriteLine( "(Nominal: 54.04V)" );

    c.SaveAsPDF( "noiseprofile.pdf" );
  }

  return 0 ;
}
