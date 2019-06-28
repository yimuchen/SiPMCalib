#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/PlotUtils/interface/Flat2DCanvas.hpp"

#include <boost/format.hpp>
#include <fstream>
#include <vector>

int
main( int argc, char* argv[] )
{
  usr::po::options_description desc( "Options for plotting head profile" );
  desc.add_options()
    ( "data", usr::po::value<std::string>(), "Data .txt file" )
    ( "output", usr::po::value<std::string>(), "Output file name" )
    ( "start", usr::po::value<unsigned>()->default_value( 0 ), "starting samples" )
    ( "end", usr::po::value<unsigned>()->default_value( 2147483647 ), "ending value" )
    ( "amax", usr::po::value<double>()->default_value( 1e6 ), "maximum area" )
    ( "amin", usr::po::value<double>()->default_value( -1e6 ), "minimum area" )
  ;

  usr::ArgumentExtender arg;
  arg.AddOptions( desc );
  arg.ParseOptions( argc, argv );

  if( !arg.CheckArg( "data" ) ){
    std::cerr << "Please provide input data" << std::endl;
    return 1;
  }
  if( !arg.CheckArg( "output" ) ){
    std::cerr << "Please provide output file name " << std::endl;
    return 1;
  }

  std::ifstream fin;
  std::string line;
  fin.open( arg.Arg( "data" ), std::ifstream::in );

  int linecount = 0;
  unsigned t, ncapture, presample, postsample;
  double convfactor;
  int8_t max = 1 << 7;
  int8_t min = ~max;
  std::vector<std::vector<int16_t> > vecprofile;

  // Getting first line
  std::getline( fin, line );
  std::istringstream linestream( line );
  linestream >> t >> ncapture >> presample >> postsample >> convfactor;

  const unsigned start = arg.Arg<unsigned>( "start" );
  const unsigned end   = std::min( arg.Arg<unsigned>( "end" ),
    presample + postsample );
  const unsigned nsample = end - start;


  // Getting all other lines
  while( std::getline( fin, line ) ){
    ++linecount;
    if( linecount % 1000 == 0 ){
      std::cout << "\rLine " <<  linecount << "..." << std::flush;
    }
    vecprofile.push_back( std::vector<int16_t>() );
    vecprofile.back().reserve( nsample );

    for( unsigned i = start; i < end; ++i ){
      const int16_t v1 = line[2*i] >= 'a' ? 10 + line[2*i] - 'a' :
                         line[2*i] >= 'A' ? 10 + line[2*i] - 'A' :
                         line[2*i] - '0';
      const int16_t v2 = line[2*i+1] >= 'a' ? 10 + line[2*i+1] - 'a' :
                         line[2*i+1] >= 'A' ? 10 + line[2*i+1] - 'A' :
                         line[2*i+1] - '0';
      const int8_t v = v1 << 4 | v2;
      vecprofile.back().push_back( v);
      max                  = std::max( v, max );
      min                  = std::min( v, min );
    }

  }

  std::cout << "Done reading!" << std::endl;

  // Creating historgram
  const double ymax = (double)( max + 1 ) * convfactor * 256;
  const double ymin = (double)( min - 1 ) * convfactor * 256;

  const double xmin = (double)t*( (double)start - presample );
  const double xmax = (double)t*( end - presample );
  std::cout << xmin << " " << xmax << std::endl;
  std::cout << ymin << " " << ymax << std::endl;

  TH2D hist( "hist", "",
             nsample, xmin, xmax,
             ( max-min ), ymin, ymax
             );

  linecount = 0;

  const double amax = arg.Arg<double>( "amax" );
  const double amin = arg.Arg<double>( "amin" );

  // Filling in the profile
  for( const auto& profile : vecprofile ){
    ++linecount;

    if( linecount % 1000 == 0 ){
      std::cout << "\rFilling profile [" << linecount << "]..." << std::flush;
    }

    double area = 0;

    for( unsigned i = 0; i < nsample; ++i ){
      area += (double)t * ( profile.at( i ) );
    }

    if( amax > area && area > amin ){

      for( unsigned i = 0; i < nsample; ++i ){
        const double x = (double)t*( ((double)start - presample) + i );
        const double y = (double)( profile.at( i ) ) * convfactor * 256.;
        if( x > xmax || x < xmin || y < ymin || y > ymax ){
          std::cout << " Weird point found!" << x << " " << y << std::endl;
        }
        hist.Fill( x, y );
      }

    }
  }

  // Setting the zero bins for aesthetics
  for( int i = 0; i < hist.GetNcells(); ++i ){
    if( hist.GetBinContent( i ) == 0 ){
      hist.SetBinContent( i, 0.3 );  // Half a order of magnitude smaller
    }
  }

  std::cout << "Done" << std::endl;
  hist.SetContour( 50 );


  usr::plt::Flat2DCanvas c;
  c.PlotHist( hist,
    usr::plt::Plot2DF( usr::plt::heat ) );

  c.DrawCMSLabel( "Preliminary", "HGCal" );
  c.DrawLuminosity( "Laser setup" );

  c.Pad().SetTextCursor( 0.05, 0.9, usr::plt::font::top_left );
  c.Pad().WriteLine(
    ( boost::format( "ADC bit=%.3lf[mV]" )%( convfactor*256 ) ).str() )
  .WriteLine( ( boost::format( "Max: #pm %d[mV]" )
                %( convfactor*256*127 ) ).str() )
  .WriteLine( "Sample rate: 2ns^{-1}" )
  .WriteLine( "Trigger rate: 600#mu s^{-1}" );

  c.Xaxis().SetTitle( "Time [ns]" );
  c.Yaxis().SetTitle( "Readout [mV]" );
  c.Zaxis().SetTitle( "Events" );

  c.SetLogz( 1 );

  c.SaveAsPDF( arg.Arg("output") );

}
