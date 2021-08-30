#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/Common/interface/Maths.hpp"
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

#include <fstream>
#include <sstream>

int
main( int argc, char* argv[] )
{
  usr::po::options_description desc(
    "Plotting instructions for ENF display" );
  desc.add_options()
    ( "input,i", usr::po::value<std::string>(),
    "Input file containing the values of ENF calculation values" )
    ( "alpha,a", usr::po::multivalue<double>(),
    "alpha value to place in plot" )
    ( "lambda,l", usr::po::multivalue<double>(),
    "lambda values to place in plot" )
    ( "sigma,s", usr::po::multivalue<double>(),
    "sigma values to place in plot" )
    ( "output,o", usr::po::defvalue<std::string>( "ENFplot.pdf" ),
    "Output file name" )
  ;

  usr::ArgumentExtender args;
  args.AddOptions( desc );
  args.ParseOptions( argc, argv );

  const std::vector<double> alpha_list  = args.ArgList<double>( "alpha" );
  const std::vector<double> lambda_list = args.ArgList<double>( "lambda" );
  const std::vector<double> sigma_list  = args.ArgList<double>( "sigma" );

  // Allocating the graphs
  std::map<uint64_t, TGraph> graphs;
  std::map<uint64_t, unsigned> npoints;

  for( const auto a : alpha_list ){
    for( const auto l : lambda_list ){
      for( const auto s : sigma_list ){
        uint64_t hash = usr::OrderedHash64( {a, l, s} );
        graphs[hash]  = TGraph( 100 );
        npoints[hash] = 0;
      }
    }
  }

  // Reading the lines
  std::ifstream infile( args.Arg<std::string>( "input" ) );
  std::string line;

  while( std::getline( infile, line ) ){
    std::istringstream iss( line );
    double mean, alpha, lambda, sigma, enf;
    iss >> mean >> lambda >> alpha >> sigma >> enf;
    const uint64_t hash = usr::OrderedHash64( {alpha, lambda, sigma} );

    if( graphs.count( hash ) ){
      std::cout << "Adding point to graph at " << hash << std::endl;
      graphs[hash].SetPoint( npoints[hash], mean, enf );
      npoints[hash] = npoints[hash] + 1;
    } else {
    }
  }

  // Plotting all the stuff.
  usr::plt::Simple1DCanvas c;
  std::vector<int> color_list = {
    usr::plt::col::blue,
    usr::plt::col::red,
    usr::plt::col::green,
    usr::plt::col::purple,
    usr::plt::col::orange,
    usr::plt::col::black,
  };

  unsigned plot_index = 0;

  for( const auto a : alpha_list ){
    for( const auto l : lambda_list ){
      for( const auto s : sigma_list ){
        uint64_t hash = usr::OrderedHash64( {a, l, s} );
        auto& graph   = graphs[hash];

        graph.Set( npoints[hash] );// Truncating the number of points
        graph.Sort();

        // Making the title stuff.
        std::string title = "";

        if( alpha_list.size() > 1 ){
          title += usr::fstr( "P_{a.p.} =%.1lf%%", a * 100 );
        }
        if( lambda_list.size() > 1 ){
          title += usr::fstr( "P_{c.t.} = %.1lf%%", l * 100 );
        }
        if( sigma_list.size() > 1 ){
          title += usr::fstr( "#sigma_{0} = %.2lf gain", s );
        }

        // The full graph
        c.PlotGraph( graph,
          usr::plt::PlotType( usr::plt::simplefunc ),
          usr::plt::LineColor( color_list[plot_index] ),
          usr::plt::TrackY( usr::plt::both ),
          usr::plt::EntryText( title ) );

        // Draw the ideal line,
        const double ideal_enf = 1./( 1-l );
        c.DrawHLine( ideal_enf,
          usr::plt::LineColor( color_list[plot_index] ),
          usr::plt::LineStyle( usr::plt::sty::lindashed ) );

        plot_index++;
      }
    }
  }

  // Additional adjustments to the plots
  c.Pad().Xaxis().SetTitle( "Mean number of photons" );
  c.Pad().Yaxis().SetTitle( "Excess Noise Factor" );
  c.Pad().DrawCMSLabel( "Simulation", "HGCAL" );
  c.Pad().SetYaxisMin( 0.9 );

  // Adding additional legends
  if( alpha_list.size() == 1 ){
    c.Pad().WriteLine( usr::fstr( "P_{a.p.} = %.1lf%%", alpha_list[0] * 100 ) );
  }
  if( lambda_list.size() == 1 ){
    c.Pad().WriteLine( usr::fstr( "P_{c.t.} = %.1lf%%", lambda_list[0] * 100 ) );
  }
  if( sigma_list.size() == 1 ){
    c.Pad().WriteLine( usr::fstr( "#sigma_{0} = %.2lf gain", sigma_list[0] ) );
  }

  c.SaveAsPDF( args.Arg<std::string>( "output" ) );

  return 0;
}
