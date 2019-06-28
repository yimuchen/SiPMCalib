#include "UserUtils/Common/interface/Maths.hpp"
#include "UserUtils/MathUtils/interface/RooFitExt.hpp"
#include "UserUtils/PlotUtils/interface/Ratio1DCanvas.hpp"

#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooRealVar.h"

#include <boost/format.hpp>
#include <vector>

int
main( int argc, char** argv )
{
  std::vector<double> regsamp;
  std::vector<double> randsamp;
  double val;

  for( unsigned i = 0; i < 1000; ++i ){
    std::cin >> val;
    regsamp.push_back( val * 1000 );
    std::cin >> val;
    randsamp.push_back( val * 1000 );
    // if( val < 0 ){ std::cout << val << std::endl;}
  }

  const double xmax = std::max(
    *std::max_element( regsamp.begin(),  regsamp.end() ),
    *std::max_element( randsamp.begin(), randsamp.end() )
    );
  const double xmin = std::min(
    *std::min_element( regsamp.begin(),  regsamp.end() ),
    *std::min_element( randsamp.begin(), randsamp.end() )
    );

  std::cout << xmax << "  " << xmin << std::endl;

  RooRealVar x( "x", "Readout", xmin, xmax, "pA" );
  RooDataSet reg( "reg", "reg", RooArgSet( x ) );
  RooDataSet rnd( "rnd", "rnd", RooArgSet( x ) );

  RooRealVar regmean( "regmean", "regmean",
                      usr::Mean( regsamp ), xmin, xmax );
  RooRealVar regsig( "regsig", "regsig",
                     usr::StdDev( regsamp ), 0.000001, xmax );
  RooGaussian reggauss( "regg", "regg", x, regmean, regsig );

  RooRealVar rndmean( "rndmean", "rndmean",
                      usr::Mean( randsamp ), xmin, xmax );
  RooRealVar rndsig( "rndsig", "rndsig",
                     usr::StdDev( randsamp ), 0.000001, xmax );
  RooGaussian rndgauss( "rndg", "rndg", x, regmean, regsig );

  for( unsigned i = 0; i < 1000; ++i ){
    x = regsamp.at( i );
    reg.add( RooArgSet( x ) );
    x = randsamp.at( i );
    rnd.add( RooArgSet( x ) );
  }

  reggauss.fitTo( reg );
  rndgauss.fitTo( rnd );


  reggauss.fitTo( reg );
  rndgauss.fitTo( rnd );

  std::cout << "KS-REG" << " " << usr::KSDistance( reg, reggauss, x ) << std::endl;
  std::cout << "KS-RND" << " " << usr::KSDistance( rnd, rndgauss, x ) << std::endl;
  return 0;
}
