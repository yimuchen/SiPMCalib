#include "SiPMCalib/SiPMCalc/interface/SiPMPdf.hpp"
#include "UserUtils/Common/interface/STLUtils/OStreamUtils.hpp"

double
mean( double mu, double lambda )
{
  double sum = 0;
  int max    = mu + 1000*std::sqrt( mu );

  for( int i = 0; i < max; ++i ){
    sum += i * SiPMPdf::GeneralPoissonProb( i, mu, lambda );
  }

  return sum;
}

double
variance( double mu, double lambda )
{
  const double std = std::sqrt( mu );
  int range_min    = std::max( mu - 5*std, 0.00 );
  int range_max    = std::max( mu + 100*std, 30.0 );

  double sum  = 0;
  double sum2 = 0;

  for( int i = range_min; i <= range_max; ++i  ){
    const double prob = SiPMPdf::GeneralPoissonProb( i, mu, lambda );
    sum  += i * prob;
    sum2 += i*i*prob;
  }

  return sum2 - sum * sum;
}

double
calcvar( double mean, double lambda )
{
  return mean / std::pow( 1-lambda, 3 );
}


double
enf( double mu, double lambda )
{
  return variance( mu, lambda )/mean( mu, lambda );
}

int
main()
{
  for( double mu = 10; mu < 15; mu += 1 ){
    usr::fout( "%5.2f | ", mu );

    for( double lambda = 0; lambda < 0.15; lambda += 0.025 ){
      usr::fout( "\t%8.05f(%8.05f)", variance( mu, lambda )/mu
               , 1.0/std::pow( 1.0-lambda, 3 )
        );

    }

    usr::fout( "\n" );
  }

  std::cout << "Descrete information" << std::endl;

  for( double mu = 10; mu < 20; mu += 1 ){
    usr::fout( "%5.2f | ", mu );

    for( double lambda = 0; lambda < 0.15; lambda += 0.025 ){
      usr::fout( "\t%8.05f", enf( mu, lambda ) );
    }

    usr::fout( "\n" );
  }

  return 0;
}
