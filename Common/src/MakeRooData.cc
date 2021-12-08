#include "SiPMCalib/Common/interface/MakeRooData.hpp"
#include "UserUtils/Common/interface/Maths.hpp"
#include "UserUtils/Common/interface/STLUtils/VectorUtils.hpp"

/**
 * @ingroup Common
 * @brief Setting a RooRealVar range based on a series of data points.
 *
 * Given a bin width, a maximum area, and a series of data values, this
 * functions
 * automatically setups up the RooRealVar var such that:
 * - The range of variable lands on a bin edge that are integer multiples of the
 *   specified bin with.
 * - The range of variable minimally covers the maximum and minimum of the data
 *   points (maximum value can be over written )
 */
extern void
SetRange( RooRealVar&                var,
          const double               binwidth,
          const double               maxarea,
          const std::vector<double>& datapoints )
{
  // Getting the maximum elements.
  const double dmax = usr::GetMaximum( datapoints );
  const double dmin = usr::GetMinimum( datapoints );
  const double m    = maxarea < 0 ?
                      dmax :
                      maxarea;

  // Additional parsing to be done of the
  const double xmin  = usr::RoundDown( dmin, binwidth );
  const double xmax  = usr::RoundUp( std::min( dmax, m ), binwidth );
  const double nbins = ( xmax-xmin ) / binwidth;

  // Setting the range
  var.setRange( xmin, xmax );
  var.setBins( nbins );
}


/**
 * @ingroup Common
 * @brief Given a formally setup of RooRealVar and a vector of data, create a
 * RooDataHist object.
 */
extern RooDataHist*
MakeData( RooRealVar&                var,
          const std::vector<double>& data,
          const double               maxarea )
{
  RooDataHist*ans = new RooDataHist( "data", "data", RooArgList( var ) );

  for( const auto x : data ){
    if( maxarea < 0 || x < maxarea ){
      var = x;
      ans->add( RooArgList( var ) );
    }
  }

  return ans;
}
