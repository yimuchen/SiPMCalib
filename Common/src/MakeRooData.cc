#include "SiPMCalib/Common/interface/MakeRooData.hpp"
#include "UserUtils/Common/interface/Maths.hpp"
#include "UserUtils/Common/interface/STLUtils/VectorUtils.hpp"

extern void
SetRange( RooRealVar&                var,
          const double               binwidth,
          const double               maxarea,
          const std::vector<double>& datapoints )
{
  // Getting the maximum elements.
  const double dmax = usr::GetMaximum( datapoints );
  const double dmin = usr::GetMinimum( datapoints );
  const double m    = maxarea < 0 ? dmax : maxarea;

  // Additional parsing to be done of the
  const double xmin  = usr::RoundDown( dmin, binwidth );
  const double xmax  = usr::RoundUp( std::min( dmax, m ), binwidth );
  const double nbins = ( xmax - xmin ) / binwidth;

  // Setting the range
  var.setRange( xmin, xmax );
  var.setBins( nbins );

}

extern RooDataHist*
MakeData( RooRealVar&                var,
          const std::vector<double>& data,
          const double               maxarea )
{
  RooDataHist* ans = new RooDataHist( "data", "data", RooArgList( var ) );

  for( const auto x : data ){
    if( maxarea < 0 || x < maxarea ){
      var = x;
      ans->add( RooArgList( var ) );
    }
  }

  return ans;
}
