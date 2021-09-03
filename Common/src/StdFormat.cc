#include "SiPMCalib/Common/interface/StdFormat.hpp"

#include <fstream>
#include <sstream>

StdFormat::StdFormat( const std::string& filename )
{
  // Getting everything row by row
  std::ifstream infile( filename );
  std::string line;

  while( std::getline( infile, line ) ){
    std::istringstream linestream( line );
    RowFormat r;
    double data = 0;
    linestream >> r.time >> r.id
    >> r.x >> r.y >> r.z
    >> r.bias >> r.ledtemp  >> r.sipmtemp;

    while( linestream >> data ){
      r.data.push_back( data );
    }

    _rows.push_back( r );
  }
}

// Macro for generating the column selector
#define COLUMN( FNAME, TYPE, MEMBER )          \
  std::vector<TYPE>                            \
  StdFormat::FNAME( RowSelect selector ) const \
  {                                            \
    std::vector<TYPE> ans;                     \
    ans.reserve( _rows.size() );               \
    for( const auto& row : _rows ){            \
      if( selector( row ) ){                   \
        ans.push_back( row.MEMBER );           \
      }                                        \
    }                                          \
    return ans;                                \
  }

COLUMN( Time,     double, time     );
COLUMN( DetId,    int,    id       );
COLUMN( X,        double, x        );
COLUMN( Y,        double, y        );
COLUMN( Z,        double, z        );
COLUMN( Bias,     double, bias     );
COLUMN( LedTemp,  double, ledtemp  );
COLUMN( SiPMTemp, double, sipmtemp );

// Special case for data manipulations
std::vector<double>
StdFormat::DataAll( RowSelect selector ) const
{
  std::vector<double> ans;

  for( const auto& row : _rows ){
    if( selector( row ) ){
      ans.insert( ans.end(), row.data.begin(), row.data.end() );
    }
  }

  return ans;
}

// Special case for data manipulations
std::vector<double>
StdFormat::DataCol( unsigned col, RowSelect selector ) const
{
  std::vector<double> ans;
  ans.reserve( _rows.size() );

  for( const auto& row : _rows ){
    if( selector( row ) ){
      ans.push_back( row.data[col] );
    }
  }

  return ans;
}
