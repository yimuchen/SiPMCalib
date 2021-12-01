#include "SiPMCalib/Common/interface/StdFormat.hpp"

#include "UserUtils/Common/interface/STLUtils/OStreamUtils.hpp"
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"

#include <exception>
#include <fstream>
#include <sstream>

/**
 * @brief Construct a new StdFormat from a file path.
 *
 * Should the file be unable to open, either because of permission issues or
 * because the file doesn't exist, this will raise an exception.
 */
StdFormat::StdFormat( const std::string& filename )
{
  // Getting everything row by row
  std::ifstream infile( filename );
  std::string line;

  if( !infile.is_open() ){
    usr::log::PrintLog( usr::log::FATAL,// Exception will be thrown
      usr::fstr( "Input file %s cannot be opened!", filename ) );
  }

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

/**
 * @brief Extracting a specific column of the data columns.
 *
 * Depending on which data collection process was invoked to collect the data,
 * different columns typically mean different data variation. This method extract
 * exactly 1 data column into the a single vector container. Notice that the
 * column starts from the 0 for the 0th data column (9th column in the data file.)
 */
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

/**
 * @brief Extracting all data columns as a single vector.
 *
 * In the case of the low-light data collection, all waveform data are placed
 * into the latter columns. This method aggregates all data columns into a single
 * vector, effectively removing all column and row structure (what is needed for
 * low-light analysis.) The user can still be specify which rows are used for
 * extraction of the data vector.
 */
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
