#ifndef SIPMCALIB_COMMON_STDFORMAT
#define SIPMCALIB_COMMON_STDFORMAT

#include <string>
#include <vector>

/**
 * @brief Standard File format data extraction class
 * @ingroup Common
 * @details
 *
 * The construction fo this class
 *
 */
class StdFormat
{
public:
  StdFormat( const std::string& );

  // Format for a single row.
  struct RowFormat
  {
public:
    double              time;
    int                 id;
    double              x;
    double              y;
    double              z;
    double              bias;
    double              ledtemp;
    double              sipmtemp;
    std::vector<double> data;
  };

private:
  // Standard format is basically a long list of RowFormats
  std::vector<RowFormat> _rows;

public:
  // Column extraction will be done directly using the method functions.
  // Here we also provide an interface for row selection.
  typedef bool (* RowSelect)( const RowFormat& );

  std::vector<double> Time( RowSelect                  = NoSelect ) const;
  std::vector<int>    DetId( RowSelect                 = NoSelect ) const;
  std::vector<double> X( RowSelect                     = NoSelect ) const;
  std::vector<double> Y( RowSelect                     = NoSelect ) const;
  std::vector<double> Z( RowSelect                     = NoSelect ) const;
  std::vector<double> Bias( RowSelect                  = NoSelect ) const;
  std::vector<double> LedTemp( RowSelect               = NoSelect ) const;
  std::vector<double> SiPMTemp( RowSelect              = NoSelect ) const;
  std::vector<double> DataAll( RowSelect               = NoSelect ) const;
  std::vector<double> DataCol( unsigned col, RowSelect = NoSelect ) const;


private:
  static bool NoSelect( const RowFormat& ){
    return true;
  }
};

#endif
