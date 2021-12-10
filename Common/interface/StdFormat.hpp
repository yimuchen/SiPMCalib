#ifndef SIPMCALIB_COMMON_STDFORMAT
#define SIPMCALIB_COMMON_STDFORMAT

#include <functional>
#include <string>
#include <vector>

/**
 * @brief Standard File format data extraction class
 * @ingroup Common
 * @details
 *
 * The standard file format expects data to be organized into columns with the
 * following format:
 *
 * - time (relative to data collection command initialization) [seconds]
 * - detector element ID.
 * - The gantry x, y, z coordinates for a single data collection. [mm]
 * - The measured bias voltage of the LED system [mV]
 * - The measured LED temperature [C]
 * - The measured SiPM temperature [C]
 * - The remaining N columns are generically called "Data" and will depend on
 *   what data collection routine is used for generate the data file.
 *
 * Each row in the data file will be organized into the a RowFormat object. And
 * the various function can be used to extract columns into a standard
 * std::vector container, to be passed to other computation routines. Users can
 * also pass a row selection function pointer to the column extraction functions
 * to limit the rows of data used for some data collection routine.
 */
class StdFormat
{
public:
  StdFormat( const std::string& );

  /**
   * @brief  Simple container for a single row of data in a standard format data
   * file.
   */
  struct RowFormat
  {
public:
    double time;
    int id;
    double x;
    double y;
    double z;
    double bias;
    double ledtemp;
    double sipmtemp;
    std::vector<double> data;
  };

private:
  std::vector<RowFormat> _rows;

public:
  typedef std::function<bool ( const RowFormat& )> RowSelect;

  /**< Data selection short hand, function should return true for the rows that
   * should be extracted*/


  /**
   * @{
   * @brief Extracting some column into a std::vector container.
   * @details User can provide a row selection function to extract specific rows
   * of interest.
   */
  std::vector<double> Time( RowSelect     = NoSelect ) const;
  std::vector<int>    DetId( RowSelect    = NoSelect ) const;
  std::vector<double> X( RowSelect        = NoSelect ) const;
  std::vector<double> Y( RowSelect        = NoSelect ) const;
  std::vector<double> Z( RowSelect        = NoSelect ) const;
  std::vector<double> Bias( RowSelect     = NoSelect ) const;
  std::vector<double> LedTemp( RowSelect  = NoSelect ) const;
  std::vector<double> SiPMTemp( RowSelect = NoSelect ) const;

  /** @} */

  std::vector<double> DataAll( RowSelect               = NoSelect ) const;
  std::vector<double> DataCol( unsigned col, RowSelect = NoSelect ) const;

  StdFormat MakeReduced( RowSelect ) const;
  void      WriteToFile( const std::string& filename ) const;

  /**
   * @brief Default row selection that does no explicit selection.
   */
  static bool NoSelect( const RowFormat& ){  return true; }

  /**
   * @{
   * @brief old school interface for looping over rows
   */
  inline std::vector<RowFormat>::const_iterator
  begin() const { return _rows.begin(); }
  inline std::vector<RowFormat>::const_iterator
  end() const { return _rows.end(); }

  /** @} */

private:
  StdFormat(); // Bare construction for reduced
};

#endif
