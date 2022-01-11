#ifndef SIPMCALIB_SIPMCALC_SIPMNONLINEARFIT_HPP
#define SIPMCALIB_SIPMCALC_SIPMNONLINEARFIT_HPP
#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/MathUtils/interface/Measurement/Measurement.hpp"

#include "SiPMCalib/Common/interface/StdFormat.hpp"

#include "Math/Interpolator.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TProfile3D.h"

#include <string>

class SiPMNonLinearFit
{
public:
  SiPMNonLinearFit();
  SiPMNonLinearFit( const std::string& file );

  // Update settings
  static usr::po::options_description FitArguments();
  static usr::po::options_description OutputArguments();

  void UpdateSettings( const usr::ArgumentExtender& );

  /**
   * @{
   * @brief Computation control flow
   */
  void ReadFiles( const std::string& zscan, const std::string& corr = "" );
  void MakeLinearGraph();
  void MakeNonLinearGraph( const double ref_n,
                           const double ref_z,
                           const double ref_bias,
                           const double gain,
                           const double ped );
  void RunNonLinearFit();

  /** @} */

  /**
   * @{
   * plotting variations
   */
  void PlotOriginal( const std::string& );
  void PlotLinearity( const std::string& );
  void PlotNonLinearity( const std::string& );

  /** @} */

private:
  std::unique_ptr<StdFormat> _raw_data;    //
  TGraphErrors               _nl_data; // Container for the zscan data
  TGraphErrors               _lin_data;
  TF1                        _nl_func;
  TF1                        _lin_func;
  TFitResult                 _nl_fit;    // Dynamically reassign
  TFitResult                 _lin_fit;    // Dynamically reassign

  // Settings for the linear-style fit
  double _lin_pmin;
  double _lin_pmax;
  double _power_z;
  double _pixel_min;
  double _pixel_max;

  // Settings for plotting
  std::string _readout_units;
  std::string _setup;
  std::string _sipm_model;
  std::string _sipm_id;

  // Classes for performing data collection
  std::unique_ptr<ROOT::Math::Interpolator> _lookup_1d;
  std::unique_ptr<TProfile3D>               _lookup_3d;


  // Helper functions
  void InitDefaults();
  void MakeInterpolator1D();
  void MakeInterpolator3D( const std::string& corr );

  double GetZPosMultiplier( const StdFormat::RowFormat& main,
                            const StdFormat::RowFormat& ref,
                            const double                ped ) const;
  double GetLookUpMultiplier( const StdFormat::RowFormat& main,
                              const StdFormat::RowFormat& ref ) const;
};


#endif
