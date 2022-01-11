#ifndef SIPMCALIB_SIPMCALC_SIPMLOWLIGHTFIT_HPP
#define SIPMCALIB_SIPMCALC_SIPMLOWLIGHTFIT_HPP

#include "SiPMCalib/SiPMCalc/interface/SiPMPdf.hpp"
#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/MathUtils/interface/Measurement/Measurement.hpp"

#include "RooAbsData.h"
#include "RooRealVar.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TSpectrum.h"

#include <iostream>
#include <memory>

/**
 * @brief Class for handling LowLight fit requestion, including batch fitting
 * for data with similar configurations, and estimation of formats.
 */
class SiPMLowLightFit
{
public:
  SiPMLowLightFit();
  SiPMLowLightFit( const std::string& config );

  // Update settings
  static usr::po::options_description DataArguments();
  static usr::po::options_description FitArguments();
  static usr::po::options_description OperationArguments();
  static usr::po::options_description OutputArguments();
  static usr::po::options_description EstArguments();

  void UpdateSettings( const usr::ArgumentExtender& );

  /**
   * @{
   * @brief main control flow functions.
   */
  void MakeBinnedData();
  void RunPDFEstimation();
  void RunFit();
  /** @} */

  /**
   * @{
   * @brief Plotting results for the SiPM analysis
   */
  void PlotPeakFind( const std::string& );
  void PlotGainFit( const std::string& );
  void PlotWidthFit( const std::string& );
  void PlotPoissonFit( const std::string& );
  void PlotSpectrumFit( const std::string& );

  /** @} */

  /**
   * @{
   * @brief Stream output for fit results
   */
  void PrintRaw( std::ostream& sout   = std::cout ) const;
  void PrintTable( std::ostream& sout = std::cout ) const;
  void PrintFitResults();

  /** @} */

  /**
   * @{
   * @brief  Access to fitting variables
   */
  RooRealVar& x()     { return *_x;      }
  RooRealVar& ped()   { return *_ped;    }
  RooRealVar& gain()  { return *_gain;   }
  RooRealVar& s0()    { return *_s0;     }
  RooRealVar& s1()    { return *_s1;     }
  RooRealVar& mean()  { return *_mean;   }
  RooRealVar& lambda(){ return *_lambda; }
  RooRealVar& alpha() { return *_alpha;  }
  RooRealVar& beta()  { return *_beta;   }
  RooRealVar& dcfrac(){ return *_dcfrac; }
  RooRealVar& eps()   { return *_eps;    }

  const RooRealVar&
  x()      const { return *_x;      }
  const RooRealVar&
  ped()    const { return *_ped;    }
  const RooRealVar&
  gain()   const { return *_gain;   }
  const RooRealVar&
  s0()     const { return *_s0;     }
  const RooRealVar&
  s1()     const { return *_s1;     }
  const RooRealVar&
  mean()   const { return *_mean;   }
  const RooRealVar&
  lambda() const { return *_lambda; }
  const RooRealVar&
  alpha()  const { return *_alpha;  }
  const RooRealVar&
  beta()   const { return *_beta;   }
  const RooRealVar&
  dcfrac() const { return *_dcfrac; }
  const RooRealVar&
  eps()    const { return *_eps;    }

  /** @} */

  /**
   * @{
   * @brief Returing physically interpreted fit paramters in measurement
   * containers
   **/
  usr::Measurement Pedestal() const;
  usr::Measurement Gain() const;
  usr::Measurement CommonNoise() const;
  usr::Measurement PixelNoise() const;
  usr::Measurement MeanPhotons() const;
  usr::Measurement MeanPhotonsFromFrac() const;
  usr::Measurement ProbCrosstalk() const;
  usr::Measurement ProbAfterpulse() const;
  usr::Measurement AfterpulseTimeNS() const;
  usr::Measurement DarkcurrentTimeNS() const;
  usr::Measurement ExcessNoiseFactor( const usr::Measurement& ) const;

  /** @} */

private:
  // Since RooFit object declaration after additional parsing to get the
  // requested range, RooFit objects must use pointer interfaces. Using
  // unique_ptr to handle memory management
  std::unique_ptr<RooRealVar> _ped;
  std::unique_ptr<RooRealVar> _gain;
  std::unique_ptr<RooRealVar> _s0;
  std::unique_ptr<RooRealVar> _s1;
  std::unique_ptr<RooRealVar> _mean;
  std::unique_ptr<RooRealVar> _lambda;
  std::unique_ptr<RooRealVar> _alpha;
  std::unique_ptr<RooRealVar> _beta;
  std::unique_ptr<RooRealVar> _dcfrac;
  std::unique_ptr<RooRealVar> _eps;
  std::unique_ptr<SiPMPdf>    _pdf;

  // The data format to be used by the used by the fit
  std::vector<double>         _arealist;
  std::unique_ptr<RooAbsData> _data;
  std::unique_ptr<RooRealVar> _x;


  // Default setting options
  void set_all_defaults();

  // data parsing settings
  void make_array_from_waveform();
  void make_array_from_sum();

  /**
   * @{
   * @brief  Objects used for fit function estimations.
   */
  std::vector<std::unique_ptr<TF1> > _peakfits;
  std::unique_ptr<TSpectrum>         _spectrum;
  std::unique_ptr<TH1D>              _est_hist;
  std::unique_ptr<TGraphErrors>      _gain_graph;
  std::unique_ptr<TGraphErrors>      _width_graph;
  std::unique_ptr<TGraphErrors>      _height_graph;
  std::unique_ptr<TF1>               _gain_fit;
  std::unique_ptr<TF1>               _width_fit;
  std::unique_ptr<TF1>               _height_fit;
  bool                               ignore_ped_est;
  bool                               ignore_gain_est;
  bool                               ignore_s0_est;
  bool                               ignore_s1_est;
  bool                               ignore_mean_est;
  bool                               ignore_lambda_est;
  double                             _est_minpeak;
  int                                _est_gausswindow;
  int                                _est_maxgausswidth;

  /** @} */

  /**
   * @{
   * @brief helper functions for estimation routines
   */
  TF1* good_local_peak_fit( double x );
  void run_gain_est();
  void run_width_est();
  void run_height_est();
  /** @} */

  // Options for reading data formats.
  std::string _inputfile;
  bool        _waveform;
  double      _binwidth;
  unsigned    _timeint;
  unsigned    _intstart;
  unsigned    _intstop;
  unsigned    _pedstart;
  unsigned    _pedstop;
  double      _pedrms;
  double      _maxarea;

  // operation parameters
  double      _intwindow;
  double      _sipmtime;
  std::string _sipmtype;
  std::string _lumitype;
  std::string _biasv;
};

#endif
