#include "SiPMCalib/Common/interface/StdFormat.hpp"
#include "SiPMCalib/Common/interface/WaveFormat.hpp"
#include "SiPMCalib/SiPMCalc/interface/SiPMLowLightFit.hpp"

#include "UserUtils/Common/interface/Maths.hpp"
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"
#include "UserUtils/MathUtils/interface/RooFitExt.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>

#include "RooDataHist.h"

/**
 * @class SiPMLowLightFit
 * @ingroup SiPMCalc
 * @brief Fitting routine class for fitting spectral data to the SiPM low light
 * PDF.
 *
 * @details This is an example.
 *
 */

void
SiPMLowLightFit::set_all_defaults()
{
  // Default options data parsing
  _inputfile = "";
  _waveform  = true;
  _binwidth  = 16;
  _intstart  = 0;
  _intstop   = -1;
  _pedstart  = -1;
  _pedstop   = -1;
  _pedrms    = 0.5;
  _maxarea   = 2147483647;

  // Fitting related options
  const double min = std::numeric_limits<double>::min();

  _x = std::make_unique<RooRealVar>( "Readout",
                                     "Readout",
                                     -1000,
                                     200000,
                                     "mV-ns" );
  _ped =
    std::make_unique<RooRealVar>(     "ped",     "ped", -300,      300 );
  _gain =
    std::make_unique<RooRealVar>(    "gain",    "gain", 0,        1000 );
  _s0 =
    std::make_unique<RooRealVar>(      "s0",      "s0", min,        15 );
  _s1 =
    std::make_unique<RooRealVar>(      "s1",      "s1", 0,           5 );
  _mean =
    std::make_unique<RooRealVar>(    "mean",    "mean", min,        50 );
  _lambda =
    std::make_unique<RooRealVar>(  "lambda",  "lambda", 0.0,       0.2 );
  _alpha =
    std::make_unique<RooRealVar>(   "alpha",   "alpha", 0,         0.2 );
  _beta =
    std::make_unique<RooRealVar>(    "beta",    "beta", 10,      20000 );
  _dcfrac =
    std::make_unique<RooRealVar>(  "dcfrac",  "dcfrac", 0,         0.4 );
  _eps =
    std::make_unique<RooRealVar>(     "eps",     "eps", 1e-5,     1e-1 );

  _pdf = std::make_unique<SiPMPdf>( "pdf",
                                    "pdf",
                                    *_x,
                                    *_ped,
                                    *_gain,
                                    *_s0,
                                    *_s1,
                                    *_mean,
                                    *_lambda,
                                    *_alpha,
                                    *_beta,
                                    *_dcfrac,
                                    *_eps );

  // Estimation options
  ignore_ped_est     = false;
  ignore_gain_est    = false;
  ignore_s0_est      = false;
  ignore_s1_est      = false;
  ignore_mean_est    = false;
  ignore_lambda_est  = false;
  _est_minpeak       = 0.05;
  _est_gausswindow   = 3;
  _est_maxgausswidth = 6;

  // operation parameters
  _intwindow = 90;
  _sipmtime  = 30;
  _sipmtype  = "Generic SiPM";
  _lumitype  = "Generic Lightsource";
  _biasv     = "Generic nominal bias";
}


SiPMLowLightFit::SiPMLowLightFit()
{
  set_all_defaults();
}


SiPMLowLightFit::SiPMLowLightFit( const std::string& config )
{
  set_all_defaults();

  usr::ArgumentExtender args;
  args.AddOptions( DataArguments() );
  args.AddOptions( FitArguments() );
  args.ParseFile( config );

  UpdateSettings( args );
}


usr::po::options_description
SiPMLowLightFit::DataArguments()
{
  usr::po::options_description desc(
    "Options for parsing the data file for low-light SiPM spectrum model" );
  desc.add_options()
    ( "input", usr::po::reqvalue<std::string>(), "Input data file" )
    ( "waveform",
    usr::po::reqvalue<bool>(),
    "Whether in input data is a raw waveform (True) are integrated (False)" )
    ( "binwidth",
    usr::po::value<double>(),
    "The bin width to use for binned data" )
    ( "maxarea", usr::po::value<double>(), "Maximum area for perform fit on" )

  /// Options unique to the waveform input format.
    ( "intstart",
    usr::po::value<unsigned>(),
    "Time slice to start the integration window" )
    ( "intstop",
    usr::po::value<unsigned>(),
    "Time slice to stop the integration window" )
    ( "pedstart",
    usr::po::value<unsigned>(),
    "Time slice to start the pedestal calculation, ignore to skip pedestal "
    "subtraction" )
    ( "pedstop",
    usr::po::value<unsigned>(),
    "Time slice to stop the pedestal calculation, ignore to skip pedestal "
    "subtraction" )
    ( "pedrms",
    usr::po::value<double>(),
    "Maximum RMS of values within the pedestal window, event is discarded if "
    "this value is surpassed" )
  ;

  return desc;
}


usr::po::options_description
SiPMLowLightFit::FitArguments()
{
  usr::po::options_description desc(
    "Options for fitting parameters adjustments."
    "For the options corresponding to the fit parameter names, the argument "
    "takes either 1, 2 or 3 doubles:\n"
    " - 1 means that the variable should be fixed at the specified values\n"
    " - 2 means that the variable range should span the specified values\n"
    " - 3 means that the variable should have some estimated value + the "
    "specified range, the estimated value will overide any sort of automated"
    " estimations methods" );
  desc.add_options()
    ( "ped",     usr::po::multivalue<double>(), "pedestal peak central value" )
    ( "gain",    usr::po::multivalue<double>(), "1-p.e gain factor" )
    ( "s0",      usr::po::multivalue<double>(), "Electrical noise" )
    ( "s1",      usr::po::multivalue<double>(), "pixel gain variation" )
    ( "mean",    usr::po::multivalue<double>(), "Mean number of p.es" )
    ( "lambda",  usr::po::multivalue<double>(), "Crosstalk proability" )
    ( "alpha",   usr::po::multivalue<double>(), "Afterpulsing probability" )
    ( "beta", usr::po::multivalue<double>(), "Afterpulting timescale factor" )
    ( "dcfrac",  usr::po::multivalue<double>(), "Dark current probability" )
    ( "epsilon",
    usr::po::multivalue<double>(),
    "Resolution factor to be used for the dark current curve" )
  ;

  return desc;
}


usr::po::options_description
SiPMLowLightFit::OperationArguments()
{
  usr::po::options_description desc(
    "Options of labeling the SiPM operation parameters that the fit have "
    "difficulty figuring out independently" );

  desc.add_options()
    ( "intwindow",
    usr::po::value<double>(),
    "The integration window required for waveform sums, this will overwrite "
    "intstart/intstop results used for waveform formats (units: ns)" )
    ( "sipmtime",
    usr::po::value<double>(),
    "SiPM recovery time, this calculation routines would not include functions "
    "to calculated this. (units: ns)" )
    ( "lumitype",
    usr::po::value<std::string>(),
    "Luminosity type display on plots" )
    ( "sipmtype",
    usr::po::value<std::string>(),
    "SiPM type to display on plots" )
    ( "biasv",
    usr::po::value<std::string>(),
    "String on plot used to indicated the bias condition" )
  ;

  return desc;
}


usr::po::options_description
SiPMLowLightFit::EstArguments()
{
  usr::po::options_description desc(
    "Options related to setting the parameters used for parameter estimation process" );
  desc.add_options()
    ( "estminpeak",
    usr::po::value<double>(),
    "The minimum fraction a peak value in the histogram have relative to the maximum bin value that it can be considered a peak" )
    ( "estgausswindow",
    usr::po::value<int>(),
    "Number of bins that that local guassian peak should follow" )
    ( "estmaxgausswidth",
    usr::po::value<int>(),
    "Maximum number of bins that the gaussian width can be before it is discarded as a primary peak candidate" )
  ;
  return desc;
}


void
SiPMLowLightFit::UpdateSettings( const usr::ArgumentExtender& args )
{
  // Updating the data parsing arguments
  _inputfile = args.ArgOpt<std::string>( "input",  _inputfile );
  _waveform  = args.ArgOpt<bool>(     "waveform",  _waveform  );
  _binwidth  = args.ArgOpt<double>(   "binwidth",  _binwidth  );
  _intstart  = args.ArgOpt<unsigned>( "intstart",  _intstart  );
  _intstop   = args.ArgOpt<unsigned>( "intstop",   _intstop   );
  _pedstart  = args.ArgOpt<unsigned>( "pedstart",  _pedstart  );
  _pedstop   = args.ArgOpt<unsigned>( "pedstop",   _pedstop   );
  _maxarea   = args.ArgOpt<double>(   "maxarea",   _maxarea   );

  // Updating the fitting arguments
  auto f1 = []( RooRealVar& x, double val ){
              x = val;
              x.setConstant( true );
            };
  auto f2 = []( RooRealVar& x, double min, double max ){
              x.setRange( min, max );
              x.setConstant( false );
            };
  auto f3 = []( RooRealVar& x, double val, double min, double max ){
              x = val;
              x.setRange( min, max );
              x.setConstant( false );
            };
  auto update_arg =
    [&args, &f1, &f2, &f3]( RooRealVar& x, const std::string& var ){
      if( args.CheckArg( var ) ){
        const auto vec = args.ArgList<double>( var );
        if( vec.size() == 1 ){
          f1( x, vec[0] );
        } else if( vec.size() == 2 ){
          f2( x, vec[0], vec[1] );
        } else if( vec.size() == 3 ){
          f3( x, vec[0], vec[1], vec[2] );
        }
      }
    };

  update_arg( *_ped,    "ped"     );
  update_arg( *_gain,   "gain"    );
  update_arg( *_s0,     "s0"      );
  update_arg( *_s1,     "s1"      );
  update_arg( *_mean,   "mean"    );
  update_arg( *_lambda, "lambda"  );
  update_arg( *_alpha,  "alpha"   );
  update_arg( *_beta,   "beta"    );
  update_arg( *_dcfrac, "dcfrac"  );
  update_arg( *_eps,    "epsilon" );

  // Option parsing for estimation related variables
  auto lock = [&args]( bool& ignore, const std::string& var ){
                if( args.CheckArg( var ) ){
                  const auto vec = args.ArgList<double>( var );
                  if( vec.size()  == 1 || vec.size() == 3 ){
                    ignore = true;
                  } else {
                    ignore = false;
                  }
                } else {
                  ignore = false;
                }
              };

  lock( ignore_ped_est,    "ped"    );
  lock( ignore_gain_est,   "gain"   );
  lock( ignore_s0_est,     "s0"     );
  lock( ignore_s1_est,     "s1"     );
  lock( ignore_mean_est,   "mean"   );
  lock( ignore_lambda_est, "lambda" );

  _est_minpeak     = args.ArgOpt<double>( "estminpeak", _est_minpeak       );
  _est_gausswindow =
    args.ArgOpt<double>( "estgausswindow", _est_gausswindow   );
  _est_maxgausswidth = args.ArgOpt<double>( "estmaxgausswidth",
                                            _est_maxgausswidth );

  // Operation parameters related variables
  _intwindow = args.ArgOpt<double>( "intwindow", _intwindow );
  _sipmtime  = args.ArgOpt<double>( "sipmtime",  _sipmtime );

  _sipmtype = args.ArgOpt<std::string>( "sipmtype", _sipmtype );
  _lumitype = args.ArgOpt<std::string>( "lumitype", _lumitype );
  _biasv    = args.ArgOpt<std::string>( "biasv",    _biasv    );
}


// Procedues to perform a low light fit. Sementing to pseudo-interactively.
void
SiPMLowLightFit::MakeBinnedData()
{
  if( _waveform ){
    make_array_from_waveform();
  } else {
    make_array_from_sum();
  }

  std::sort( _arealist.begin(), _arealist.end() );

  // Parsing for converting into data.
  const double xmin = usr::RoundDown( _arealist.front(), _binwidth );
  const double xmax = usr::RoundUp( std::min(
                                      _arealist.back(),
                                      _maxarea ),
                                    _binwidth );
  const double nbins = ( xmax-xmin ) / _binwidth;

  x().setRange( xmin, xmax );
  x().setBins( nbins );

  _data.reset( new RooDataHist( "data", "data", RooArgList( x() ) ) );

  for( const auto a : _arealist ){
    if( a < xmax ){
      x() = a;
      _data->add( RooArgList( x() ) );
    }
  }
}


void
SiPMLowLightFit::make_array_from_waveform()
{
  // Using the common SiPM waveformat class to get wave from.
  WaveFormat wformat( _inputfile );

  // Looping over the value,
  for( unsigned i = 0; i < wformat.NWaveforms(); ++i ){
    if( wformat.PedRMS( i, _pedstart, _pedstop ) > _pedrms ){
      continue;
    }
    const double a = wformat.WaveformSum( i,
                                          _intstart,
                                          _intstop,
                                          _pedstart,
                                          _pedstop );
    _arealist.push_back( a );
  }

  // Additional parsing required for plotting
  const unsigned start = std::min( _intstart, wformat.NSamples() );
  const unsigned stop  = std::min( _intstop, wformat.NSamples() );
  _intwindow = wformat.Time() * ( stop-start );
}


void
SiPMLowLightFit::make_array_from_sum()
{
  // Using the Common format class for get summed results
  StdFormat  sformat( _inputfile );
  const auto results = sformat.DataAll();
  _arealist.insert( _arealist.begin(), results.begin(), results.end() );
}


void
SiPMLowLightFit::RunFit()
{
  usr::ConvergeFitPDFToData( *_pdf, *_data, usr::MaxFitIteration( 3 ) );

  // Limiting to 3 to save runtime.
}


usr::Measurement
SiPMLowLightFit::Pedestal() const
{
  return usr::Measurement( ped().getVal(), ped().getError() );
}


usr::Measurement
SiPMLowLightFit::Gain() const
{
  return usr::Measurement( gain().getVal(), gain().getError() );
}


usr::Measurement
SiPMLowLightFit::CommonNoise() const
{
  return usr::Measurement( s0().getVal(), s0().getError() );
}


usr::Measurement
SiPMLowLightFit::PixelNoise() const
{
  return usr::Measurement( s1().getVal(), s1().getError() );
}


usr::Measurement
SiPMLowLightFit::MeanPhotons() const
{
  return usr::Measurement( mean().getVal(), mean().getError() );
}


usr::Measurement
SiPMLowLightFit::MeanPhotonsFromFrac() const
{
  // Getting the mean photon using the fractional method:
  const double den = _arealist.size();
  const double num = std::count_if( _arealist.begin(),
                                    _arealist.end(), [this]( double x )->bool {
    return x < ( this->ped().getVal()+this->gain().getVal() / 2 );
  } );

  const double frac = num / den;
  const double err  = TMath::Sqrt( frac * ( 1-frac ) / den );

  const double cen = -TMath::Log( frac );
  const double unc = fabs( -TMath::Log( frac )+TMath::Log( frac+err ) );

  return usr::Measurement( cen, unc );
}


usr::Measurement
SiPMLowLightFit::ProbCrosstalk() const
{
  const double lval      = lambda().getVal();
  const double lerr      = lambda().getError();
  const double xtalk_cen = 1-TMath::Exp( -lval  );
  const double xtalk_hi  = 1-TMath::Exp( -( lval+lerr ) );
  const double xtalk_lo  = 1-TMath::Exp( -( lval-lerr ) );
  const double xtalk_err = ( xtalk_hi-xtalk_lo ) / 2;
  return usr::Measurement( xtalk_cen, xtalk_err );
}


usr::Measurement
SiPMLowLightFit::ProbAfterpulse() const
{
  return usr::Measurement( alpha().getVal(), alpha().getError() );
}


usr::Measurement
SiPMLowLightFit::AfterpulseTimeNS() const
{
  const usr::Measurement mb( beta().getVal(), beta().getError() );
  return _sipmtime * mb / Gain();
}


usr::Measurement
SiPMLowLightFit::DarkcurrentTimeNS() const
{
  const usr::Measurement prob( dcfrac().getVal(), dcfrac().getError() );
  return _intwindow / prob;
}


usr::Measurement
SiPMLowLightFit::ExcessNoiseFactor( const usr::Measurement& n ) const
{
  const double stddev = usr::StdDev( _arealist );
  const double mean   = usr::Mean( _arealist )-ped().getVal();
  const double enf    = n * ( stddev / mean ) * ( stddev / mean );
  return usr::Measurement( enf, enf * n.RelAvgError() );
}
