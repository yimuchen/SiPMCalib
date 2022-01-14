#include "SiPMCalib/Common/interface/StdFormat.hpp"
#include "SiPMCalib/InvSqCalc/interface/InvSqFunc.hpp"
#include "SiPMCalib/SiPMCalc/interface/NonLinearModel.hpp"
#include "SiPMCalib/SiPMCalc/interface/SiPMNonLinearFit.hpp"

#include "UserUtils/Common/interface/Maths.hpp"
#include "UserUtils/Common/interface/STLUtils.hpp"
#include "UserUtils/MathUtils/interface/RootMathTools.hpp"
#include "UserUtils/PlotUtils/interface/Flat2DCanvas.hpp"
#include "UserUtils/PlotUtils/interface/PlotCommon.hpp"
#include "UserUtils/PlotUtils/interface/Ratio1DCanvas.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

/**
 * @class SiPMNonLinearFit
 * @ingroup SiPMCalc
 *
 * @brief Fitting Routine class for SiPM nonlinearity.
 *
 * @details The class takes the function of the input of a the zscan command
 * from the [control program][control-zscan], calculates the number of estimated
 * incoming photons for each of the data rows collected via the inverse square
 * estimate and a corrector file for the LED bias settings, to get a graph for
 * the readout (in units of fired pixels) vs the estimated number of incoming
 * photons. The output of this file can be used to investigate the nonlinearity
 * between the light input and the entire signal output chain (SiPM and readout
 * amplification/digitization chain) in the format of presentation ready plots.
 *
 * As in the final analysis chain when we are running bulk analysis over many,
 * similarly parameterized SiPMs, certain inputs to the nonlinearity fit will
 * only be available after addition analysis routines of files other than the
 * output of the `zscan` command , these we leave as function input arguments.
 * All other arugments shall be defined as program options so that they can
 * either be specified with a UNIX-cfg file or via command line arguments. For
 * details of these arguments, see the following methods:
 *
 * - `FitArguments()` : For augmenting how the data within the zscan data should
 *   be manipulated.
 * - `OutputArguments()` : For modifying the output plot display string.
 *
 * This class the responsible main for morphing the data according the inverse
 * square law and some corrections for the LED bias to get the number of
 * estimated photons for data row collected. A typically analysis routine will
 * run the following functions in order:
 *
 * - `ReadFiles()`: Read the data collected via a `zscan` command, and some LED
 *   bias correction file.
 * - `MakeLinearGraph()`: Getting the luminosity multiplier as a function of the
 *   data collection z cooridnates, assuming that for some small LED bias range
 *   the SiPM is sufficient linear and the deviation from a perfect point light
 *   source for the given data collection range is small.
 * - `MakeNonLinearGraph()`: using the the geometric luminosity multiplier and
 *   the LED bias corrector file, make the graph of a function of estimated
 *   incoming photons vs the readout value.
 * - `RunNonLinearFit()`: routine for calling the fitter and fitting the data to
 *   the results.
 *
 * Look at each of the detailed documentation of the functions above to see how
 * each step is done. In addition to the analysis functions, the various `Plot`
 * methods can be called to get presentation ready-plots at each stage of the
 * analysis chain. For an example of how the class can be used, please refer to
 * the FitNonLinear.cc file.
 *
 * [control-zscan]:
 * https://umdcms.github.io/SiPMCalibControl/classctlcmd_1_1motioncmd_1_1zscan.html
 *
 */


SiPMNonLinearFit::SiPMNonLinearFit()
{
  InitDefaults();
}


SiPMNonLinearFit::SiPMNonLinearFit( const std::string& config )
{
  InitDefaults();
}


/**
 * @brief Arguments for updating how the linear and nonlinear part of the fit
 * should be augmented.
 *
 * The allow options are:
 *
 * - `linpmin`, `linpmax`: These are the bias value limits of where the SiPM
 *   output is expected to be sufficiently linear and the readout should follow a
 *   near-perfect inverse square with the data collection z value (with z-offset
 *   and pedestal of course.) Notice that these values should effective bracket a
 *   single pwm duty cycle settings used in the zscan command, except that the
 *   the bias values would have a finite range due to DC drifts in the system.
 * - `powerz`: In the case the a bias correction file is not provided. one can be
 *   generated on the fly using the readout values of data points closest to the
 *   z value provided by this option. Notice that because this method lacks
 *   proper pedestal subtraction, the output of this method can potentially be
 *   wrong. Notice this will be ignored regardless of input value if a bias
 *   corrector file was provided for reading.
 * - `npixels`: The number of pixels to use in the fit. This can either be a
 *   fixed value (in which case, simply use this options with a single input) or
 *   can be used to indicate the range of the fit.
 */
usr::po::options_description
SiPMNonLinearFit::FitArguments()
{
  usr::po::options_description desc( "Fitting arguments for non-linear fit" );
  desc.add_options()
    ( "linpmin",
    usr::po::value<double>(),
    "Minimum bias voltage for inverse square fit" )
    ( "linpmax",
    usr::po::value<double>(),
    "Maximum bias voltage for the inverse square fit" )
    ( "linzmin",
    usr::po::value<double>(),
    "Minimum coordinate z value to use for the inverse square fit" )
    ( "powerz",
    usr::po::value<double>(),
    "z value to use for power level look up, (ignored if a correction file is provided)" )
    ( "npixels",
    usr::po::multivalue<double>(),
    "The number of pixels in the SiPM of interest" )
  ;
  return desc;
}


/**
 * @brief Options for modifying the output plots
 *
 * @details The following options are implemented:
 * - `units` The string to be placed in the units of the y axis in the non-linear
 *   plot.
 * - `setup` The string to be placed in the top right of the pad (where the
 *   luminosity data is usually placed in CMS plots)
 * - `sipmmodel` The string representing the SiPM model (written underneath the
 *   "HGCAL" standard line at the top left of the canvase)
 * - `sipmid` The string representing the unique ID of the SiPM (written
 *   underneath the SiPM model string.)
 *
 */
usr::po::options_description
SiPMNonLinearFit::OutputArguments()
{
  usr::po::options_description desc(
    "Additional plotting arguments for non-linear fit" );
  desc.add_options()
    ( "units",
    usr::po::value<std::string>(),
    "Units of the y axis in final plot" )
    ( "setup",
    usr::po::value<std::string>(),
    "Setup string to display on the top of the plot" )
    ( "sipmmodel", usr::po::value<std::string>(), "SiPM model display" )
    ( "sipmid", usr::po::value<std::string>(), "SiPM model index" )
  ;
  return desc;
}


/**
 * @brief  Initializing the setting parameters
 *
 */
void
SiPMNonLinearFit::InitDefaults()
{
  _lin_pmax  = 0;
  _lin_pmin  = 0;
  _lin_zmin  = 0;
  _power_z   = 500;
  _pixel_min = 500;
  _pixel_max = 7000;

  _readout_units = "[A.U.]";
  _setup         = "UMD Gantry";
  _sipm_model    = "HDR2 (7500 pix.)";
  _sipm_id       = "";
}


/**
 * @brief Updating the fitting routine parameters according to the input
 * arguments.
 */
void
SiPMNonLinearFit::UpdateSettings( const usr::ArgumentExtender& args )
{
  _lin_pmin = args.ArgOpt<double>( "linpmin", _lin_pmin );
  _lin_pmax = args.ArgOpt<double>( "linpmax", _lin_pmax );
  _lin_zmin = args.ArgOpt<double>( "linzmin", _lin_pmin );
  _power_z  = args.ArgOpt<double>( "powerz", _power_z );

  if( args.CheckArg( "npixels" ) ){
    const auto npix = args.ArgList<double>( "npixels" );
    if( npix.size() == 1 ){
      _pixel_min = _pixel_max = npix.at( 0 );
    } else if( npix.size() == 2 ){
      _pixel_min = usr::GetMinimum( npix );
      _pixel_max = usr::GetMaximum( npix );
    } // Do nothing
  }

  _readout_units = args.ArgOpt<std::string>( "units", _readout_units );
  _setup         = args.ArgOpt<std::string>( "setup", _setup );
  _sipm_model    = args.ArgOpt<std::string>( "sipmmodel", _sipm_model );
  _sipm_id       = args.ArgOpt<std::string>( "sipmid", _sipm_id );
}


/**
 * @brief Reading the raw data file and the data correction file, storing results
 * into the class.
 *
 * At this point no data manipulation is used. Only data correction is made. For
 * details on how the corrector file is handled, see the MakeInterpolator3D()
 * method.
 */
void
SiPMNonLinearFit::ReadFiles( const std::string& zscan, const std::string& corr )
{
  // Making the common data container for the raw data.
  _raw_data.reset( new StdFormat( zscan ));

  // Clearing and updating the exiting look up class
  _lookup_1d.reset( nullptr );
  _lookup_3d.reset( nullptr );

  if( corr == "" ){
    MakeInterpolator1D();
  } else {
    MakeInterpolator3D( corr );
  }
}


/**
 * @brief In the case that the external correction table is not available. Create
 * a simplified lookup table using correcting data readout values.
 *
 * The look-up table here is generated using the data collection values at the
 * specified z value, and will only be a correction value for the bias voltage.
 * (temperature drift correction will be omitted). Notice that this method
 * effectively ignores any sort of imperfect pedestal subtraction results. and
 * thus can be inaccurate, but works well for a quick and dirty evaluation of
 * results.
 */
void
SiPMNonLinearFit::MakeInterpolator1D()
{
  const auto   z         = _raw_data->Z();
  const double closest_z = *min_element(
    z.begin(),
    z.end(), [this]( const double l,  const double r )->bool {
    return fabs( l-this->_power_z ) < fabs( r-this->_power_z );
  } );
  auto at_z = [closest_z]( const StdFormat::RowFormat& row )->bool {
                return row.z == closest_z;
              };
  std::vector<double> bias    = _raw_data->Bias( at_z );
  std::vector<double> readout = _raw_data->DataCol( 0, at_z );

  bias.insert( bias.begin(), -10000000 );
  readout.insert( readout.begin(), readout.front() );
  bias.push_back( 100000000 );
  readout.push_back( readout.back() );

  // Making the 1D interpolator
  _lookup_1d.reset(
    new ROOT::Math::Interpolator( bias,
                                  readout,
                                  ROOT::Math::Interpolation::kLINEAR ));
}


/**
 * @brief Creating the data correction lookup table from an external correction
 * file.
 *
 * The corrector file is a file containing numbers in 2 columns, the first
 * columns is the LED bias value, the second number is the relative luminosity
 * output (in arbitrary units) of the corresponding LED bias value (assuming 0
 * temperature drift for now). All values are then used to construct the
 * interpolator.
 */
void
SiPMNonLinearFit::MakeInterpolator3D( const std::string& corr )
{
  std::ifstream infile( corr );
  std::string   line;

  std::vector<double> bias;
  std::vector<double> readout;

  while( std::getline( infile, line ) ){
    std::istringstream linestream( line );
    double             bias_value;
    double             readout_value;
    linestream >> bias_value >> readout_value;
    bias.push_back( bias_value );
    readout.push_back( readout_value );
  }

  _lookup_1d.reset(
    new ROOT::Math::Interpolator( bias,
                                  readout,
                                  ROOT::Math::Interpolation::kLINEAR ));
}


/**
 * @brief Making the linear graph as well as fitting to the inverse square law
 * for this function. Since the function will be plotted
 *
 * TODO: proper evaluation of data corrections.
 */
void
SiPMNonLinearFit::MakeLinearGraph()
{
  auto lin_power = [this]( const StdFormat::RowFormat& row )->bool {
                     return row.bias >= this->_lin_pmin &&
                            row.bias <= this->_lin_pmax;
                   };
  const std::vector<double> z    = _raw_data->Z( lin_power );
  const std::vector<double> lumi = _raw_data->DataCol( 0, lin_power );
  const std::vector<double> unc  = _raw_data->DataCol( 1, lin_power );
  const std::vector<double> bias = _raw_data->Bias( lin_power );
  const std::vector<double> stmp = _raw_data->SiPMTemp( lin_power );
  const std::vector<double> ptmp = _raw_data->LedTemp( lin_power );
  const std::vector<double> zero( z.size(), 0 );

  assert( z.size() > 0  );

  // Making the TGraph object
  _lin_data = TGraphErrors( z.size(),
                            z.data(),
                            lumi.data(),
                            zero.data(),
                            unc.data() );

  const double zmin = usr::GetMinimum( z );
  const double zmax = usr::GetMaximum( z );
  const double lmin = usr::GetMinimum( lumi );
  const double lmax = usr::GetMaximum( lumi );

  _lin_func = TF1( usr::RandomString( 12 ).c_str(), InvSq_Z, zmin, zmax, 4 );

  // Setting function limits
  _lin_func.SetParLimits( 1, 0.0, 1e2 ); // Horizontal offset
  _lin_func.SetParLimits( 3, 0.0, 1e2 ); // Pedestal

  // Making estimates
  _lin_func.SetParameter( 0, 0.0 );
  _lin_func.SetParameter( 1, 0.0 );
  _lin_func.SetParameter( 2, lmax );
  _lin_func.SetParameter( 3, lmin );

  // _lin_func.SetParLimits( 1, 0, 5.0 );
  usr::fit::FitGraph( _lin_data,
                      _lin_func,
                      usr::fit::GraphXRange( std::max( zmin, _lin_zmin ),
                                             zmax ) );
}


/**
 * @brief Creating the graph for the number of expected photons (x axis) vs the
 * number of fired pixels as seen in the readout given the reference point for a
 * number of given photons.
 *
 * The function expects the following inputs, the ref_n, ref_z, and ref_bias
 * refers to that for this SiPM, we expect ref_n effective photons for when the
 * SiPM is at position ref_z and LED bias value ref_bias. As ref_n is going to
 * be obtain by a low-light spectral fit (external to the data collected via the
 * zscan data), these values need to be provided externally as function inputs.
 * As bias voltages can drift between the zscan command the data collection for
 * a low light spectrum, the closest value data row will be used as the
 * "reference row".
 *
 * Once the reference row has been determined, the photon multiplier according
 * to the z position (via the fit to inverse square) and the bias/temperature
 * drifts via the look up table, can be used to morph the data according to the
 * two multipliers to get the number of incoming photons vs the readout value.
 *
 * In addition to the 3 inputs above, the gain and the pedestal should also be
 * given to this function to further reduce the readout value to be equivalent
 * to the average number of pixels fired. One can of course, opt to leave these
 * two values at 1 and 0 if one wants to compare the raw readout value v.s. the
 * number of incoming photons.
 *
 * Notice here we will still be saving the bias voltage for each of the data
 * points, this is useful should the user need to debug the corrector file.
 */
void
SiPMNonLinearFit::MakeNonLinearGraph( const double ref_n,
                                      const double ref_z,
                                      const double ref_bias,
                                      const double gain,
                                      const double ped  )
{
  // Finding the new reference point data power
  const auto ref_row = *std::min_element(
    _raw_data->begin(),
    _raw_data->end(),
    [ref_z, ref_bias]( const StdFormat::RowFormat& l,
                       const StdFormat::RowFormat& r )->bool {
    if( fabs( l.z-ref_z ) != fabs( r.z-ref_z ) ){
      return fabs( l.z-ref_z ) < fabs( r.z-ref_z );
    } else {
      return fabs( l.bias-ref_bias ) < fabs( r.bias-ref_bias );
    }
  } );

  std::vector<double> n_in;
  std::vector<double> n_out;
  std::vector<double> in_unc;
  std::vector<double> out_unc;
  std::vector<double> bias;
  std::vector<double> zeros;

  // Looping over all data

  // TODO: the positional estimator should also include uncertainties. Find a
  // good way to properly propagate this.
  for( const auto row : *_raw_data ){
    const double z_mult = GetZPosMultiplier( row, ref_row, ped );
    const double l_mult = GetLookUpMultiplier( row,  ref_row );

    n_in.push_back( ref_n * z_mult * l_mult );
    n_out.push_back( row.data.at( 0 ) / gain );
    in_unc.push_back( 0 );  // FIXME!!
    out_unc.push_back( row.data.at( 1 ) / gain );

    bias.push_back( row.bias );
    zeros.push_back( 0 );
  }

  _nl_data = TGraph2DErrors( n_in.size(),
                             n_in.data(),
                             n_out.data(),
                             bias.data(),
                             in_unc.data(),
                             out_unc.data(),
                             zeros.data() );
}


/**
 * @brief Running the non linear fit.
 *
 */
void
SiPMNonLinearFit::RunNonLinearFit()
{
  const double xmin = usr::plt::GetXmin( _nl_data );
  const double xmax = usr::plt::GetXmax( _nl_data );

  _nl_func = TF1( usr::RandomString( 12 ).c_str(), NLOModel, xmin, xmax, 5 );
  _nl_func.SetParameter( 0, 1.0 ); // Gain
  _nl_func.SetParLimits( 0, 0, 1000 );

  // _nl_func.FixParameter( 2, 1.0 );
  _nl_func.SetParameter( 2, 0.0 );  // Pedestal
  _nl_func.SetParLimits( 2, 0.0, 10 );

  _nl_func.SetParameter( 3, 0.0 ); // Alpha
  _nl_func.SetParLimits( 3, 0, 1.0 );
  _nl_func.SetParameter( 4, 0.0 ); // Beta
  _nl_func.SetParLimits( 4, 0, 10.0 );

  if( _pixel_min == _pixel_max ){
    _nl_func.FixParameter( 1, _pixel_min );
  } else {
    _nl_func.SetParameter( 1, _pixel_min ); // Total Number of pixels
    _nl_func.SetParLimits( 1, _pixel_min, _pixel_max );
  }

  TGraphErrors _nl_g = TGraphErrors( _nl_data.GetN(),
                                     _nl_data.GetX(),
                                     _nl_data.GetY(),
                                     _nl_data.GetEX(),
                                     _nl_data.GetEY() );

  // Running the simple
  _nl_fit = usr::fit::FitGraph( _nl_g, _nl_func );

  // Always printing the results for the non-linear fit.
  _nl_fit.Print();
}


/**
 * @brief Generating the plot for the inverse square fit.
 *
 * To make the results easier to interpret. We are going to shift the plotting
 * results. Such that the z offset is 0.
 *
 */
void
SiPMNonLinearFit::PlotLinearity( const std::string& outfile )
{
  usr::plt::Ratio1DCanvas c;

  TF1                 lfunc = _lin_func;
  TGraphErrors        ldata = _lin_data;
  static const double z0    = _lin_func.GetParameter( 0 );

  for( int i = 0 ; i <  ldata.GetN() ; ++i ){
    ldata.GetX()[i] = ldata.GetX()[i]-z0;
  }

  const double zmin = usr::plt::GetXmin( ldata );
  const double zmax = usr::plt::GetXmax( ldata );
  lfunc.SetParameter( 0, 0.0 );
  lfunc.SetRange( zmin, zmax );
  auto fittemp =
    usr::fit::FitGraph( ldata,
                        lfunc,
                        usr::fit::GraphXRange( zmin+_lin_zmin,
                                               zmax ) );

  c.PlotGraph(
    ldata,
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
    usr::plt::MarkerSize( 0.1 ),
    usr::plt::ExtendXRange( false ) );
  auto& fitg = c.PlotFunc(
    lfunc,
    usr::plt::PlotType( usr::plt::fittedfunc ),
    usr::plt::EntryText(  "Fitted Data" ),
    usr::plt::VisualizeError( fittemp ),
    usr::plt::FillColor( usr::plt::col::cyan ),
    usr::plt::TrackY( usr::plt::tracky::max ),
    usr::plt::LineColor( usr::plt::col::blue ));
  auto& datag = c.PlotGraph(
    ldata,
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::EntryText( "Readout" ),
    usr::plt::MarkerColor(  usr::plt::col::black ),
    usr::plt::MarkerStyle(  usr::plt::sty::mkrcircle ),
    usr::plt::MarkerSize( 0.2 ),
    usr::plt::LineColor( usr::plt::col::gray, 1.0 ) );

  c.TopPad().DrawVLine( zmin+_lin_zmin,
                        usr::plt::LineStyle( usr::plt::sty::lindashed ),
                        usr::plt::LineColor( usr::plt::col::gray, 0.5 ) );

  c.PlotScale(
    fitg,
    fitg,
    usr::plt::PlotType( usr::plt::fittedfunc ),
    usr::plt::FillColor( usr::plt::col::cyan ) );
  c.PlotScale(
    datag,
    fitg,
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::MarkerSize( 0.2 ) );
  c.BottomPad().DrawVLine( zmin+_lin_zmin,
                           usr::plt::LineStyle( usr::plt::sty::lindashed ),
                           usr::plt::LineColor( usr::plt::col::gray, 0.5 ) );

  c.TopPad().Yaxis().SetTitle( "Luminosity [mV-ns]" );
  c.BottomPad().Xaxis().SetTitle( "Gantry z + z_{0} [mm]" );
  c.BottomPad().Yaxis().SetTitle( "Data/Fit" );

  c.TopPad().SetLogy( kTRUE );
  c.TopPad().SetLogx( kTRUE );
  c.TopPad().SetYaxisMin( usr::plt::GetYmin( _lin_data ) / 30  );
  c.BottomPad().SetLogx( kTRUE );

  c.DrawLuminosity( "LED setup" );
  c.DrawCMSLabel( "Preliminary", "HGCal" );
  c.TopPad()
  .SetTextCursor( c.TopPad().InnerTextLeft(), 0.45 )
  .WriteLine( "#frac{L_{0}(z+z_{0})}{(x_{0}^{2} + (z+z_{0})^{2})^{3/2}} + Ped" );
  c.TopPad()
  .SetTextCursor( c.TopPad().InnerTextLeft(), 0.25 )
  .WriteLine( usr::fstr( "Fit Ped = %.2lf_{#pm%.3lf}",
                         _lin_func.GetParameter( 3 ),
                         _lin_func.GetParError( 3 ) ) )
  .WriteLine( usr::fstr( "Fit z_{0} = %.2lf_{#pm%.3lf}",
                         -_lin_func.GetParameter( 0 ),
                         _lin_func.GetParError( 0 ) ) )
  .WriteLine( usr::fstr( "Fit x_{0} = %.2lf_{#pm%.3lf}",
                         _lin_func.GetParameter( 1 ),
                         _lin_func.GetParError( 1 ) ) )
  .WriteLine( usr::fstr( "#chi^{2}/D.o.F = %.2lf/%d",
                         _lin_fit.Chi2(),
                         _lin_fit.Ndf() ) );

  c.SaveAsPDF( outfile );
}


/**
 * @brief  Plotting the original data plots with the various bias voltages
 * represented using different color points.
 *
 * Notice that this functions needs to be ran *after* the linear fit has been
 * ran, as we will need the z off set subtraction to avoid unnecessary confusions
 * when interpreting the results of the plots.
 */
void
SiPMNonLinearFit::PlotOriginal( const std::string& outfile )
{
  usr::plt::Flat2DCanvas c;

  auto zval    = _raw_data->Z();
  auto readout = _raw_data->DataCol( 0 );
  auto bias    = _raw_data->Bias();

  std::transform( zval.begin(),
                  zval.end(),
                  zval.begin(), [this]( const double z )->double {
    return z-this->_lin_func.GetParameter( 0 );
  } );

  TGraph2D g( zval.size(), zval.data(), readout.data(), bias.data() );

  c.PlotColGraph( g,
                  usr::plt::MarkerSize( 0.5 ),
                  usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ) );
  c.Pad().Xaxis().SetTitle( "Gantry coordinates" );
  c.Pad().Yaxis().SetTitle( "Readout value" );
  c.Pad().Zaxis().SetTitle( "LED Bias [mV]" );
  c.Pad().DrawLuminosity( _setup );
  c.Pad().DrawCMSLabel( "", "HGCAL" );
  c.Pad().SetLogx( kTRUE );
  c.Pad().SetLogy( kTRUE );
  c.SaveAsPDF( outfile );
}


/**
 * @brief Plotting the morphed data with the bias voltage as colors.
 *
 * This can be done before the data points nonlinear fit has been performed to
 * help debug with issues with the bias corrector files.
 */
void
SiPMNonLinearFit::PlotMorphed( const std::string& outfile )
{
  usr::plt::Flat2DCanvas c;
  c.PlotColGraph( _nl_data,
                  usr::plt::MarkerSize( 0.5 ),
                  usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ) );

  c.Pad().Xaxis().SetTitle( "Estimated number of incident photons" );
  c.Pad().Yaxis().SetTitle( ( "Readout "+_readout_units ).c_str() );
  c.Pad().Zaxis().SetTitle( "LED Bias [mV]" );
  c.Pad().DrawLuminosity( _setup );
  c.Pad().DrawCMSLabel( "", "HGCAL" );
  c.Pad().SetLogx( kTRUE );
  c.Pad().SetLogy( kTRUE );
  c.SaveAsPDF( outfile );
}


/**
 * @brief Ploting the fit results of the zscan data used for inverse square
 * estimation.
 */
void
SiPMNonLinearFit::PlotNonLinearity( const std::string& outfile )
{
  usr::plt::Ratio1DCanvas c;

  TGraphErrors _nl_g = TGraphErrors( _nl_data.GetN(),
                                     _nl_data.GetX(),
                                     _nl_data.GetY(),
                                     _nl_data.GetEX(),
                                     _nl_data.GetEY() );

  auto& fitg = c.PlotFunc(
    _nl_func,
    usr::plt::PlotType( usr::plt::fittedfunc ),
    usr::plt::VisualizeError( _nl_fit ),
    usr::plt::EntryText( "Fit" ),
    usr::plt::LineColor( usr::plt::col::blue, 1.0 ),
    usr::plt::FillColor( usr::plt::col::cyan, 1.0 ),
    usr::plt::Precision( 0.001, true ),
    usr::plt::TrackY( usr::plt::tracky::max ));
  auto& datag = c.PlotGraph(
    _nl_g,
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::EntryText( "Data" ),
    usr::plt::TrackY( usr::plt::tracky::both ),
    usr::plt::MarkerColor(  usr::plt::col::black ),
    usr::plt::MarkerStyle(  usr::plt::sty::mkrcircle ),
    usr::plt::MarkerSize( 0.2 ),
    usr::plt::LineColor( usr::plt::col::black, 1.0 ),
    usr::plt::ExtendXRange( false ));

  // Plotting the ratio canvas
  c.PlotScale( fitg, fitg, usr::plt::PlotType( usr::plt::fittedfunc ) );
  c.PlotScale( datag,
               fitg,
               usr::plt::PlotType (
                 usr::plt::scatter ),
               usr::plt::MarkerSize( 0.2 ) );

  // Axis decorators
  c.TopPad().Yaxis().SetTitle( ( "Readout "+_readout_units ).c_str() );
  c.BottomPad().Xaxis().SetTitle( "Estimated number of incident photons" );
  c.BottomPad().Yaxis().SetTitle( "Data/Fit" );
  c.TopPad().SetLogy( kTRUE );
  c.TopPad().SetLogx( kTRUE );
  c.BottomPad().SetLogx( kTRUE );

  /// Adding additional decorator string
  c.TopPad().DrawLuminosity( _setup );
  c.TopPad().DrawCMSLabel( "", "HGCAL" );
  c.TopPad().WriteLine( _sipm_model ).WriteLine( _sipm_id );

  // Saving output file
  c.SaveAsPDF( outfile );
}


/**
 * @brief Getting the photon multiplier factor based on the z coordinates
 * relative to some reference point.
 *
 * Ths function will only ever the use z cooridnates of the data rows of
 * interest. In addition, the use can provide a custom pedestal value which will
 * make it so that the pedestal value evaluated from the inverse square fit is
 * ignored.
 */
double
SiPMNonLinearFit::GetZPosMultiplier( const StdFormat::RowFormat& main,
                                     const StdFormat::RowFormat& ref,
                                     const double                ped ) const
{
  const double op_ped =  ped == 0 ?
                        _lin_func.GetParameter( 3 ) :
                        ped;

  auto eval_row = [this, op_ped]( const StdFormat::RowFormat& r )->double {
                    return ( this->_lin_func.Eval( r.z ))-op_ped;
                  };
  return eval_row( main ) / eval_row( ref );
}


/**
 * @brief Getting the multiplier factor that should be applied given the bias
 * voltage and the temperature of the row compared with the reference data
 * point.
 *
 * This automatically detects whether the look up for the 3D interpolator
 * exists for bias and temperature lookup; otherwise simply use the 1D bias
 * interpolator constructed from the data scan.
 */
double
SiPMNonLinearFit::GetLookUpMultiplier( const StdFormat::RowFormat& main,
                                       const StdFormat::RowFormat& ref ) const
{
  if( _lookup_3d != nullptr ){
    //TODO: implement the temperature compensation version of this look up table.
    return 0;
  } else {
    return _lookup_1d->Eval( main.bias ) / _lookup_1d->Eval( ref.bias );
  }
}
