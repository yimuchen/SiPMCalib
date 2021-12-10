#include "SiPMCalib/Common/interface/StdFormat.hpp"
#include "SiPMCalib/InvSqCalc/interface/InvSqFunc.hpp"
#include "SiPMCalib/SiPMCalc/interface/NonLinearModel.hpp"
#include "SiPMCalib/SiPMCalc/interface/SiPMNonLinearFit.hpp"

#include "UserUtils/Common/interface/Maths.hpp"
#include "UserUtils/Common/interface/STLUtils.hpp"
#include "UserUtils/PlotUtils/interface/PlotCommon.hpp"
#include "UserUtils/PlotUtils/interface/Ratio1DCanvas.hpp"

void
SiPMNonLinearFit::InitDefaults()
{
  _lin_pmax  = 0;
  _lin_pmin  = 0;
  _power_z   = 500;
  _pixel_min = 1000;
  _pixel_max = 7000;

  _readout_units = "[A.U.]";
  _setup         = "UMD Gantry";
  _sipm_model    = "HDR2 (7500 pix.)";
  _sipm_id       = "";
}


SiPMNonLinearFit::SiPMNonLinearFit()
{
  InitDefaults();
}


SiPMNonLinearFit::SiPMNonLinearFit( const std::string& config )
{
  InitDefaults();
}


usr::po::options_description
SiPMNonLinearFit::FitArguments()
{
  usr::po::options_description desc( "Fitting arguments for non-linear fit" );
  desc.add_options()
    ( "linpmin",
    usr::po::value<double>(),
    "Minimum power level for linear power level fit" )
    ( "linpmax",
    usr::po::value<double>(),
    "Maximum power level for linear power level fit" )
    ( "powerz",
    usr::po::value<double>(),
    "z value to use for power level look up, will be ignored if a correction file is provided" )

    ( "npixels",
    usr::po::multivalue<double>(),
    "The number of pixels in the SiPM of interest, you can either enter a single number so that "
    "this number is fixed in during the fitting, or 2 numbers to indicate the range." )
  ;
  return desc;
}


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
 * @brief Updating the results according to the input arguments.
 */
void
SiPMNonLinearFit::UpdateSettings( const usr::ArgumentExtender& args )
{
  _lin_pmin = args.ArgOpt<double>( "linpmin", _lin_pmin );
  _lin_pmax = args.ArgOpt<double>( "linpmax", _lin_pmax );
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
 * @brief Reading the raw data file and the data correction file, storing
 * results into the class.
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


void SiPMNonLinearFit::PlotLinearity( const std::string& ) {}


/**
 * @brief In the case that the external correction table is not available.
 * Create a simplified lookup table using correcting data readout values.
 *
 * The look-up table here is generated using the data collection values at the
 * maximum z value, and will only be a correction value for the bias voltage.
 * (temperature drift correction will be omitted)
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
  _lookup_1d.reset( new ROOT::Math::Interpolator( bias,
                                                  readout,
                                                  ROOT::Math::Interpolation::
                                                  kLINEAR ));
}


/**
 * @brief Creating the data correction lookup table from an external correction
 * file.
 *
 * TODO: Properly implement
 */
void SiPMNonLinearFit::MakeInterpolator3D( const std::string& corr )
{}

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
                            zero.data() );

  const double zmin = usr::GetMinimum( z );
  const double zmax = usr::GetMaximum( z );
  const double lmin = usr::GetMinimum( lumi );
  const double lmax = usr::GetMaximum( lumi );

  _lin_func = TF1( usr::RandomString( 12 ).c_str(), InvSq_Z, zmin, zmax, 5 );
  _lin_func.FixParameter( 4, 1.0 ); // Fixing the redundent scaling paramter
  _lin_data.Fit( &_lin_func, "Q EX0 M E N 0" ); // Running first fit

  // Scaling scaling parameter according to uncertainty ratio
  const double scale = _lin_func.GetParError( 0 ) / _lin_func.GetParError( 2 );
  _lin_func.FixParameter( 4, 1 / scale  );
  _lin_func.SetParameter( 2, _lin_func.GetParameter( 2 ) *  scale );

  _lin_fit = *( _lin_data.Fit( &_lin_func, "Q EX0 M E N 0 S" ).Get()); // Saving fit result
}


/**
 * @brief Creating the graph for the number of expected photons (x axis) vs the
 * number of fired pixels as seen in the readout given the reference point for
 * a number of given photons.
 *
 * As the reference point might not exist perfectly, here we find the closest
 * point.
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
      return fabs( l.bias-ref_bias ) < fabs( r.z-ref_bias );
    }
  } );

  std::vector<double> n_in;
  std::vector<double> n_out;
  std::vector<double> in_unc;
  std::vector<double> out_unc;

  // Looping over all data
  for( const auto row : *_raw_data ){
    const double z_mult = GetZPosMultiplier( row, ref_row, ped );
    const double l_mult = GetLookUpMultiplier( row,  ref_row );

    n_in.push_back( ref_n * z_mult * l_mult );
    n_out.push_back( row.data.at( 0 ) / gain );
    in_unc.push_back( 0 );  // FIXME!!
    out_unc.push_back( row.data.at( 1 ) / gain );
  }

  _nl_data = TGraphErrors( n_in.size(),
                           n_in.data(),
                           n_out.data(),
                           in_unc.data(),
                           out_unc.data() );
}


/**
 * @brief Running the non linear fit. Definition of the SiPM will be determined by
 *
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
  _nl_func.SetParameter( 2, 0.0 );  // Pedestal
  _nl_func.SetParLimits( 2, 0.0, 100 );
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

  // Running the simple
  _nl_fit = *( _nl_data.Fit( &_nl_func, "Q EX0 M E N 0 S" ).Get());

  _nl_fit.Print();
}


/**
 * @brief Ploting the fit results of the zscan data used for inverse square estimation.
 */
void
SiPMNonLinearFit::PlotNonLinearity( const std::string& outfile )
{
  usr::plt::Ratio1DCanvas c;

  auto& fitg = c.PlotFunc(
    _nl_func,
    usr::plt::PlotType( usr::plt::fittedfunc ),
    usr::plt::VisualizeError( _nl_fit ),
    usr::plt::EntryText( "Fit" ),
    usr::plt::LineColor( usr::plt::col::blue, 1.0 ),
    usr::plt::FillColor( usr::plt::col::cyan, 1.0 ),
    RooFit::Precision( 0.0001 ),
    usr::plt::TrackY( usr::plt::tracky::max ));
  auto& datag = c.PlotGraph(
    _nl_data,
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::EntryText( "Data" ),
    usr::plt::TrackY( usr::plt::tracky::max ),
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
  c.TopPad().Yaxis().SetTitle( ( "Readout"+_readout_units ).c_str() );
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
 * @brief TODO: UPDATE DOCUMENTATION
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
    return 0;
  } else {
    return _lookup_1d->Eval( main.bias ) / _lookup_1d->Eval( ref.bias );
  }
}
