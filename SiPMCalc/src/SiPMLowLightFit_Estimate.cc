// ------------------------------------------------------------------------------
// Functions for running variable estimation from data set.
// ------------------------------------------------------------------------------
#include "SiPMCalib/SiPMCalc/interface/SiPMLowLightFit.hpp"

#include "UserUtils/Common/interface/STLUtils/OStreamUtils.hpp"
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"
#include "UserUtils/MathUtils/interface/RooFitExt.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

#include "TError.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TMath.h"
#include "TSpectrum.h"

void
SiPMLowLightFit::RunPDFEstimation()
{
  // Resetting the objects used for estimation
  _peakfits.clear();
  _gain_graph.reset( nullptr );
  _gain_fit.reset( nullptr );
  _width_graph.reset( nullptr );
  _width_fit.reset( nullptr );
  _height_graph.reset( nullptr );
  _height_fit.reset( nullptr );

  _est_hist.reset( usr::TH1DFromRooData( *_data, *_x ) );
  _spectrum.reset( new TSpectrum( 20 ) );// 20 peaks should be plenty
  _spectrum->Search( _est_hist.get(), 1, "nobackground" );

  for( int i = 0; i < _spectrum->GetNPeaks(); ++i ){
    const double x   = _spectrum->GetPositionX()[i];
    TF1*         fit = good_local_peak_fit( x );
    if( fit ){
      _peakfits.emplace_back( fit );
    }
  }

  std::sort( _peakfits.begin(),
             _peakfits.end(), []( auto& left, auto& right ){
    return left->GetParameter( 1 ) < right->GetParameter( 1 );
  } );

  // Early exit if peak finding failed to more than 1 peak.
  if( _peakfits.size() < 2 ){
    usr::log::PrintLog( usr::log::WARNING,
                        "Peak finder found less than two peaks! Results from the estimations " "routine will not be used. Check to see if data is behaving properly "
                        "running as expected" );
    return;
  }

  run_gain_est();
  run_width_est();
  run_height_est();
}


TF1*
SiPMLowLightFit::good_local_peak_fit( const double x )
{
  const int    bin      = _est_hist->FindBin( x );
  const double binwidth = _est_hist->GetBinWidth( bin );
  const int    maxbin   = _est_hist->GetMaximumBin();
  const double max      = _est_hist->GetBinContent( maxbin );

  TF1*fit = new TF1( usr::RandomString(
                       6 ).c_str(),
                     "gaus",
                     x-2 * binwidth,
                     x+2 * binwidth );
  fit->SetParameter( 1, x );
  fit->SetParameter( 2, 2 * binwidth );
  _est_hist->Fit( fit, "QN0LR" );

  // Checking if the peak is "good"
  bool is_good = true;

  // Ignoring peaks that are too small
  if( _est_hist->GetBinContent( bin ) < ( max / 20. ) ){ is_good = false; }

  // Skipping peaks that have maxmia outside of fit rage (not really a peak)
  if( std::fabs( x-fit->GetParameter( 1 ) ) > 2 * binwidth ){
    is_good = false;
  }

  // Skipping peaks that are too wide (typically intermidiate structure
  // misidentified as primary peaks)
  if( std::fabs( fit->GetParameter( 2 ) ) > 6 * binwidth ){ is_good = false; }

  if( is_good ){
    return fit;
  } else {
    delete fit;
    return nullptr;
  }
}


void
SiPMLowLightFit::run_gain_est()
{
  auto f = []( const double*xx, const double*par ){
             const double x = xx[0];
             const double a = par[0];
             const double b = par[1];
             return x * a+b;
           };

  const unsigned np = _peakfits.size();
  _gain_fit.reset( new TF1( usr::RandomString( 6 ).c_str(), f, 0, np, 2 ) );
  _gain_graph.reset( new TGraphErrors( np ) );

  for( unsigned i = 0; i < np; ++i ){
    const double x     = i;
    const double x_err = 0;
    const double y     = _peakfits.at( i )->GetParameter( 1 );
    const double y_err = _peakfits.at( i )->GetParError( 1 );

    _gain_graph->SetPoint( i, x, y );
    _gain_graph->SetPointError( i, x_err, y_err );
  }

  _gain_fit->SetParameter( 0, _gain_graph->GetY()[1]-_gain_graph->GetY()[0] );
  _gain_fit->SetParameter( 1, _gain_graph->GetY()[0] );
  _gain_graph->Fit( _gain_fit.get(), "QRN0 EX0" );

  if( _gain_fit ){
    if( !ignore_ped_est ){
      ped() = _gain_fit->GetParameter( 1 );
    }
    if( !ignore_gain_est ){
      gain() = _gain_fit->GetParameter( 0 );
    }
  }
}


void
SiPMLowLightFit::run_width_est()
{
  auto f = []( const double*xx, const double*par ){
             const double x  = xx[0];
             const double s0 = par[0];
             const double s1 = par[1];
             return sqrt( s0 * s0+x * s1 * s1 );
           };
  const unsigned np = _peakfits.size();
  _width_fit.reset( new TF1( usr::RandomString( 6 ).c_str(), f, 0, np, 2 ) );
  _width_graph.reset( new TGraphErrors( np ) );

  for( unsigned i = 0; i < np; ++i ){
    const double x     = i;
    const double x_err = 0;
    const double y     = _peakfits.at( i )->GetParameter( 2 );
    const double y_err = _peakfits.at( i )->GetParError( 2 );

    _width_graph->SetPoint( i, x, y );
    _width_graph->SetPointError( i, x_err, y_err );
  }

  _width_fit->SetParameter( 0, _width_graph->GetY()[0] );
  _width_fit->SetParameter( 1, 0 );
  _width_graph->Fit( _width_fit.get(), "QRN0 W EX0" );

  if( _width_fit ){
    if( !ignore_s0_est ){
      s0() = _width_fit->GetParameter( 0 );
    }
    if( !ignore_s1_est ){
      s1() = _width_fit->GetParameter( 1 );
    }
  }
}


void
SiPMLowLightFit::run_height_est()
{
  auto f = []( const double*xx, const double*par ){
             const double x    = xx[0];
             const double N    = par[0];
             const double mu   = par[1];
             const double l    = par[2];
             const double y    = ( mu+x * l );
             double       prod = 1;

             for( int i = 1; i <= x; ++i ){
               prod *= y;
               prod /= (double)( i );
             }

             prod *= mu / y;
             return N * prod * TMath::Exp( -y );
           };

  const unsigned np = _peakfits.size();
  _height_fit.reset( new TF1( usr::RandomString( 6 ).c_str(), f, 0, np, 3 ) );
  _height_graph.reset( new TGraphErrors( np ) );

  for( unsigned i = 0; i < np; ++i ){
    const double x     = i;
    const double x_err = 0;
    const double y     = _peakfits.at( i )->GetParameter( 0 );
    const double y_err = _peakfits.at( i )->GetParError( 0 );

    _height_graph->SetPoint( i, x, y );
    _height_graph->SetPointError( i, x_err, y_err );
  }

  _height_fit->SetParameter( 0, _height_graph->Integral() );
  _height_fit->SetParameter( 1, 1 );
  _height_fit->SetParameter( 2, 0.01 );

  _height_graph->Fit( _height_fit.get(), "QRN0 EX0" );

  if( _height_fit ){
    if( !ignore_mean_est ){
      mean() = _height_fit->GetParameter( 1 );
    }
    if( !ignore_lambda_est ){
      lambda() = _height_fit->GetParameter( 2 );
    }
  }
}
