#include "UserUtils/Common/interface/Maths.hpp"
#include "UserUtils/Common/interface/STLUtils.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

#include "TGraphErrors.h"
#include "TMinuit.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

std::vector<double> bias_data;
std::vector<double> pulse_data;
std::vector<double> sipm_data;
std::vector<double> readout_data;
std::vector<double> uncertainty_data;

double bias_min;
double bias_max;
double pulse_min;
double pulse_max;
double sipm_min;
double sipm_max;

double model( double bias, double pulse, double sipm, double*param );
void   model_fcn( int& nparam, double*gin, double & f, double*param, int flag );

TGraph* MakeBiasProfile( const double pulse_center, const double sipm_center );
TGraph* MakeBiasFit( const double pulse_center,
                     const double sipm_center,
                     double*      param );

TGraph* MakePulseProfile( const double bias_center, const double sipm_center );
TGraph* MakePulseFit( const double bias_center,
                      const double sipm_center,
                      double*      param );

TGraph* MakeSiPMProfile( const double bias_center, const double sipm_center );
TGraph* MakeSiPMFit( const double bias_center,
                     const double sipm_center,
                     double*      param );

int
main( int argc, char*argv[] )
{
  std::ifstream infile( argv[1] );
  std::string   line;

  double time, chipid,  x, y, z, bias, pulse, sipm, readout, unc;

  while( std::getline( infile, line ) ){
    std::istringstream linestream( line );
    linestream >> time >> chipid >>  x >> y >> z >> bias >> pulse  >> sipm >>
    readout >> unc;

    bias_data.push_back( bias / 1000 );
    pulse_data.push_back( pulse );
    sipm_data.push_back( sipm );
    readout_data.push_back( -readout );
    uncertainty_data.push_back( unc / std::sqrt( 2000 ) );
  }

  bias_min  = usr::GetMinimum( bias_data );
  bias_max  = usr::GetMaximum( bias_data );
  pulse_min = usr::GetMinimum( pulse_data );
  pulse_max = usr::GetMaximum( pulse_data );
  sipm_min  = usr::GetMinimum( sipm_data );
  sipm_max  = usr::GetMaximum( sipm_data );

  // Shifting data for easier fitting
  for( unsigned i = 0; i < bias_data.size(); ++i ){
    bias_data[i]  -= bias_min;
    pulse_data[i] -= pulse_min;
    sipm_data[i]  -= sipm_min;
  }

  TMinuit minimizer( 8 );
  minimizer.SetFCN( model_fcn );
  double arglist[10];
  int    ierflg = 0;

  // minimizer.mnexcm( "SET ERR", arglist, 1, ierflg );
  minimizer.mnparm( 0, "c0", 0.1, 0.1, -1e6, 1e6, ierflg );
  minimizer.mnparm( 1, "p1", 0.1, 0.1, -1e6, 1e6, ierflg );
  minimizer.mnparm( 2, "s1", 0.1, 0.1, -1e6, 1e6, ierflg );
  minimizer.mnparm( 3, "b0", 0.1, 0.1, -1e6, 1e6, ierflg );
  minimizer.mnparm( 4, "b1", 0.1, 0.1, -1e6, 1e6, ierflg );
  minimizer.mnparm( 5, "b2", 0.1, 0.1, -1e6, 1e6, ierflg );
  minimizer.mnparm( 6, "b3", 0.1, 0.1, -1e6, 1e6, ierflg );
  minimizer.mnparm( 7, "b4", 0.1, 0.1, -1e6, 1e6, ierflg );

  arglist[0] = 500000;
  arglist[1] = 1.;
  minimizer.mnseek();
  minimizer.mnexcm( "MIGRAD", arglist, 2, ierflg );

  // Getting all fit results
  double parameter[8];
  minimizer.GetParameter( 0, parameter[0], unc );
  minimizer.GetParameter( 1, parameter[1], unc );
  minimizer.GetParameter( 2, parameter[2], unc );
  minimizer.GetParameter( 3, parameter[3], unc );
  minimizer.GetParameter( 4, parameter[4], unc );
  minimizer.GetParameter( 5, parameter[5], unc );
  minimizer.GetParameter( 6, parameter[6], unc );
  minimizer.GetParameter( 7, parameter[7], unc );

  std::cout << bias_min << " " << pulse_min << " " << sipm_min << std::endl;
  std::cout << parameter[0] << std::endl;
  std::cout << parameter[1] << std::endl;
  std::cout << parameter[2] << std::endl;
  std::cout << parameter[3] << std::endl;
  std::cout << parameter[4] << std::endl;
  std::cout << parameter[5] << std::endl;
  std::cout << parameter[6] << std::endl;
  std::cout << parameter[7] << std::endl;

  Double_t amin, edm, errdef;
  Int_t    nvpar, nparx, icstat;
  minimizer.mnstat( amin, edm, errdef, nvpar, nparx, icstat );

  // gMinuit->mnprin(3,amin);

  // Shifting back the data for plotting
  for( unsigned i = 0; i < bias_data.size(); ++i ){
    bias_data[i]  += bias_min;
    pulse_data[i] += pulse_min;
    sipm_data[i]  += sipm_min;
  }

  {
    TGraph*bg_p1 = MakeBiasProfile( 21.25, 24 );
    TGraph*bg_p2 = MakeBiasProfile( 23.5, 24 );
    TGraph*bg_p3 = MakeBiasProfile( 26.125, 24 );
    TGraph*m_p1  = MakeBiasFit( 21.25, 24, parameter );
    TGraph*m_p2  = MakeBiasFit( 23.5, 24, parameter );
    TGraph*m_p3  = MakeBiasFit( 26.125, 24, parameter );

    usr::plt::Simple1DCanvas c;

    c.PlotGraph( m_p1,
                 usr::plt::TrackY(
                   usr::plt::tracky::max ),
                 usr::plt::PlotType(
                   usr::plt::simplefunc ),
                 usr::plt::LineColor( usr::plt::col::cyan ) );

    c.PlotGraph( bg_p1,
                 usr::plt::TrackY(
                   usr::plt::tracky::max ),
                 usr::plt::PlotType(
                   usr::plt::scatter ),
                 usr::plt::EntryText(
                   "T_{LED}=21.25#circ C, T_{SiPM}=24#circ C" ),
                 usr::plt::MarkerColor( usr::plt::col::blue ),
                 usr::plt::LineColor( usr::plt::col::blue ),
                 usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
                 usr::plt::MarkerSize(
                   0.2 ) );
    c.PlotGraph( m_p2,
                 usr::plt::TrackY(
                   usr::plt::tracky::max ),
                 usr::plt::PlotType(
                   usr::plt::simplefunc ),
                 usr::plt::LineColor( usr::plt::col::pink ) );
    c.PlotGraph( bg_p2,
                 usr::plt::TrackY(
                   usr::plt::tracky::max ),
                 usr::plt::PlotType(
                   usr::plt::scatter ),
                 usr::plt::EntryText(
                   "T_{LED}=23.5#circ C, T_{SiPM}=24#circ C" ),
                 usr::plt::MarkerColor( usr::plt::col::red ),
                 usr::plt::LineColor( usr::plt::col::red ),
                 usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
                 usr::plt::MarkerSize(
                   0.2 ) );
    c.PlotGraph( m_p3,
                 usr::plt::TrackY(
                   usr::plt::tracky::max ),
                 usr::plt::PlotType(
                   usr::plt::simplefunc ),
                 usr::plt::LineColor( usr::plt::col::limegreen ) );
    c.PlotGraph( bg_p3,
                 usr::plt::TrackY(
                   usr::plt::tracky::max ),
                 usr::plt::PlotType(
                   usr::plt::scatter ),
                 usr::plt::EntryText(
                   "T_{LED}=26.13#circ C, T_{SiPM}=24#circ C" ),
                 usr::plt::MarkerColor( usr::plt::col::green ),
                 usr::plt::LineColor( usr::plt::col::green ),
                 usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
                 usr::plt::MarkerSize(
                   0.2 ) );
    c.Pad().Xaxis().SetTitle( "Bias voltage [V]" );
    c.Pad().Yaxis().SetTitle( "Readout [V-nS]" );

    c.SaveAsPDF( "Readout_v_Bias.pdf" );
  }
  {
    TGraph*bg_p1 = MakePulseProfile( 1.85, 24 );
    TGraph*bg_p2 = MakePulseProfile( 1.88, 24 );
    TGraph*bg_p3 = MakePulseProfile( 1.96, 24 );
    TGraph*m_p1  = MakePulseFit( 1.85, 24, parameter );
    TGraph*m_p2  = MakePulseFit( 1.88, 24, parameter );
    TGraph*m_p3  = MakePulseFit( 1.96, 24, parameter );

    usr::plt::Simple1DCanvas c;

    c.PlotGraph( m_p1,
                 usr::plt::TrackY(
                   usr::plt::tracky::max ),
                 usr::plt::PlotType(
                   usr::plt::simplefunc ),
                 usr::plt::LineColor( usr::plt::col::cyan ) );

    c.PlotGraph( bg_p1,
                 usr::plt::TrackY(
                   usr::plt::tracky::max ),
                 usr::plt::PlotType(
                   usr::plt::scatter ),
                 usr::plt::EntryText(
                   "Bias=1.85V, T_{SiPM}=24#circ C" ),
                 usr::plt::MarkerColor( usr::plt::col::blue ),
                 usr::plt::LineColor( usr::plt::col::blue ),
                 usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
                 usr::plt::MarkerSize(
                   0.2 ) );
    c.PlotGraph( m_p2,
                 usr::plt::TrackY(
                   usr::plt::tracky::max ),
                 usr::plt::PlotType(
                   usr::plt::simplefunc ),
                 usr::plt::LineColor( usr::plt::col::pink ) );
    c.PlotGraph( bg_p2,
                 usr::plt::TrackY(
                   usr::plt::tracky::max ),
                 usr::plt::PlotType(
                   usr::plt::scatter ),
                 usr::plt::EntryText(
                   "Bias=1.88V, T_{SiPM}=24#circ C" ),
                 usr::plt::MarkerColor( usr::plt::col::red ),
                 usr::plt::LineColor( usr::plt::col::red ),
                 usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
                 usr::plt::MarkerSize(
                   0.2 ) );
    c.PlotGraph( m_p3,
                 usr::plt::TrackY(
                   usr::plt::tracky::max ),
                 usr::plt::PlotType(
                   usr::plt::simplefunc ),
                 usr::plt::LineColor( usr::plt::col::limegreen ) );
    c.PlotGraph( bg_p3,
                 usr::plt::TrackY(
                   usr::plt::tracky::max ),
                 usr::plt::PlotType(
                   usr::plt::scatter ),
                 usr::plt::EntryText(
                   "Bias=1.96V, T_{SiPM}=24#circ C" ),
                 usr::plt::MarkerColor( usr::plt::col::green ),
                 usr::plt::LineColor( usr::plt::col::green ),
                 usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
                 usr::plt::MarkerSize(
                   0.2 ) );
    c.Pad().Xaxis().SetTitle( "LED Temperature [#circ C]" );
    c.Pad().Yaxis().SetTitle( "Readout [V-nS]" );

    c.SaveAsPDF( "Readout_v_Pulse.pdf" );
  }
  {
    TGraph*bg_p1 = MakeSiPMProfile( 1.85, 21.5 );
    TGraph*bg_p2 = MakeSiPMProfile( 1.88, 21.5 );
    TGraph*bg_p3 = MakeSiPMProfile( 1.96, 21.5 );
    TGraph*m_p1  = MakeSiPMFit( 1.85, 21.5, parameter );
    TGraph*m_p2  = MakeSiPMFit( 1.88, 21.5, parameter );
    TGraph*m_p3  = MakeSiPMFit( 1.96, 21.5, parameter );

    usr::plt::Simple1DCanvas c;

    c.PlotGraph( m_p1,
                 usr::plt::TrackY(
                   usr::plt::tracky::max ),
                 usr::plt::PlotType(
                   usr::plt::simplefunc ),
                 usr::plt::LineColor( usr::plt::col::cyan ) );

    c.PlotGraph( bg_p1,
                 usr::plt::TrackY(
                   usr::plt::tracky::max ),
                 usr::plt::PlotType(
                   usr::plt::scatter ),
                 usr::plt::EntryText(
                   "Bias=1.85V, T_{LED}=21.5#circ C" ),
                 usr::plt::MarkerColor( usr::plt::col::blue ),
                 usr::plt::LineColor( usr::plt::col::blue ),
                 usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
                 usr::plt::MarkerSize(
                   0.2 ) );
    c.PlotGraph( m_p2,
                 usr::plt::TrackY(
                   usr::plt::tracky::max ),
                 usr::plt::PlotType(
                   usr::plt::simplefunc ),
                 usr::plt::LineColor( usr::plt::col::pink ) );
    c.PlotGraph( bg_p2,
                 usr::plt::TrackY(
                   usr::plt::tracky::max ),
                 usr::plt::PlotType(
                   usr::plt::scatter ),
                 usr::plt::EntryText(
                   "Bias=1.88V, T_{LED}=21.5#circ C" ),
                 usr::plt::MarkerColor( usr::plt::col::red ),
                 usr::plt::LineColor( usr::plt::col::red ),
                 usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
                 usr::plt::MarkerSize(
                   0.2 ) );
    c.PlotGraph( m_p3,
                 usr::plt::TrackY(
                   usr::plt::tracky::max ),
                 usr::plt::PlotType(
                   usr::plt::simplefunc ),
                 usr::plt::LineColor( usr::plt::col::limegreen ) );
    c.PlotGraph( bg_p3,
                 usr::plt::TrackY(
                   usr::plt::tracky::max ),
                 usr::plt::PlotType(
                   usr::plt::scatter ),
                 usr::plt::EntryText(
                   "Bias=1.96V, T_{LED}=21.5#circ C" ),
                 usr::plt::MarkerColor( usr::plt::col::green ),
                 usr::plt::LineColor( usr::plt::col::green ),
                 usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
                 usr::plt::MarkerSize(
                   0.2 ) );
    c.Pad().Xaxis().SetTitle( "SiPM Temperature [#circ C]" );
    c.Pad().Yaxis().SetTitle( "Readout [V-nS]" );

    c.SaveAsPDF( "Readout_v_SiPM.pdf" );
  }
  return 0;
}


double
model( double x, double y, double z, double*param )
{
  return ( param[0]+param[1] * y+param[2] * z )
         * (  param[3]+param[4] * x+param[5] * x * x+param[6] * x * x * x
              +param[7]
              * x * x * x * x );
}


void
model_fcn( int& nparam, double*gin, double & f, double*param, int flag )
{
  double chi_sq = 0;
  double delta;

  for( unsigned i = 0; i < bias_data.size(); ++i ){
    delta =
      ( readout_data[i]
        -model( bias_data[i],
                pulse_data[i],
                sipm_data[i],
                param ) ) / uncertainty_data[i];
    chi_sq += delta * delta;
  }

  f = chi_sq;
}


TGraph*
MakeBiasProfile( const double pulse_center, const double sipm_center )
{

  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> ey;
  std::vector<double> z;

  for( unsigned i = 0; i < bias_data.size(); ++i  ){
    if( pulse_data[i] < pulse_center+0.25 &&
        pulse_data[i] > pulse_center-0.25 && sipm_data[i] < sipm_center+0.25 &&
        sipm_data[i] > sipm_center-0.25 ){
      x.push_back( bias_data[i] );
      y.push_back( readout_data[i] );
      ey.push_back( uncertainty_data[i] );
      z.push_back( 0 );
    }
  }

  return new TGraphErrors( x.size(), x.data(), y.data(), z.data(), ey.data() );
}


TGraph*
MakeBiasFit( const double pulse_center, const double sipm_center, double*param )
{
  double bias_g_max = usr::RoundUp( usr::GetMaximum( bias_data ), 0.005 );
  double bias_g_min = usr::RoundDown( usr::GetMinimum( bias_data ), 0.005 );

  std::vector<double> x;
  std::vector<double> y;

  for( double bias = bias_g_min; bias < bias_g_max; bias += 0.001 ){
    x.push_back( bias );
    y.push_back( model( bias-bias_min,
                        pulse_center-pulse_min,
                        sipm_center-sipm_min,
                        param ) );
  }

  return new TGraph( x.size(), x.data(), y.data() );
}


TGraph*
MakePulseProfile( const double bias_center, const double sipm_center )
{

  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> ey;
  std::vector<double> z;

  for( unsigned i = 0; i < bias_data.size(); ++i  ){
    if( bias_data[i] < bias_center+0.005 && bias_data[i] > bias_center-0.005 &&
        sipm_data[i] < sipm_center+0.125 && sipm_data[i] > sipm_center-0.125 ){
      x.push_back( pulse_data[i] );
      y.push_back( readout_data[i] );
      ey.push_back( uncertainty_data[i] );
      z.push_back( 0 );
    }
  }

  return new TGraphErrors( x.size(), x.data(), y.data(), z.data(), ey.data() );
}


TGraph*
MakePulseFit( const double bias_center, const double sipm_center, double*param )
{
  double pulse_g_max = usr::RoundUp( usr::GetMaximum( pulse_data ), 0.01 );
  double pulse_g_min = usr::RoundDown( usr::GetMinimum( pulse_data ), 0.01 );

  std::vector<double> x;
  std::vector<double> y;

  for( double pulse = pulse_g_min; pulse < pulse_g_max; pulse += 0.01 ){
    x.push_back( pulse );
    y.push_back( model( bias_center-bias_min,
                        pulse-pulse_min,
                        sipm_center-sipm_min,
                        param ) );
  }

  return new TGraph( x.size(), x.data(), y.data() );
}


TGraph*
MakeSiPMProfile( const double bias_center, const double pulse_center )
{

  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> ey;
  std::vector<double> z;

  for( unsigned i = 0; i < bias_data.size(); ++i  ){
    if( bias_data[i] < bias_center+0.005 && bias_data[i] > bias_center-0.005 &&
        pulse_data[i] < pulse_center+0.125 &&
        pulse_data[i] > pulse_center-0.125 ){
      x.push_back( sipm_data[i] );
      y.push_back( readout_data[i] );
      ey.push_back( uncertainty_data[i] );
      z.push_back( 0 );
    }
  }

  return new TGraphErrors( x.size(), x.data(), y.data(), z.data(), ey.data() );
}


TGraph*
MakeSiPMFit( const double bias_center, const double pulse_center, double*param )
{
  double sipm_g_max = usr::RoundUp( usr::GetMaximum( sipm_data ), 0.01 );
  double sipm_g_min = usr::RoundDown( usr::GetMinimum( sipm_data ), 0.01 );

  std::vector<double> x;
  std::vector<double> y;

  for( double sipm = sipm_g_min; sipm < sipm_g_max; sipm += 0.01 ){
    x.push_back( sipm );
    y.push_back( model( bias_center-bias_min,
                        pulse_center-pulse_min,
                        sipm-sipm_min,
                        param ) );
  }

  return new TGraph( x.size(), x.data(), y.data() );
}
