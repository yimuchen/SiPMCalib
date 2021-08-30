#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

#include "TGraphErrors.h"
#include <fstream>
#include <iostream>
#include <sstream>

const double bias_min   = 1.83491;
const double pulse_min  = 21.064;
const double sipm_min   = 23.712;
double fit_parameter[8] = {
  122.165, 1.89808,        -2.61087
  ,        3.67959,        28.5222,        409.016
  ,        -1854.9,        9999.79};


double  Model( double, double, double, double* param );
TGraph* MakeGraph( const std::map<double, std::string>&, bool correct );

int
main( int argc, char* argv[] )
{
  std::map<double, std::string> z20;
  z20[0]  = "data/202006017/zmotion_stabilitysipm_center_z20_20200617_1400.txt";
  z20[10] = "data/202006017/zmotion_stabilitysipm_z10_20200617_1400.txt";
  z20[20] = "data/202006017/zmotion_stabilitysipm_z20_20200617_1400.txt";
  z20[30] = "data/202006017/zmotion_stabilitysipm_z30_20200617_1400.txt";
  z20[40] = "data/202006017/zmotion_stabilitysipm_z40_20200617_1400.txt";
  z20[50] = "data/202006017/zmotion_stabilitysipm_z50_20200617_1400.txt";
  z20[60] = "data/202006017/zmotion_stabilitysipm_z60_20200617_1400.txt";

  std::map<double, std::string> z50;
  z50[0]   = "data/202006017/zmotion_stabilitysipm_center_z50.txt";
  z50[-10] = "data/202006017/zmotion_stabilitysipm_z10_m.txt";
  z50[-20] = "data/202006017/zmotion_stabilitysipm_z20_m.txt";
  z50[-30] = "data/202006017/zmotion_stabilitysipm_z30_m.txt";
  z50[10]  = "data/202006017/zmotion_stabilitysipm_z10_p.txt";
  z50[20]  = "data/202006017/zmotion_stabilitysipm_z20_p.txt";
  z50[30]  = "data/202006017/zmotion_stabilitysipm_z30_p.txt";

  std::map<double, std::string> z80;
  z80[0]   = "data/202006017/zmotion_stabilitysipm_center_z80_zbase80.txt";
  z80[-10] = "data/202006017/zmotion_stabilitysipm_z10_zbase80.txt";
  z80[-20] = "data/202006017/zmotion_stabilitysipm_z20_zbase80.txt";
  z80[-30] = "data/202006017/zmotion_stabilitysipm_z30_zbase80.txt";
  z80[-40] = "data/202006017/zmotion_stabilitysipm_z40_zbase80.txt";
  z80[-50] = "data/202006017/zmotion_stabilitysipm_z50_zbase80.txt";
  z80[-60] = "data/202006017/zmotion_stabilitysipm_z60_zbase80.txt";
  z80[10]  = "data/202006017/zmotion_stabilitysipm_zu10_zbase80.txt";
  z80[20]  = "data/202006017/zmotion_stabilitysipm_zu20_zbase80.txt";

  TGraph* z20_orig = MakeGraph( z20, false );
  TGraph* z20_corr = MakeGraph( z20, true );
  TGraph* z50_orig = MakeGraph( z50, false );
  TGraph* z50_corr = MakeGraph( z50, true );
  TGraph* z80_orig = MakeGraph( z80, false );
  TGraph* z80_corr = MakeGraph( z80, true );


  {
    usr::plt::Simple1DCanvas c;

    c.PlotGraph( z20_orig,
      usr::plt::PlotType( usr::plt::scatter ),
      usr::plt::TrackY( usr::plt::tracky::both ),
      usr::plt::MarkerColor( usr::plt::col::black ),
      usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
      usr::plt::EntryText( "Uncorrected" ),
      usr::plt::MarkerSize( 0.6 ) );
    c.PlotGraph( z20_corr,
      usr::plt::PlotType( usr::plt::scatter ),
      usr::plt::TrackY( usr::plt::tracky::both ),
      usr::plt::MarkerColor( usr::plt::col::red ),
      usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
      usr::plt::EntryText( "Corrected" ),
      usr::plt::MarkerSize( 0.6 ) );

    c.Pad().SetYaxisMax( 1.05 );
    c.Pad().SetYaxisMin( 0.95 );
    c.Pad().DrawHLine( 1.00
                     , usr::plt::LineColor( usr::plt::col::gray )
                     , usr::plt::LineStyle( usr::plt::sty::lindashed ) );
    c.Pad().DrawHLine( 1.01
                     , usr::plt::LineColor( usr::plt::col::lightgray )
                     , usr::plt::LineStyle( usr::plt::sty::lindashed ) );
    c.Pad().DrawHLine( 0.99
                     , usr::plt::LineColor( usr::plt::col::lightgray )
                     , usr::plt::LineStyle( usr::plt::sty::lindashed ) );

    c.Pad().DrawLuminosity( "Base position z=20mm" );
    c.Pad().Xaxis().SetTitle( "Relative z motion [mm]" );
    c.Pad().Yaxis().SetTitle( "Relative readout" );

    c.SaveAsPDF( "Zstability_20.pdf" );
  }
  {
    usr::plt::Simple1DCanvas c;

    c.PlotGraph( z50_orig,
      usr::plt::PlotType( usr::plt::scatter ),
      usr::plt::TrackY( usr::plt::tracky::both ),
      usr::plt::MarkerColor( usr::plt::col::black ),
      usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
      usr::plt::EntryText( "Uncorrected" ),
      usr::plt::MarkerSize( 0.6 ) );
    c.PlotGraph( z50_corr,
      usr::plt::PlotType( usr::plt::scatter ),
      usr::plt::TrackY( usr::plt::tracky::both ),
      usr::plt::MarkerColor( usr::plt::col::red ),
      usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
      usr::plt::EntryText( "Corrected" ),
      usr::plt::MarkerSize( 0.6 ) );

    c.Pad().SetYaxisMax( 1.05 );
    c.Pad().SetYaxisMin( 0.95 );
    c.Pad().DrawHLine( 1.00
                     , usr::plt::LineColor( usr::plt::col::gray )
                     , usr::plt::LineStyle( usr::plt::sty::lindashed ) );
    c.Pad().DrawHLine( 1.01
                     , usr::plt::LineColor( usr::plt::col::lightgray )
                     , usr::plt::LineStyle( usr::plt::sty::lindashed ) );
    c.Pad().DrawHLine( 0.99
                     , usr::plt::LineColor( usr::plt::col::lightgray )
                     , usr::plt::LineStyle( usr::plt::sty::lindashed ) );

    c.Pad().DrawLuminosity( "Base position z=50mm" );
    c.Pad().Xaxis().SetTitle( "Relative z motion [mm]" );
    c.Pad().Yaxis().SetTitle( "Relative readout" );

    c.SaveAsPDF( "Zstability_50.pdf" );
  }
  {
    usr::plt::Simple1DCanvas c;

    c.PlotGraph( z80_orig,
      usr::plt::PlotType( usr::plt::scatter ),
      usr::plt::TrackY( usr::plt::tracky::both ),
      usr::plt::MarkerColor( usr::plt::col::black ),
      usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
      usr::plt::EntryText( "Uncorrected" ),
      usr::plt::MarkerSize( 0.6 ) );
    c.PlotGraph( z80_corr,
      usr::plt::PlotType( usr::plt::scatter ),
      usr::plt::TrackY( usr::plt::tracky::both ),
      usr::plt::MarkerColor( usr::plt::col::red ),
      usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
      usr::plt::EntryText( "Corrected" ),
      usr::plt::MarkerSize( 0.6 ) );

    c.Pad().SetYaxisMax( 1.05 );
    c.Pad().SetYaxisMin( 0.95 );
    c.Pad().DrawHLine( 1.00
                     , usr::plt::LineColor( usr::plt::col::gray )
                     , usr::plt::LineStyle( usr::plt::sty::lindashed ) );
    c.Pad().DrawHLine( 1.01
                     , usr::plt::LineColor( usr::plt::col::lightgray )
                     , usr::plt::LineStyle( usr::plt::sty::lindashed ) );
    c.Pad().DrawHLine( 0.99
                     , usr::plt::LineColor( usr::plt::col::lightgray )
                     , usr::plt::LineStyle( usr::plt::sty::lindashed ) );

    c.Pad().DrawLuminosity( "Base position z=80mm" );
    c.Pad().Xaxis().SetTitle( "Relative z motion [mm]" );
    c.Pad().Yaxis().SetTitle( "Relative readout" );

    c.SaveAsPDF( "Zstability_80.pdf" );
  }
  {
    usr::plt::Simple1DCanvas c;
    c.PlotGraph( z80_corr,
      usr::plt::PlotType( usr::plt::scatter ),
      usr::plt::TrackY( usr::plt::tracky::both ),
      usr::plt::MarkerColor( usr::plt::col::red ),
      usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
      usr::plt::EntryText( "Corrected" ),
      usr::plt::MarkerSize( 0.6 ) );
    c.Pad().SetYaxisMax( 1.05 );
    c.Pad().SetYaxisMin( 0.95 );
    c.Pad().DrawHLine( 1.00
                     , usr::plt::LineColor( usr::plt::col::gray )
                     , usr::plt::LineStyle( usr::plt::sty::lindashed ) );
    c.Pad().DrawHLine( 1.01
                     , usr::plt::LineColor( usr::plt::col::lightgray )
                     , usr::plt::LineStyle( usr::plt::sty::lindashed ) );
    c.Pad().DrawHLine( 0.99
                     , usr::plt::LineColor( usr::plt::col::lightgray )
                     , usr::plt::LineStyle( usr::plt::sty::lindashed ) );

    c.Pad().DrawLuminosity( "Base position z=80mm" );
    c.Pad().Xaxis().SetTitle( "Relative z motion [mm]" );
    c.Pad().Yaxis().SetTitle( "Relative readout" );
    c.SaveAsPDF( "Zstability_80_corr.pdf" );
  }
  {
    usr::plt::Simple1DCanvas c;
    c.PlotGraph( z50_corr,
      usr::plt::PlotType( usr::plt::scatter ),
      usr::plt::TrackY( usr::plt::tracky::both ),
      usr::plt::MarkerColor( usr::plt::col::red ),
      usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
      usr::plt::EntryText( "Corrected" ),
      usr::plt::MarkerSize( 0.6 ) );
    c.Pad().SetYaxisMax( 1.05 );
    c.Pad().SetYaxisMin( 0.95 );
    c.Pad().DrawHLine( 1.00
                     , usr::plt::LineColor( usr::plt::col::gray )
                     , usr::plt::LineStyle( usr::plt::sty::lindashed ) );
    c.Pad().DrawHLine( 1.01
                     , usr::plt::LineColor( usr::plt::col::lightgray )
                     , usr::plt::LineStyle( usr::plt::sty::lindashed ) );
    c.Pad().DrawHLine( 0.99
                     , usr::plt::LineColor( usr::plt::col::lightgray )
                     , usr::plt::LineStyle( usr::plt::sty::lindashed ) );

    c.Pad().DrawLuminosity( "Base position z=50mm" );
    c.Pad().Xaxis().SetTitle( "Relative z motion [mm]" );
    c.Pad().Yaxis().SetTitle( "Relative readout" );
    c.SaveAsPDF( "Zstability_50_corr.pdf" );
  }
  {
    usr::plt::Simple1DCanvas c;
    c.PlotGraph( z20_corr,
      usr::plt::PlotType( usr::plt::scatter ),
      usr::plt::TrackY( usr::plt::tracky::both ),
      usr::plt::MarkerColor( usr::plt::col::red ),
      usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
      usr::plt::EntryText( "Corrected" ),
      usr::plt::MarkerSize( 0.6 ) );
    c.Pad().SetYaxisMax( 1.05 );
    c.Pad().SetYaxisMin( 0.95 );
    c.Pad().DrawHLine( 1.00
                     , usr::plt::LineColor( usr::plt::col::gray )
                     , usr::plt::LineStyle( usr::plt::sty::lindashed ) );
    c.Pad().DrawHLine( 1.01
                     , usr::plt::LineColor( usr::plt::col::lightgray )
                     , usr::plt::LineStyle( usr::plt::sty::lindashed ) );
    c.Pad().DrawHLine( 0.99
                     , usr::plt::LineColor( usr::plt::col::lightgray )
                     , usr::plt::LineStyle( usr::plt::sty::lindashed ) );

    c.Pad().DrawLuminosity( "Base position z=20mm" );
    c.Pad().Xaxis().SetTitle( "Relative z motion [mm]" );
    c.Pad().Yaxis().SetTitle( "Relative readout" );
    c.SaveAsPDF( "Zstability_20_corr.pdf" );
  }
}


TGraph*
MakeGraph( const std::map<double, std::string>& input_files, bool correct )
{
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> ey;
  std::vector<double> z;


  // Setting up the reference point
  unsigned n_points        = 0;
  double reference_readout = 0;
  double reference_bias    = 0;
  double reference_pulse   = 0;
  double reference_sipm    = 0;

  std::ifstream infile( input_files.at( 0 ) );
  std::string line;

  double time, chipid,  x_val, y_val, z_val, bias, pulse, sipm, readout, unc;

  while( std::getline( infile, line ) ){
    std::istringstream linestream( line );
    linestream >> time >> chipid
    >>  x_val >> y_val >> z_val
    >> bias >> pulse  >> sipm
    >> readout >> unc;

    n_points++;
    reference_readout += readout;
    reference_bias    += bias;
    reference_pulse   += pulse;
    reference_sipm    += sipm;
  }

  reference_readout /= n_points;
  reference_bias    /= n_points;
  reference_pulse   /= n_points;
  reference_sipm    /= n_points;

  for( const auto& p : input_files  ){
    const double motion_z  = p.first;
    const std::string file = p.second;
    std::ifstream infile( file );
    std::string line;


    while( std::getline( infile, line ) ){
      std::istringstream linestream( line );
      linestream >> time >> chipid
      >>  x_val >> y_val >> z_val
      >> bias >> pulse  >> sipm
      >> readout >> unc;
      ;

      // Setting first point as reference
      x.push_back( motion_z );
      y.push_back( readout / reference_readout );
      ey.push_back( unc/std::sqrt( 2000 )/reference_readout );
      z.push_back( 0 );

      if( correct ){
        double scale = Model( reference_bias-bias_min
                            , reference_pulse-pulse_min
                            , reference_sipm-sipm_min
                            , fit_parameter );
        scale /= Model( bias-bias_min
                      , pulse- pulse_min
                      , sipm-sipm_min
                      , fit_parameter );
        y.back()  *= scale;
        ey.back() *= scale;
        x.back()  -= 0.5;
      }
    }
  }

  return new TGraphErrors( x.size(), x.data(), y.data(), z.data(), ey.data() );
}


double
Model( double x, double y, double z, double* param )
{
  return ( param[0] +
           param[1] * y +
           param[2] * z
           ) *
         (  param[3] +
            param[4]*x +
            param[5]*x*x +
            param[6]*x*x*x +
            param[7]*x*x*x*x
         );
}
