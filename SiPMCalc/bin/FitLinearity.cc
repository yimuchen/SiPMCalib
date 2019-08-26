#include "SiPMCalib/SiPMCalc/interface/NonLinearModel.hpp"

#include "UserUtils/Common/interface/STLUtils/OStreamUtils.hpp"
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"
#include "UserUtils/MathUtils/interface/Measurement/Measurement.hpp"
#include "UserUtils/PlotUtils/interface/Flat2DCanvas.hpp"
#include "UserUtils/PlotUtils/interface/Ratio1DCanvas.hpp"

#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TStyle.h"

#include <boost/format.hpp>
#include <cmath>
#include <iostream>

// Global parameters
const unsigned N = 1584;
// const double gain = 156;

usr::Measurement A( const unsigned idx );
usr::Measurement B( const unsigned idx );

double LinearModel( const double*, const double* );
double LOModel( const double*, const double* );
double NLOModel( const double*, const double* );

int
main( int argc, char* argv[] )
{
  const double en0 = 5.5 / ( A( 1 )*B( 6 ) ).CentralValue();

  for( int i = 1; i <= 6; ++i ){
    usr::fout( "%d | %.5lf %.6lf %.3lf | %.5lf %.6lf %.3lf\n"
             , i
             ,   A( i ).CentralValue()
             ,   A( i ).AbsAvgError()
             ,   fabs( std::log10( A( i ).CentralValue() ) )
             , B( i ).CentralValue()
             , B( i ).AbsAvgError()
             , fabs( std::log10( B( i ).CentralValue() ) )
      );
  }

  struct LumiPoint
  {
    unsigned Aidx;
    unsigned Bidx;
    double   lumival;
    double   lumiunc;
  };

  std::vector<LumiPoint> raw;

  raw.push_back( {1, 1,  21176.583, 404.0374} );
  raw.push_back( {1, 2,  12220.202, 494.3889} );
  raw.push_back( {1, 3,  5899.950, 362.7447} );
  raw.push_back( {1, 4,  2022.905, 250.4348} );
  raw.push_back( {1, 5,  352.929, 106.1599} );
  raw.push_back( {1, 6,  99.291, 50.8222} );

  raw.push_back( {2, 1,  23572.173, 439.5601} );
  raw.push_back( {2, 2,  14964.365, 412.0734} );
  raw.push_back( {2, 3,  7681.354, 409.9390} );
  raw.push_back( {2, 4,  2918.226, 290.4511} );
  raw.push_back( {2, 5,  504.682, 127.1172} );
  raw.push_back( {2, 6,  114.440, 57.0368} );

  raw.push_back( {3, 1,  24788.129, 481.3429} );
  raw.push_back( {3, 2,  16191.926, 516.0712} );
  raw.push_back( {3, 3,  8618.587, 405.7989} );
  raw.push_back( {3, 4,  3416.504, 310.4952} );
  raw.push_back( {3, 5,  595.080, 139.5651} );
  raw.push_back( {3, 6,  135.586, 63.1189} );

  raw.push_back( {4, 1,  27282.977, 397.2383} );
  raw.push_back( {4, 2,  18394.781, 489.2027} );
  raw.push_back( {4, 3,  10467.488, 409.0632} );
  raw.push_back( {4, 4,  4553.770, 341.6897} );
  raw.push_back( {4, 5,  818.038, 163.5040} );
  raw.push_back( {4, 6,  173.014, 72.2192} );

  raw.push_back( {5, 1,  28511.370, 551.3734} );
  raw.push_back( {5, 2,  20079.564, 440.2232} );
  raw.push_back( {5, 3,  12081.811, 418.8330} );
  raw.push_back( {5, 4,  5423.645, 365.2495} );
  raw.push_back( {5, 5,  1011.776, 183.0310} );
  raw.push_back( {5, 6,  191.629, 73.1937} );

  raw.push_back( {6, 1,  30130.076, 465.7346} );
  raw.push_back( {6, 2,  22153.314, 418.2185} );
  raw.push_back( {6, 3,  13517.337, 453.9808} );
  raw.push_back( {6, 4,  6564.868, 383.1526} );
  raw.push_back( {6, 5,  1205.762, 197.0262} );
  raw.push_back( {6, 6,  236.283, 85.3561} );

  std::sort( raw.begin(), raw.end(),
    []( const LumiPoint& l, const LumiPoint& r ) -> bool {
    return ( A( l.Aidx ).CentralValue()*B( l.Bidx ).CentralValue() )
    < ( A( r.Aidx ).CentralValue()*B( r.Bidx ).CentralValue() );
  } );

  std::vector<TGraphErrors> glist;
  glist.push_back( TGraphErrors( raw.size() ) );

  for( unsigned i = 0; i < 6; i++ ){
    glist.push_back( TGraphErrors( 6 ) );
  }

  for( unsigned i = 0; i < raw.size(); ++i ){
    const unsigned idx1 = raw[i].Aidx;
    const unsigned idx2 = raw[i].Bidx;
    const double xval   = ( A( idx1 ).CentralValue()*B( idx2 ).CentralValue() );
    const double xunc   = ( A( idx1 )*B( idx2 ) ).AbsAvgError();
    // Scaling of luminosity is required to keep fit in a reasonable value
    const double lumi = raw[i].lumival / 1e3;
    const double lunc = raw[i].lumiunc / 1e3;
    glist[0].SetPoint( i, en0*xval /N, lumi );
    glist[0].SetPointError( i, en0*xunc/ N, lunc );

    glist[ idx2 ].SetPoint( idx1-1, en0*xval /N, lumi );
    glist[ idx2 ].SetPointError( idx1-1, en0*xunc /N, lunc );
  }

  const int colorlist[7] = {
    kBlue,
    kBlack,
    kGreen + 2,
    kOrange -1,
    kRed +2,
    kPink + 7,
    kViolet + 1
  };

  const double xmin = glist[0].GetX()[0];
  const double xmax = glist[0].GetX()[ raw.size() -1 ];

  // Running the fit for a linear model
  TF1 f( "f", LinearModel, xmin, xmax, 2 );
  f.SetRange( xmin, 0.3  );
  auto fit = glist[0].Fit( &f, "QN0 EX0 R S" );

  std::cout << "==== LINEAR FIT ====================" << std::endl;
  std::cout << "PDE x Gain: "
            << f.GetParameter( 0 ) << "+-"
            << f.GetParError( 0 )  << std::endl;

  // Fitting leading order model
  TF1 expf( "ExpF", LOModel,  xmin, 1.3, 3 );
  expf.SetParameter( 0, 1 );
  expf.SetParameter( 1, 1 );
  expf.SetParameter( 2, 0 );
  glist[0].Fit( &expf, "QN0 EX0 R" );
  glist[0].Fit( &expf, "QN0 EX0 R" );
  glist[0].Fit( &expf, "QN0 EX0 R" );
  auto fitexp = glist[0].Fit( &expf, "QN0 EX0 R S" );

  std::cout << "==== LEADING ORDER ====================" << std::endl;
  std::cout << "PDE :"
            << expf.GetParameter( 1 )  << " +- " << expf.GetParError( 1 )
            << std::endl;
  std::cout << "Gain: "
            << expf.GetParameter( 0 ) << " +- " << expf.GetParError( 0 )
            << std::endl;

  TF1 nlof( "NLOF", NLOModel, xmin, xmax, 5 );
  nlof.SetParameter( 0, expf.GetParameter( 0 ) );
  nlof.SetParameter( 1, expf.GetParameter( 1 ) );
  nlof.SetParameter( 2, 0 );
  nlof.SetParameter( 3, 17 );
  nlof.SetParameter( 4, expf.GetParameter( 2 ) );
  glist[0].Fit( &nlof, "QN0 EX0 R" );
  glist[0].Fit( &nlof, "QN0 EX0 R" );
  glist[0].Fit( &nlof, "QN0 EX0 R" );
  auto fitnlo = glist[0].Fit( &nlof, "QN0 EX0 R S" );

  std::cout << "==== NEXT TO LEADING ORDER ====================" << std::endl;
  std::cout << "PDE :"
            << nlof.GetParameter( 1 )  << " +- " << nlof.GetParError( 1 )
            << std::endl;
  std::cout << "Gain: "
            << nlof.GetParameter( 0 ) << " +- " << nlof.GetParError( 0 )
            << std::endl;
  std::cout << "alpha :"
            << nlof.GetParameter( 2 )  << " +- " << nlof.GetParError( 2 )
            << std::endl;
  std::cout << "beta: "
            << nlof.GetParameter( 3 ) << " +- " << nlof.GetParError( 3 )
            << std::endl;

  usr::plt::Ratio1DCanvas c;

  c.PlotGraph( glist[0],
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::TrackY( usr::plt::TrackY::both ) );
  glist[0].SetMarkerColor( kGray+2 );
  glist[0].SetLineColor( kGray );
  glist[0].SetMarkerStyle( usr::plt::sty::mkrcircle );
  glist[0].SetMarkerSize( 0.2 );

  for( int i = 1; i < 7; ++i ){
    auto& g = c.PlotGraph( glist[i],
      usr::plt::PlotType( usr::plt::scatter ),
      usr::plt::TrackY( usr::plt::TrackY::both ) );
    g.SetMarkerStyle( usr::plt::sty::mkrcircle );
    g.SetMarkerSize( 0.2 );
    g.SetMarkerColor( colorlist[i] );
    g.SetLineColor( colorlist[i] );
  }

  auto& fg = c.PlotFunc( f,
    usr::plt::PlotType( usr::plt::fittedfunc ),
    usr::plt::VisualizeError( fit ),
    usr::plt::EntryText( usr::fstr( "Linear_{GOF=%.1lf/%d}"
                                  , fit->Chi2(), fit->Ndf() ) ),
    RooFit::Precision( 1e-4 )
    );
  fg.SetLineColor( kBlue );
  fg.SetFillColorAlpha( kCyan, 0.3 );

  auto& fexpg = c.PlotFunc( expf,
    usr::plt::PlotType( usr::plt::fittedfunc ),
    usr::plt::VisualizeError( fitexp ),
    usr::plt::EntryText( usr::fstr( "LO_{GOF=%.1lf/%d}"
                                  , fitexp->Chi2(), fitexp->Ndf() ) ),
    RooFit::Precision( 1e-4 )
    );
  fexpg.SetLineColor( usr::plt::col::red );
  fexpg.SetFillColorAlpha( usr::plt::col::pink, 0.3 );


  auto& fnlog = c.PlotFunc( nlof,
    usr::plt::PlotType( usr::plt::fittedfunc ),
    usr::plt::VisualizeError( fitnlo ),
    usr::plt::EntryText( usr::fstr( "NLO_{GOF=%.1lf/%d}"
                                  , fitnlo->Chi2(), fitnlo->Ndf() ) ),
    RooFit::Precision( 1e-4 )
    );
  fnlog.SetLineColor( usr::plt::col::green );
  fnlog.SetFillColorAlpha( usr::plt::col::green, 0.3 );

  c.PlotScale( fnlog, fnlog, usr::plt::PlotType( usr::plt::fittedfunc ) );
  c.PlotScale( fg,    fnlog, usr::plt::PlotType( usr::plt::simplefunc ) );
  c.PlotScale( fexpg, fnlog, usr::plt::PlotType( usr::plt::simplefunc ) );

  c.BottomPad().DrawVLine(
    f.GetXmax(),
    usr::plt::LineColor( usr::plt::col::darkgray ),
    usr::plt::LineStyle( usr::plt::sty::linshortdash ) );
  c.BottomPad().DrawVLine(
    expf.GetXmax(),
    usr::plt::LineColor( usr::plt::col::darkgray ),
    usr::plt::LineStyle( usr::plt::sty::lindensedot ) );
  c.BottomPad().DrawVLine( 1,
    usr::plt::LineColor( usr::plt::col::darkgray ),
    usr::plt::LineStyle( usr::plt::sty::linsolid ) );


  for( int i = 1; i < 7; ++i ){
    c.PlotScale( glist[i], fnlog,
      usr::plt::PlotType( usr::plt::scatter ) );
  }


  c.DrawLuminosity( "Laser Setup" );
  c.DrawCMSLabel( "", "Linearity Test" );
  c.TopPad()
  .WriteLine( "Hamamatsu S13360-2050VE" )
  .WriteLine( "N_{pixel} = 1584" )
  .WriteLine( "V_{bias} = 51.5 V" )
  .WriteLine( usr::fstr( "P_{sec} = %.1lf#pm%.1lf%%"
                       , nlof.GetParameter( 2 )*100
                       , nlof.GetParError( 2 )*100 ) )
  .WriteLine( usr::fstr( "#tau_{SiPM}/#tau_{pulse} = %.1lf#pm%.1lf"
                       , nlof.GetParameter( 3 )
                       , nlof.GetParError( 3 ) ) )
  ;


  for( int i = 1; i < 7; ++i ){
    const double nd   = TMath::Log( B( i ).CentralValue() ) / TMath::Log( 10 );
    const double xmin = glist[i].GetX()[0];
    const double xmax = glist[i].GetX()[glist[i].GetN()-1];
    const double xcen = sqrt( xmin * xmax );
    c.BottomPad().SetTextAlign( usr::plt::font::bottom_center );
    c.BottomPad().WriteAtData(
      xcen, 1.5,
      ( boost::format( "ND = %.2lf" )% fabs( nd ) ).str(),
      usr::plt::TextColor( colorlist[i] ),
      usr::plt::TextSize( c.Font().tiny() ) );
  }

  c.TopPad().SetDataMin( glist[0].GetY()[0] * 0.3 );
  c.TopPad().SetLogx( 1 );
  c.BottomPad().SetLogx( 1 );
  c.SetLogy( 1 );
  c.BottomPad().Xaxis().SetTitle( "#bar{N}_{photon} / N_{pix}" );
  c.TopPad().Yaxis().SetTitle( "Readout [V-#mu s]" );
  c.BottomPad().Yaxis().SetTitle( "Data/NLO Fit" );
  c.TopPad().FinalizeLegend( usr::plt::align::bottom_right );
  c.SaveAsPDF( "Lineartest.pdf" );


  {
    const auto corrmatrix = fitnlo->GetCorrelationMatrix();
    const unsigned size   = corrmatrix.GetNrows();

    TH2D corrhist( "corrhist", "", size, 0, size,
                   size, 0, size  );

    for( unsigned i = 0; i < size; ++i ){
      for( unsigned j = 0; j < size; ++j ){
        corrhist.Fill( i, j, corrmatrix( i, j ) );
      }
    }

    usr::plt::Flat2DCanvas c;
    c.PlotHist( corrhist, usr::plt::Plot2DF( usr::plt::heattext ) );
    corrhist.SetMaximum( 1 );
    corrhist.SetMinimum( -1 );

    gStyle->SetPalette( kBlackBody );
    const char binlabel[5][10] = {
      "Gain", "PDE", "P_{sec}", "#beta", "Ped"
    };

    for( unsigned i = 0; i < size; ++i ){
      c.Xaxis().SetBinLabel( i+1, binlabel[i] );
      c.Yaxis().SetBinLabel( i+1, binlabel[i] );
    }

    c.Zaxis().SetTitle( "Correlation" );
    gStyle->SetPaintTextFormat( ".2lf" );

    c.SaveAsPDF( "LinearFit_Correlation.pdf" );
  }
  return 0;
}


usr::Measurement
A( const unsigned idx )
{
  static const double Araw_lo[7] = {
    0.7081,
    1668.4,
    2508.0,
    2977.1,
    4195.3,
    5273.5,
    6578.1
  };

  static const double Araw_hi[7] = {
    0.7015,
    1670.6,
    2511.8,
    2983.0,
    4200.5,
    5286.7,
    6604.4
  };
  if( idx < 1 || idx > 6 ){ return usr::Measurement( 0, 0, 0 ); }
  if( idx == 6 ){ return usr::Measurement( 1, 0, 0 ); }
  const double lo   = ( Araw_lo[idx] - Araw_hi[0] )/( Araw_hi[6] - Araw_lo[0] );
  const double hi   = ( Araw_hi[idx] - Araw_lo[0] )/( Araw_lo[6] - Araw_hi[0] );
  const double diff = fabs( hi-lo );

  return usr::Measurement( ( lo+hi )/2, diff, diff );
}

usr::Measurement
B( const unsigned idx )
{
  static const double Braw_lo[7] = {
    0.696,
    31069,
    8633.1,
    2480.4,
    701.10,
    97.002,
    20.379
  };

  static const double Braw_hi[7] = {
    0.700,
    31328,
    8691.0,
    2497.4,
    704.95,
    97.330,
    20.470
  };
  if( idx < 1 || idx > 6 ){ return usr::Measurement( 0, 0, 0 ); }
  if( idx == 1 ){ return usr::Measurement( 1, 0, 0 ); }
  const double lo   = ( Braw_lo[idx] - Braw_hi[0] )/( Braw_hi[1] - Braw_lo[0] );
  const double hi   = ( Braw_hi[idx] - Braw_lo[0] )/( Braw_lo[1] - Braw_hi[0] );
  const double diff = fabs( hi-lo )/2;

  return usr::Measurement( ( lo+hi )/2, diff, diff );
}
