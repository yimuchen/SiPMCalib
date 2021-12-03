#include "SiPMCalib/SiPMCalc/interface/SiPMDarkPdf.hpp"
#include "TMath.h"

SiPMDarkPdf::SiPMDarkPdf( const char* name,
                          const char* title,
                          RooAbsReal& _x,
                          RooAbsReal& _ped,
                          RooAbsReal& _gain,
                          RooAbsReal& _s0,
                          RooAbsReal& _s1,
                          RooAbsReal& _dcfrac,
                          RooAbsReal& _epsilon ) :
  RooAbsPdf( name, title ),
  x        ( "x", "obs", this, _x ),
  ped      ( "ped", "pedestal", this, _ped ),
  gain     ( "gain", "gain", this, _gain ),
  s0       ( "s0", "comnoise", this, _s0 ),
  s1       ( "s1", "pixnoise", this, _s1 ),
  dcfrac   ( "dcfrac1", "dcfrac1", this, _dcfrac ),
  epsilon  ( "epsilon", "epsilon", this, _epsilon ),
  mdistro  ( ped, ped+gain, epsilon, sqrt( s0 * s0+s1 * s1 ) )
{}


SiPMDarkPdf::SiPMDarkPdf( const SiPMDarkPdf& other, const char* name ) :
  RooAbsPdf( other, name ),
  x        (       "x",        this, other.x    ),
  ped      (     "ped",      this, other.ped  ),
  gain     (    "gain",     this, other.gain ),
  s0       (      "s0",       this, other.s0   ),
  s1       (      "s1",       this, other.s1   ),
  dcfrac   (  "dcfrac1",  this, other.dcfrac ),
  epsilon  ( "epsilon",  this, other.epsilon ),
  mdistro  ( ped, gain, epsilon, sqrt( s0 * s0+s1 * s1 ) )
{}


TObject*
SiPMDarkPdf::clone( const char* name ) const
{
  return new SiPMDarkPdf( *this, name );
}


SiPMDarkPdf::~SiPMDarkPdf(){}

double
SiPMDarkPdf::evaluate() const
{
  mdistro.SetParam( ped, ped+gain, epsilon, sqrt( s0 * s0+s1 * s1 ) );
  return ( 1-dcfrac ) * TMath::Gaus( x, ped, s0, kTRUE )
         +dcfrac * mdistro.Evaluate( x );
}


// ------------------------------------------------------------------------------
// Running estimation
// ------------------------------------------------------------------------------
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

#include "TSpectrum.h"

void
SiPMDarkPdf::RunEstimate( const RooAbsData& data, const std::string& plot )
{
  TH1* h = data.createHistogram( usr::RandomString( 6 ).c_str(),
                                 dynamic_cast<const RooAbsRealLValue&>( x.arg() ) );
  TSpectrum spec( 5 );// 5 peaks should be plenty

  spec.Search( h, 1, "nobackground" );

  const double binwidth = h->GetBinWidth( 1 );

  // Fitting first peak ;
  TF1 g0( usr::RandomString( 6 ).c_str(), "gaus",
          spec.GetPositionX()[0]-2 * binwidth,
          spec.GetPositionX()[0]+2 * binwidth );
  g0.SetParameter( 0, h->Integral() );
  g0.SetParameter( 1, spec.GetPositionX()[0] );
  g0.SetParameter( 2, 2 * binwidth );
  h->Fit( &g0, "QN0LR" );

  // Fitting the second histogram
  TF1 g1( usr::RandomString( 6 ).c_str(), "gaus",
          spec.GetPositionX()[1]-2 * binwidth,
          spec.GetPositionX()[1]+2 * binwidth );
  g1.SetParameter( 1, spec.GetPositionX()[1] );
  g1.SetParameter( 2, 2 * binwidth );
  h->Fit( &g1, "QN0LR" );

  ped = g0.GetParameter( 1 );
  s0  = g0.GetParameter( 2 );

  s1     = sqrt( g1.GetParameter( 2 ) * g1.GetParameter( 2 )-s0 * s0 );
  gain   = g1.GetParameter( 1 )-g0.GetParameter( 1 );
  dcfrac = g1.GetParameter( 0 ) / g0.GetParameter( 0 );

  if( plot != "" ){
    TH1D hist( usr::RandomString( 6 ).c_str(), "",
               h->GetNbinsX(),
               h->GetXaxis()->GetXmin(),
               h->GetXaxis()->GetXmax() );

    for( int i = 1; i <= h->GetNcells(); ++i ){
      hist.SetBinContent( i, h->GetBinContent( i ) );
    }

    usr::plt::Simple1DCanvas c;
    c.PlotHist( hist,
                usr::plt::EntryText( "All data" ),
                usr::plt::LineColor( usr::plt::col::black ) );
    c.PlotFunc( g0,
                usr::plt::EntryText( "Pedestal peak fit" ),
                usr::plt::LineColor( usr::plt::col::red ) );
    c.PlotFunc( g1,
                usr::plt::EntryText( "1Geiger peak fit" ),
                usr::plt::LineColor( usr::plt::col::green ) );

    c.Pad().Xaxis().SetTitle( "Readout" );
    c.Pad().Yaxis().SetTitle( "Events" );

    c.SaveAsPDF( plot+"_peakfind.pdf" );
  }

  delete h;
}
