#include "SiPMCalib/SiPMCalc/interface/SiPMPdf.hpp"

#include "RooAbsData.h"
#include "RooConstVar.h"
#include "RooRealVar.h"
#include "Rtypes.h"

#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TMath.h"
#include "TSpectrum.h"

#include <algorithm>
#include <iostream>

// ClassImp( SiPMPdf );

// Full model construction
SiPMPdf::SiPMPdf( const char* name, const char* title,
                  RooAbsReal& _x,
                  RooAbsReal& _ped,
                  RooAbsReal& _gain,
                  RooAbsReal& _s0,
                  RooAbsReal& _s1,
                  RooAbsReal& _mean,
                  RooAbsReal& _lambda,
                  RooAbsReal& _alpha,
                  RooAbsReal& _beta,
                  RooAbsReal& _dcfrac,
                  RooAbsReal& _epsilon
                  ) :
  RooAbsPdf( name, title ),
  x(          "x",      "obs",          this, _x ),
  ped(        "ped",    "pedestal",     this, _ped ),
  gain(       "gain",   "gain",         this, _gain ),
  s0(         "s0",     "comnoise",     this, _s0 ),
  s1(         "s1",     "pixnoise",     this, _s1 ),
  mean(       "mean",   "mean",         this, _mean ),
  lambda(     "lambda", "crosstalk",    this, _lambda ),
  alpha(      "alpha",  "alpha",        this, _alpha ),
  beta(       "beta",   "beta",         this, _beta ),
  dcfraction( "dcfrac", "darkfraction", this, _dcfrac ),
  epsilon(    "eps",    "epsilon",      this, _epsilon ),
  mdistro( ped, ped+gain, epsilon, sqrt( s0*s0 + s1*s1 ) )
{
}

// No Dark current model
SiPMPdf::SiPMPdf( const char* name, const char* title,
                  RooAbsReal& _x,
                  RooAbsReal& _ped,
                  RooAbsReal& _gain,
                  RooAbsReal& _s0,
                  RooAbsReal& _s1,
                  RooAbsReal& _mean,
                  RooAbsReal& _lambda,
                  RooAbsReal& _alpha,
                  RooAbsReal& _beta
                  ) :
  RooAbsPdf( name, title ),
  x(          "x",      "obs",          this, _x ),
  ped(        "ped",    "pedestal",     this, _ped ),
  gain(       "gain",   "gain",         this, _gain ),
  s0(         "s0",     "comnoise",     this, _s0 ),
  s1(         "s1",     "pixnoise",     this, _s1 ),
  mean(       "mean",   "mean",         this, _mean ),
  lambda(     "lambda", "crosstalk",    this, _lambda ),
  alpha(      "alpha",  "alpha",        this, _alpha ),
  beta(       "beta",   "beta",         this, _beta ),
  dcfraction( "dcfrac", "darkfraction", this, RooFit::RooConst( 0 ) ),
  epsilon(    "eps",    "epsilon",      this, RooFit::RooConst( 0.01 ) ),
  mdistro( ped, ped+gain, epsilon, sqrt( s0*s0 + s1*s1 ) )
{}

SiPMPdf::SiPMPdf( const char* name, const char* title,
                  RooAbsReal& _x,
                  RooAbsReal& _ped,
                  RooAbsReal& _gain,
                  RooAbsReal& _s0,
                  RooAbsReal& _s1,
                  RooAbsReal& _mean,
                  RooAbsReal& _lambda
                  ) :
  RooAbsPdf( name, title ),
  x(          "x",      "obs",        this, _x ),
  ped(        "ped",    "pedestal",   this, _ped ),
  gain(       "gain",   "gain",       this, _gain ),
  s0(         "s0",     "comnoise",   this, _s0 ),
  s1(         "s1",     "pixnoise",   this, _s1 ),
  mean(       "mean",   "mean",       this, _mean ),
  lambda(     "lambda", "crosstalk",  this, _lambda ),
  alpha(      "alpha",  "alpha",      this, RooFit::RooConst( 0 ) ),
  beta(       "beta",   "beta",       this, RooFit::RooConst( 1000 ) ),
  dcfraction( "dcfrac", "dcfraction", this, RooFit::RooConst( 0 ) ),
  epsilon(    "eps",    "epsilon",    this, RooFit::RooConst( 0.01 ) ),
  mdistro( ped, ped+gain, epsilon, sqrt( s0*s0 + s1*s1 ) )
{}


SiPMPdf::SiPMPdf( const SiPMPdf& other, const char* name ) :
  RooAbsPdf( other, name ),
  x( "x", this, other.x ),
  ped( "ped", this,  other.ped ),
  gain( "gain", this, other.gain ),
  s0( "s0", this, other.s0 ),
  s1( "s1", this, other.s1 ),
  mean( "mean", this, other.mean ),
  lambda( "lambda", this, other.lambda ),
  alpha( "alpha", this, other.alpha ),
  beta( "beta",  this, other.beta ),
  dcfraction( "dcfrac", this, other.dcfraction ),
  epsilon( "eps", this, other.epsilon ),
  mdistro( ped, ped+gain, epsilon, sqrt( s0*s0 + s1*s1 ) )
{}

SiPMPdf::~SiPMPdf(){}

TObject*
SiPMPdf::clone( const char* name ) const
{
  return new SiPMPdf( *this, name );
}

// ------------------------------------------------------------------------------
// Implementation of function evaluation
// ------------------------------------------------------------------------------

double
SiPMPdf::evaluate() const
{
  double prob = gen_poisson( 0 ) * gauss_k( 0 );

  for( int k = 1; k < mean + 2*TMath::Sqrt( mean ) + 15; ++k ){
    double probk = binomial_prob( k, 0 ) * gauss_k( k );

    for( int i = 1; i <= k; ++i ){
      probk += binomial_prob( k, i ) * ap_eff( k, i );
    }

    prob += gen_poisson( k ) * probk;
  }

  if( prob <= 0 ){
    prob = std::numeric_limits<double>::min();
    // Forcing non-zero to avoid fit crashing.
  }

  return prob;
}

double
SiPMPdf::gen_poisson( const int k ) const
{
  if( lambda == 0 ){
    return TMath::Poisson( k, mean );
  }

  const double y = ( mean + k * lambda );
  double prod    = 1;

  for( int i = 1; i <= k; ++i ){
    prod *= y;
    prod /= (double)( i );
  }

  prod *= mean / y;
  return prod * TMath::Exp( -y );
}

double
SiPMPdf::ap_eff( const int k, const int i ) const
{
  const double pk = ped + gain * k;
  const double sk = TMath::Sqrt( s0*s0 + k*s1*s1 );
  const double y  = x - pk;

  if( i > 1 ){
    if( y < 0 ){ return 0; }

    double ans = TMath::Exp( -y / beta ) / beta;

    for( int j = 1; j <= i-1; ++j ){
      ans *= ( y / beta );
      ans /= double(j);
    }

    return ans;
  } else {
    const double ap    = TMath::Exp( -y / beta ) / beta;
    const double smear = 0.5 * ( 1+TMath::Erf( y/TMath::Sqrt( 2*sk*sk ) ) );
    return ap * smear;
  }
}

double
SiPMPdf::gauss_k( const int k ) const
{
  const double pk = ped + gain * k;
  const double sk = TMath::Sqrt( s0*s0 + k*s1*s1 );

  // Adding dark current to all Geiger discharge peaks
  if( dcfraction == 0. || k > 0  ){
    return TMath::Gaus( x, pk, sk, true );
  } else {
    mdistro.SetParam( ped, ped+gain, epsilon, TMath::Sqrt( s0*s0+s1*s1 ) );
    return ( 1-dcfraction ) * TMath::Gaus( x, pk, sk, true )
           + dcfraction * mdistro.Evaluate( x-pk );
  }
}

double
SiPMPdf::binomial_prob( const int k, const int i ) const
{
  if( alpha > 0 ){
    return TMath::BinomialI( alpha, k, i )
           - TMath::BinomialI( alpha, k, i+1 );
  } else {
    return i == 0 ? 1 : 0;
  }
}

// ------------------------------------------------------------------------------
// Functions for running variable estimation from data set.
// ------------------------------------------------------------------------------
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"

static TGraphErrors IntXTGraphErrors(
  const std::vector<double>& val,
  const std::vector<double>& error );

static TF1 FitGain( TGraphErrors& graph );
static TF1 FitWidth( TGraphErrors& graph );
static TF1 FitPoisson( TGraphErrors& graph );

static void PlotPeakFind(   const std::string& plot,
                            TH1*, TSpectrum&, std::vector<TF1>& );
static void PlotGainFit(    const std::string& plot, TGraphErrors&, TF1& );
static void PlotWidthFit(   const std::string& plot, TGraphErrors&, TF1& );
static void PlotPoissonFit( const std::string& plot, TGraphErrors&, TF1& );

void
SiPMPdf::RunEstimate( const RooAbsData& data, const std::string& plot )
{
  TH1* hist = data.createHistogram( usr::RandomString( 6 ).c_str(),
    dynamic_cast<const RooAbsRealLValue&>( x.arg() ) );
  TSpectrum spec( 20 );// 20 peaks should be plenty

  spec.Search( hist, 1, "nobackground" );

  std::vector<TF1> fitlist;
  std::vector<double> peaklist;
  std::vector<double> peakunc;
  std::vector<double> widthlist;
  std::vector<double> widthunc;
  std::vector<double> heightlist;
  std::vector<double> heightunc;

  const int maxbin     = hist->FindBin( spec.GetPositionX()[0] );
  const double histmax = hist->GetBinContent( maxbin );

  for( int i = 0; i < spec.GetNPeaks(); ++i ){
    const double x        = spec.GetPositionX()[i];
    const int bin         = hist->FindBin( spec.GetPositionX()[i] );
    const double binwidth = hist->GetBinWidth( bin );

    // Skipping over bin's with too small a content
    if( hist->GetBinContent( bin ) < ( histmax / 20. ) ){ continue; }

    fitlist.push_back( TF1( usr::RandomString( 6 ).c_str(), "gaus",
      x - 2*binwidth, x + 2*binwidth ) );
    auto& f = fitlist.back();
    f.SetParameter( 1, x );
    f.SetParameter( 2, 2*binwidth );
    hist->Fit( &f, "QN0LR" );

    // Skipping over peaks that are very non-Gaussian
    if( std::fabs( x-f.GetParameter( 1 ) ) > 2*binwidth ){
      fitlist.pop_back();
      continue;
    }
    if( std::fabs( f.GetParameter( 2 ) ) > 4*binwidth ){
      fitlist.pop_back();
      continue;
    }

    std::cout << i << " "
              << x << " "
              << f.GetParameter( 1 ) << " "
              << f.GetParameter( 2 ) << std::endl;

    peaklist.push_back( f.GetParameter( 1 ) );
    peakunc.push_back( f.GetParError( 1 ) );
    widthlist.push_back( f.GetParameter( 2 ) );
    widthunc.push_back( f.GetParError( 2 ) );
    heightlist.push_back( f.GetParameter( 0 ) );
    heightunc.push_back( std::sqrt( f.GetParameter( 0 ) ) );
  }

  std::sort( peaklist.begin(), peaklist.end() );

  TGraphErrors gaingraph = IntXTGraphErrors( peaklist, peakunc );
  TGraphErrors s1graph   = IntXTGraphErrors( widthlist, widthunc  );
  TGraphErrors poigraph  = IntXTGraphErrors( heightlist, heightunc );

  TF1 fg = FitGain( gaingraph );
  TF1 fs = FitWidth( s1graph );
  TF1 fp = FitPoisson( poigraph );

  ped    = fg.GetParameter( 1 );
  gain   = fg.GetParameter( 0 );
  s0     = fs.GetParameter( 0 );
  s1     = fs.GetParameter( 1 );
  mean   = fp.GetParameter( 1 );
  lambda = fp.GetParameter( 2 );

  if( plot != "" ){
    PlotPeakFind( plot, hist, spec, fitlist );
    PlotGainFit( plot, gaingraph, fg );
    PlotWidthFit( plot, s1graph, fs );
    PlotPoissonFit( plot, poigraph, fp  );
  }

  delete hist;
}

static TGraphErrors
IntXTGraphErrors( const std::vector<double>& val,
                  const std::vector<double>& error )
{
  TGraphErrors ans( val.size() );

  for( unsigned i = 0; i < val.size(); ++i ){
    ans.SetPoint( i, i, val.at( i ) );
    ans.SetPointError( i, 0, error.at( i ) );
  }

  return ans;
}

static TF1
FitGain( TGraphErrors& graph )
{
  // Fitting function the peak positions
  auto lf = []( const double* xx, const double* par ){
              const double x = xx[0];
              const double a = par[0];
              const double b = par[1];
              return x*a + b;
            };

  TF1 f( usr::RandomString( 6 ).c_str(), lf, 0, graph.GetN(), 2 );
  f.SetParameter( 0, graph.GetY()[1]-graph.GetY()[0] );
  f.SetParameter( 1, graph.GetY()[0] );

  graph.Fit( &f, "QRN0 EX0" );

  return f;
}

static TF1
FitWidth( TGraphErrors& graph )
{
  auto lf = []( const double* xx, const double* par ){
              const double x  = xx[0];
              const double s0 = par[0];
              const double s1 = par[1];
              return sqrt( s0*s0 + x * s1 *s1 );
            };
  TF1 f( usr::RandomString( 6 ).c_str(), lf, 0, graph.GetN(), 2 );
  f.SetParameter( 0, graph.GetY()[0] );
  f.SetParameter( 1, 0 );

  graph.Fit( &f, "QRN0 EX0" );

  return f;
}

static TF1
FitPoisson( TGraphErrors& graph )
{
  auto lf = []( const double* xx, const double* par ){
              const double x  = xx[0];
              const double N  = par[0];
              const double mu = par[1];
              const double l  = par[2];
              const double y  = ( mu + x * l );
              double prod     = 1;

              for( int i = 1; i <= x; ++i ){
                prod *= y;
                prod /= (double)( i );
              }

              prod *= mu / y;
              return N * prod * TMath::Exp( -y );
            };

  TF1 f( usr::RandomString( 6 ).c_str(), lf, 0, graph.GetN(),  3 );
  f.SetParameter( 0, graph.Integral() );
  f.SetParameter( 1, 1 );
  f.SetParameter( 2, 0.1 );

  graph.Fit( &f, "QRN0 EX0" );

  return f;
}

// ------------------------------------------------------------------------------
// Plotting Estimation results
// ------------------------------------------------------------------------------
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

void
PlotPeakFind( const std::string& plot,
              TH1* h, TSpectrum& spec, std::vector<TF1>& fitlist )
{
  usr::plt::Simple1DCanvas c;

  TH1D hist( usr::RandomString( 6 ).c_str(), ""
           , h->GetNbinsX()
           , h->GetXaxis()->GetXmin()
           , h->GetXaxis()->GetXmax() );

  for( int i = 1; i <= h->GetNcells(); ++i ){
    hist.SetBinContent( i, h->GetBinContent( i ) );
  }

  c.PlotHist( hist,
    usr::plt::PlotType( usr::plt::hist ),
    usr::plt::EntryText( "Readout" ),
    usr::plt::LineColor( usr::plt::col::black ) );

  for( auto& f : fitlist ){
    c.PlotFunc( f,
      usr::plt::PlotType( usr::plt::simplefunc ),
      usr::plt::LineColor( usr::plt::col::red ),
      &f == &fitlist.front() ?
      usr::plt::EntryText( "Local Gaussian fit" ) :
      RooCmdArg::none() );
  }

  for( int i = 0; i < spec.GetNPeaks(); ++i ){
    auto& line = c.Pad().DrawVLine(
      spec.GetPositionX()[i],
      usr::plt::LineColor( usr::plt::col::darkgray ),
      usr::plt::LineStyle( usr::plt::sty::lindotted ) );
    if( i == 0 ){
      c.Pad().AddLegendEntry( line, "Peak search results", "L" );
    }
  }

  c.DrawCMSLabel( "Peak finding", "Spectral Fit" );
  c.Pad().Xaxis().SetTitle( "Readout" );
  c.Pad().Yaxis().SetTitle( "Events" );
  c.SaveAsPDF( plot + "_peakfind.pdf" );
}

void
PlotGainFit( const std::string& plot, TGraphErrors& graph, TF1& fit )
{
  usr::plt::Simple1DCanvas c;

  c.PlotGraph( graph,
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::TrackY( usr::plt::TrackY::both ),
    usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
    usr::plt::MarkerSize( 0.5 ) );

  c.PlotFunc( fit,
    usr::plt::PlotType( usr::plt::simplefunc ),
    usr::plt::LineColor( usr::plt::col::red ) );

  c.Xaxis().SetTitle( "Discharge peak" );
  c.Yaxis().SetTitle( "Fitted peak position [ADC #times ns]" );
  c.DrawCMSLabel( "Spectral Fit" );
  c.Pad().WriteLine( usr::fstr( "Est. gain=%.0lf#pm%.1f"
                              , fit.GetParameter( 0 )
                              , fit.GetParError( 0 ) ) );

  c.SaveAsPDF( plot + "_gainfit.pdf" );
}

void
PlotWidthFit( const std::string& plot, TGraphErrors& graph, TF1& fit )
{
  usr::plt::Simple1DCanvas c;

  c.PlotGraph( graph,
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::TrackY( usr::plt::TrackY::both ),
    usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
    usr::plt::MarkerSize( 0.5 ) );
  c.PlotFunc( fit,
    usr::plt::PlotType( usr::plt::simplefunc ),
    usr::plt::LineColor( usr::plt::col::red ) );

  c.Xaxis().SetTitle( "Discharge peak" );
  c.Yaxis().SetTitle( "Fitted peak width [ADC #times ns]" );
  c.DrawCMSLabel( "Spectral Fit" );
  c.Pad().WriteLine( usr::fstr( "Est. s_{0}=%.0lf#pm%.1f"
                              , fit.GetParameter( 0 )
                              , fit.GetParError( 0 ) ) );
  c.Pad().WriteLine( usr::fstr( "Est. s_{1}=%.0lf#pm%.1f"
                              , fit.GetParameter( 1 )
                              , fit.GetParError( 1 ) ) );

  c.SaveAsPDF( plot + "_widthfit.pdf" );
}

void
PlotPoissonFit( const std::string& plot, TGraphErrors& graph, TF1& fit )
{
  usr::plt::Simple1DCanvas c;

  TH1D hist( usr::RandomString( 5 ).c_str(), "",
             graph.GetN(), -0.5, graph.GetN()-0.5 );

  for( int i = 0; i < graph.GetN(); ++i ){
    hist.Fill( i, fit.Eval( i ) );
  }

  c.PlotHist( hist,
    usr::plt::PlotType( usr::plt::hist ),
    usr::plt::LineColor( usr::plt::col::red ) );

  c.PlotGraph( graph,
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::TrackY( usr::plt::TrackY::both ),
    usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
    usr::plt::MarkerSize( 1 ) );


  c.Xaxis().SetTitle( "Discharge peak" );
  c.Yaxis().SetTitle( "Fitted peak height [Events]" );
  c.DrawCMSLabel( "Spectral Fit" );
  c.Pad().WriteLine( usr::fstr( "Est. mean=%.2lf#pm%.3f"
                              , fit.GetParameter( 1 )
                              , fit.GetParError( 1 ) ) );
  c.Pad().WriteLine( usr::fstr( "Est. #lambda=%.3lf#pm%.4f"
                              , fit.GetParameter( 2 )
                              , fit.GetParError( 2 ) ) );

  c.SaveAsPDF( plot + "_poissonfit.pdf" );
}
