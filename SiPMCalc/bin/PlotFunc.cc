#include "SiPMCalib/SiPMCalc/interface/SiPMDarkFunc.hpp"
#include "SiPMCalib/SiPMCalc/interface/SiPMPdf.hpp"
#include "UserUtils/Common/interface/Maths.hpp"
#include "UserUtils/Common/interface/Format.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

#include <boost/format.hpp>
#include <iostream>

int
main( int argc, char const* argv[] )
{
  const double width = 16;
  const double xmin  = usr::RoundDown( 100, width );
  const double xmax  = usr::RoundUp( 1500, width );
  const double nbins = ( xmax - xmin  ) / width;
  RooRealVar x( "x", "Readout (mV#times ns)", -64, 1500 );
  RooRealVar ped( "ped", "ped", 0, -100, 1500 );
  RooRealVar gain( "gain", "gain", 300, 0, 10000 );
  RooRealVar s0( "s0", "s0", 20, 10, 100000 );
  RooRealVar s1( "s1", "s1",  5, 0.001, 100 );
  RooRealVar mean( "mean", "mean", 1.2, 0.0001, 50 );
  RooRealVar lambda( "lambda", "lambda", 0.10, 0, 0.50 );
  RooRealVar alpha( "alpha", "alpha", 0.1, 0, 1 );
  RooRealVar beta( "beta", "beta", 30, 5, 1000 );
  RooRealVar dcfrac( "dcfrac", "dcfrac", 0.05, 0, 0.2 );
  RooRealVar epsilon( "eps", "eps", 0.01, 1e-6, 1e-1 );
  SiPMPdf p0( "p0", "p0",
              x, ped, gain, s0, s1, mean, lambda,
              alpha, beta,
              dcfrac, epsilon );
  x.setBin( nbins );

  const usr::fs::path output_dir( "results/pdf_debug/" );

  {
    usr::plt::Simple1DCanvas c( x );

    std::vector<int> color    = {kBlue, kGreen, kRed };
    std::vector<double> alist = {0.01, 0.05, 0.2};
    beta = 50;

    for( unsigned i = 0; i < alist.size(); ++i ){
      alpha = alist.at( i );

      auto& g = c.PlotPdf( p0,
        usr::plt::PlotType( usr::plt::simplefunc ),
        usr::plt::TrackY( usr::plt::TrackY::both ),
        RooFit::Precision( 1e-5 ),
        usr::plt::LineColor( color.at( i % color.size() ) ),
        usr::plt::EntryText( usr::fstr( "#alpha = %.3lf, #beta = %.1lf",
          alpha.getVal(), beta.getVal() ) )
        );
    }

    c.SetLogy( 1 );
    c.Pad().SetYaxisMin( 1e-12 );

    c.SaveAsPDF( output_dir/"compare_alpha.pdf" );
  }

  {
    usr::plt::Simple1DCanvas c( x );

    std::vector<int> color    = {kBlue, kGreen, kRed };
    std::vector<double> blist = {30, 50, 100 };

    for( unsigned i = 0; i < blist.size(); ++i ){
      beta = blist.at( i );

      c.PlotPdf( p0,
        usr::plt::PlotType( usr::plt::simplefunc ),
        usr::plt::TrackY( usr::plt::TrackY::both ),
        usr::plt::LineColor( color.at( i % color.size() ) ),
        usr::plt::EntryText( usr::fstr( "#alpha = %.3lf, #beta = %.1lf",
          alpha.getVal(), beta.getVal() ) )
        );
    }

    c.Pad().SetDataMin( 1e-12 );
    c.SetLogy( 1 );

    c.SaveAsPDF( output_dir/"compare_beta.pdf" );
  }


  {
    ped    = 25;
    gain   = 276;
    s0     = 15.87;
    s1     = 14.52;
    mean   = 1.33;
    lambda = 0.164;
    alpha  = 0.17;
    beta   = 81;
    const double xmin = x.getMin();
    const double xmax = x.getMax();
    const double sep  = 0.1;
    const unsigned np = ( xmax - xmin ) / sep;
    usr::plt::Simple1DCanvas c( x );

    TGraph p1ap( np );
    TGraph p2ap( np );
    TGraph tot( np );

    for( unsigned i = 0; i < np; ++i ){
      x = xmin + sep * i;

      const double p1y = p0.gen_poisson( 1 ) *
                         p0.binomial_prob( 1, 1 ) * p0.ap_eff( 1, 1 );
      const double p2y = p0.gen_poisson( 2 ) * (
        p0.binomial_prob( 2, 1 ) * p0.ap_eff( 2, 1 )
        +p0.binomial_prob( 2, 2 ) * p0.ap_eff( 2, 2 ) );
      const double ty = p0.Eval();

      p1ap.SetPoint( i, x.getVal(), p1y );
      p2ap.SetPoint( i, x.getVal(), p2y );
      tot.SetPoint( i, x.getVal(), ty );

    }

    std::cout << p1ap.Integral() << " " << p2ap.Integral() << std::endl;

    c.PlotGraph( p1ap, usr::plt::EntryText( "1 Geiger AP" ) );
    c.PlotGraph( p2ap, usr::plt::EntryText( "2 Geiger AP" ) );
    c.PlotGraph( tot,  usr::plt::EntryText( "Total distribution" ) );

    p1ap.SetLineColor( kBlue );
    p2ap.SetLineColor( kRed  );
    tot.SetLineColor( kBlack );

    c.SetLogy( true );
    c.Pad().SetYaxisMin( 1e-8 );

    c.SaveAsPDF( output_dir/"compare_afterpulse.pdf" );
  }

  {
    ped    = 25;
    gain   = 276;
    s0     = 15.87;
    s1     = 14.52;
    mean   = 1.33;
    lambda = 0.164;
    alpha  = 0.17;
    beta   = 81;
    const double xmin = x.getMin();
    const double xmax = x.getMax();
    const double sep  = 0.1;
    const unsigned np = ( xmax - xmin ) / sep;
    usr::plt::Simple1DCanvas c( x );

    TGraph p0m( np );
    TGraph p1m( np );
    TGraph p2m( np );
    TGraph tot( np );

    for( unsigned i = 0; i < np; ++i ){
      x = xmin + sep * i;

      const double p0y = p0.gen_poisson( 0 ) * p0.gauss_k( 0 );
      const double p1y = p0.gen_poisson( 1 ) * p0.gauss_k( 1 );
      const double p2y = p0.gen_poisson( 2 ) * p0.gauss_k( 2 );
      const double ty  = p0.Eval();

      p0m.SetPoint( i, x.getVal(), p0y );
      p1m.SetPoint( i, x.getVal(), p1y );
      p2m.SetPoint( i, x.getVal(), p2y );
      tot.SetPoint( i, x.getVal(), ty );
    }

    std::cout << p1m.Integral() << " " << p2m.Integral() << std::endl;

    c.PlotGraph( p0m, usr::plt::EntryText( "Pedestal" ) );
    c.PlotGraph( p1m, usr::plt::EntryText( "1 Geiger Discharge" ) );
    c.PlotGraph( p2m, usr::plt::EntryText( "2 Geiger Discharge" ) );
    c.PlotGraph( tot, usr::plt::EntryText( "Total distribution" ) );

    p0m.SetLineColor( kGreen );
    p1m.SetLineColor( kBlue );
    p2m.SetLineColor( kRed  );
    tot.SetLineColor( kBlack );

    c.SetLogy( true );
    c.Pad().SetYaxisMin( 1e-8 );

    c.SaveAsPDF( output_dir/"compare_geiger.pdf" );
  }

  {
    const double loedge    = 0;
    const double hiedge    = 1;
    const double width     = 0.1;
    const unsigned nsample = 300;
    MDistro md( loedge, hiedge, 1e-3, width );

    std::vector<TGraph> unsmeared;
    std::vector<TGraph> smeared;
    const std::vector<double> eplist = {0.99*1e-2, 0.99*1e-3, 0.99*1e-4};
    const std::vector<int> collist   = {kBlue, kRed, kGreen };

    for( const double ep : eplist ){
      unsmeared.push_back( TGraph( nsample ) );
      smeared.push_back( TGraph( nsample ) );

      md.SetParam( loedge, hiedge, ep, width );

      for( unsigned i = 0; i < nsample; ++i ){
        const double x = -0.5 + ( 2.0 * i ) / nsample;
        smeared.back().SetPoint( i, x, md.Evaluate( x ) );
        unsmeared.back().SetPoint( i, x, md.MFuncEval( x ) );
      }
    }

    usr::plt::Simple1DCanvas uc;
    usr::plt::Simple1DCanvas c;

    for( unsigned i = 0; i < eplist.size(); ++i ){
      c.PlotGraph( smeared.at( i ),
        usr::plt::EntryText( (
            boost::format( "#epsilon=10^{%d}" )
            % usr::GetExponent( eplist.at( i ) ) ).str() ),
        usr::plt::TrackY( usr::plt::TrackY::both )
        );
      uc.PlotGraph( unsmeared.at( i ),
        usr::plt::EntryText( (
            boost::format( "#epsilon=10^{%d}" )
            % usr::GetExponent( eplist.at( i ) ) ).str() ),
        usr::plt::TrackY( usr::plt::TrackY::both )
        );

      unsmeared.at( i ).SetLineColor( collist.at( i ) );
      smeared.at( i ).SetLineColor( collist.at( i ) );
    }

    std::cout << unsmeared.at( 0 ).Integral() << " "
              << smeared.at( 0 ).Integral() << std::endl;
    std::cout << unsmeared.at( 1 ).Integral() << " "
              << smeared.at( 1 ).Integral() << std::endl;
    std::cout << unsmeared.at( 2 ).Integral() << " "
              << smeared.at( 2 ).Integral() << std::endl;

    c.Xaxis().SetTitle( "Geiger discharge" );
    c.Yaxis().SetTitle( "Prob." );
    uc.Xaxis().SetTitle( "Geiger discharge" );
    uc.Yaxis().SetTitle( "Prob." );


    c.SaveAsPDF( output_dir/"compare_epsilon.pdf" );
    uc.SaveAsPDF( output_dir/"compare_epsilon_unsmeared.pdf" );
  }

  return 0;
}
