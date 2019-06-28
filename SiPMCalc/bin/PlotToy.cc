#include "SiPMCalib/SiPMCalc/interface/ToyRunCommon.hpp"

#include "UserUtils/MathUtils/interface/RooFitExt.hpp"
#include "UserUtils/PlotUtils/interface/Ratio1DCanvas.hpp"

#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooRealVar.h"

#include <boost/format.hpp>
#include <iostream>
#include <sstream>


int
main( int argc, char** argv )
{
  usr::ArgumentExtender arg( "data/modelcfg.json" );
  arg.AddOptions( ToyOptions() );
  arg.ParseOptions( argc, argv );

  // Setting up variables for reading.
  double status;
  double ped, pederr;
  double gain, gainerr;
  double s0, s0err;
  double s1, s1err;
  double mean, meanerr;
  double lambda, lambdaerr;
  double dcfrac, dcfracerr;
  double alpha, alphaerr;
  double beta, betaerr;
  double fittime;

  // Setting up variables for fitting;
  RooRealVar pull( "pull", "Pull value", -5, 5 );
  pull.setRange( "full",   -5, 5 );
  pull.setRange( "narrow", -2, 2 );

  RooDataSet pedpull(   "ped",    "Pedestal",        RooArgSet( pull ) );
  RooDataSet gainpull(  "gain",   "Gain",            RooArgSet( pull ) );
  RooDataSet s0pull(    "s0",     "Common noise",    RooArgSet( pull ) );
  RooDataSet s1pull(    "s1",     "Pixel noise",     RooArgSet( pull ) );
  RooDataSet meanpull(  "mean",   "Mean",            RooArgSet( pull ) );
  RooDataSet lambdapull( "lambda", "Cross talk rate", RooArgSet( pull ) );
  RooDataSet dcfracpull( "dcfrac", "Dark count rate", RooArgSet( pull ) );
  RooDataSet alphapull( "alpha",  "alpha",           RooArgSet( pull ) );
  RooDataSet betapull(  "beta",   "beta",            RooArgSet( pull ) );

  std::ifstream fin;
  std::string line;
  unsigned linecount = 0;
  double pullval;
  fin.open( filename( arg ), std::ifstream::in );

  std::cout << filename( arg ) << std::endl;

  while( std::getline( fin, line ) ){

    std::istringstream linestream( line );
    linestream >> status
    >> ped >> pederr
    >> gain >> gainerr
    >> s0 >> s0err
    >> s1 >> s1err
    >> mean >> meanerr
    >> lambda >> lambdaerr
    >> dcfrac >> dcfracerr
    >> alpha >> alphaerr
    >> beta >> betaerr
    >> fittime;

    std::cout << "\rLine: " << linecount  << "  "
              << ped << " " << pederr
              << std::flush;
    ++linecount;

    // Skipping if excessively large pull value is obtained.
    pullval = ( ped - Pedestal( arg ) )/ pederr;
    if( pederr > 1e-6 ){
      if( pullval < pull.getMin()  || pullval > pull.getMax() ){continue;}
    }

    pullval = ( gain - Gain( arg ) )/ gainerr;
    if( gainerr > 1e-6 ){
      if( pullval < pull.getMin()  || pullval > pull.getMax() ){continue;}}

    pullval = ( s0 - ComNoise( arg ) )/ s0err;
    if( s0err > 1e-6 ){
      if( pullval < pull.getMin()  || pullval > pull.getMax() ){continue;}}

    pullval = ( s1 - PixNoise( arg ) )/ s1err;
    if( s1err > 1e-6 ){
      if( pullval < pull.getMin()  || pullval > pull.getMax() ){continue;}}

    pullval = ( mean - Mean( arg ) )/ meanerr;
    if( meanerr > 1e-6 ){
      if( pullval < pull.getMin()  || pullval > pull.getMax() ){continue;}}

    pullval = ( lambda - Lambda( arg ) )/ lambdaerr;
    if( lambdaerr > 1e-6 ){
      if( pullval < pull.getMin()  || pullval > pull.getMax() ){continue;}}

    pullval = ( dcfrac - DCFrac( arg ) )/ dcfracerr;
    if( dcfracerr > 1e-6 ){
      if( pullval < pull.getMin()  || pullval > pull.getMax() ){continue;}}

    pullval = ( alpha - Alpha( arg ) )/ alphaerr;
    if( alphaerr > 1e-6 ){
      if( pullval < pull.getMin()  || pullval > pull.getMax() ){continue;}}

    pullval = ( beta - Beta( arg ) )/ betaerr;
    if( betaerr > 1e-6 ){
      if( pullval < pull.getMin()  || pullval > pull.getMax() ){continue;}}

    // Filling in datasets
    pull = ( ped - Pedestal( arg ) )/ pederr;
    pedpull.add( RooArgSet( pull ) );

    pull = ( gain - Gain( arg ) )/ gainerr;
    gainpull.add( RooArgSet( pull ) );

    pull = ( s0 - ComNoise( arg ) )/ s0err;
    s0pull.add( RooArgSet( pull ) );

    pull = ( s1 - PixNoise( arg ) )/ s1err;
    s1pull.add( RooArgSet( pull ) );

    pull = ( mean - Mean( arg ) )/ meanerr;
    meanpull.add( RooArgSet( pull ) );

    pull = ( lambda - Lambda( arg ) )/ lambdaerr;
    lambdapull.add( RooArgSet( pull ) );

    pull = ( dcfrac - DCFrac( arg ) )/ dcfracerr;
    dcfracpull.add( RooArgSet( pull ) );

    pull = ( alpha - Alpha( arg ) )/ alphaerr;
    alphapull.add( RooArgSet( pull ) );

    pull = ( beta - Beta( arg ) )/ betaerr;
    betapull.add( RooArgSet( pull ) );

    std::cout << " " << pedpull.sumEntries() <<  " ";
  }

  std::cout << "...Done!" << std::endl;


  RooRealVar pullmean( "mean", "mean", 0, -2, 2 );
  RooRealVar pullsig( "sig", "sig", 1, 0.01, 4 );
  RooGaussian pdf( "pdf", "pdf", pull, pullmean, pullsig );
  RooRealVar npullmean( "nmean", "nmean", 0, -2, 2 );
  RooRealVar npullsig( "nsig", "nsig", 1, 0.01, 2 );
  RooGaussian npdf( "npdf", "npdf", pull, npullmean, npullsig );

  std::ofstream pullfile( pullfilename( arg ), std::ofstream::out );

  for( RooDataSet* data : { &pedpull,
                            &gainpull,
                            &s0pull,
                            &s1pull,
                            &meanpull,
                            &lambdapull,
                            &dcfracpull,
                            &alphapull,
                            &betapull } ){
    pullmean  = 0;
    pullsig   = 1;
    npullmean = 0;
    npullsig  = 1;

    RooDataSet* ndata = dynamic_cast<RooDataSet*>(
      data->reduce( RooFit::CutRange( "narrow" ) ) );

    auto fit = pdf.fitTo(
      *data, RooFit::Save(),
      RooFit::Warnings( false ),
      RooFit::PrintLevel( -100 ),
      RooFit::PrintEvalErrors( -100 ),
      RooFit::Verbose( false )
      );

    auto nfit = npdf.fitTo(
      *ndata, RooFit::Save(),
      RooFit::Range( "narrow" ),
      RooFit::Warnings( false ),
      RooFit::PrintLevel( -100 ),
      RooFit::PrintEvalErrors( -100 ),
      RooFit::Verbose( false )
      );

    const double ksprob  = usr::KSProb( *data, pdf, pull );
    const double nksprob = usr::KSProb( *ndata, npdf, pull );

    usr::plt::Ratio1DCanvas c(
      usr::plt::RangeByVar( pull ),
      usr::plt::Ratio1DCanvas::default_height );

    // c.PlotData( data,
    //   RooFit::Invisible() );

    auto& gpdf = c.PlotPdf( pdf,
      usr::plt::PlotType( usr::plt::fittedfunc ),
      RooFit::VisualizeError( *fit, 1, false ),
      RooFit::Normalization( data->sumEntries(), RooAbsReal::NumEvent ),
      RooFit::Range( "full" ),
      RooFit::NormRange( "full" ),
      usr::plt::EntryText( ( boost::format(
        "#mu=%.2lf_{#pm%.3lf} #sigma=%.2lf_{#pm%.3lf}\n(G.o.F p-val.=%.3lf)" )
                             % pullmean.getVal() % pullmean.getError()
                             % pullsig.getVal()  % pullsig.getError()
                             % ksprob ).str() )
      );
    auto& gnpdf = c.PlotPdf( npdf,
      usr::plt::PlotType( usr::plt::fittedfunc ),
      RooFit::VisualizeError( *nfit, 1, false ),
      RooFit::Normalization( ndata->sumEntries(), RooAbsReal::NumEvent ),
      RooFit::Range( "narrow" ),
      RooFit::NormRange( "narrow" ),
      usr::plt::EntryText( ( boost::format(
        "#mu=%.2lf_{#pm%.3lf} #sigma=%.2lf_{#pm%.3lf}\n(G.o.F p-val=%.3lf)" )
                             % npullmean.getVal() % npullmean.getError()
                             % npullsig.getVal()  % npullsig.getError()
                             % nksprob ).str() )
      );
    auto& gdata = c.PlotData( data,
      usr::plt::PlotType( usr::plt::scatter ),
      usr::plt::EntryText( "Toy MC values" )
      );


    // Styling
    gpdf.SetLineColor( kBlue );
    gpdf.SetFillStyle( usr::plt::sty::fillsolid );
    gpdf.SetFillColorAlpha( kCyan, 0.3 );

    gnpdf.SetLineColor( kRed );
    gnpdf.SetFillStyle( usr::plt::sty::fillsolid );
    gnpdf.SetFillColorAlpha( kPink, 0.3 );

    gdata.SetMarkerSize( 0.2 );


    auto& rdata = c.PlotScale( gdata, gpdf,
      usr::plt::PlotType( usr::plt::scatter ) );
    auto& rndata = c.PlotScale( gdata, gnpdf,
      usr::plt::PlotType( usr::plt::scatter ) );

    rdata.SetLineColor( kBlue );
    rdata.SetMarkerColor( kBlue );
    rndata.SetLineColor( kRed );
    rndata.SetMarkerColor( kRed );

    c.BottomPad().Yaxis().SetTitle( "Data/Fit" );
    c.TopPad().DrawCMSLabel( usr::plt::cap::sim, "HGCal" );
    c.TopPad().DrawLuminosity( ( boost::format( "Num. Events = %ld" )
                                 % arg.Arg<int>( "nEvents" ) ).str() );
    c.BottomPad().Xaxis().SetTitle(
      ( std::string( "Pull value (" )
        + data->GetTitle()
        + std::string( ")" ) ).c_str()
      );


    // Saving
    // arg.SetFilePrefix( usr::resultpath( "SiPMCalib", "SiPMCalc" )/"pullgraph" );
    // arg.SetNameScheme( {
    //   {"model", ""},
    //   {"nEvents", "nEvt"},
    //   {"fit", ""},
    //   {"rate", "r" },
    //   {"xtalkrate", "x"},
    //   {"darkrate", "dc"},
    //   {"alpha", "a"}
    // } );

    c.SaveAsPDF( pullplotfilename( arg,
      std::string( "pull" ) + data->GetName() ) );

    pullfile << data->GetName() << " "
             << pullmean.getVal() << " " << pullmean.getError() << " "
             << pullsig.getVal() << " " << pullsig.getError() << " "
             << npullmean.getVal() << " " << npullmean.getError() << " "
             << npullsig.getVal() << " " << npullsig.getError() << " "
             << std::endl;

    delete fit;
    delete nfit;
    delete ndata;
  }

  pullfile.close();

  return 0;
}
