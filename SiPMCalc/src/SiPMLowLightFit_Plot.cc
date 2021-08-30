#include "SiPMCalib/SiPMCalc/interface/SiPMLowLightFit.hpp"
#include "UserUtils/PlotUtils/interface/Ratio1DCanvas.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"


void
SiPMLowLightFit::PlotPeakFind( const std::string& output )
{
  usr::plt::Simple1DCanvas c;

  c.PlotHist( _est_hist.get(),
    usr::plt::PlotType( usr::plt::hist ),
    usr::plt::EntryText( "Readout" ),
    usr::plt::LineColor( usr::plt::col::black ) );

  for( auto& f : _peakfits ){
    c.PlotFunc( f.get(),
      usr::plt::PlotType( usr::plt::simplefunc ),
      usr::plt::LineColor( usr::plt::col::red ),
      f == _peakfits.front() ? usr::plt::EntryText( "Local Gaussian fit" ) :
      RooCmdArg::none() );
  }

  for( int i = 0; i < _spectrum->GetNPeaks(); ++i ){
    auto& line = c.Pad().DrawVLine(
      _spectrum->GetPositionX()[i],
      usr::plt::LineColor( usr::plt::col::darkgray ),
      usr::plt::LineStyle( usr::plt::sty::lindotted ) );
    if( i == 0 ){
      c.Pad().AddLegendEntry( line, "Peak search results", "L" );
    }
  }

  c.DrawCMSLabel( "Peak finding", "Spectral Fit" );
  c.Pad().Xaxis().SetTitle( "Readout" );
  c.Pad().Yaxis().SetTitle( "Events" );
  c.SaveAsPDF( output );
}

void
SiPMLowLightFit::PlotGainFit( const std::string& output )
{
  usr::plt::Simple1DCanvas c;

  c.PlotGraph( _gain_graph.get(),
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::TrackY( usr::plt::tracky::both ),
    usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
    usr::plt::MarkerSize( 0.5 ) );

  c.PlotFunc( _gain_fit.get(),
    usr::plt::PlotType( usr::plt::simplefunc ),
    usr::plt::LineColor( usr::plt::col::red ) );

  c.Xaxis().SetTitle( "Discharge peak" );
  c.Yaxis().SetTitle( "Fitted peak position [ADC #times ns]" );
  c.DrawCMSLabel( "Spectral Fit" );
  c.Pad().WriteLine( usr::fstr( "Est. gain=%.0lf#pm%.1f"
                              , _gain_fit->GetParameter( 0 )
                              , _gain_fit->GetParError( 0 ) ) );
  c.Pad().WriteLine( usr::fstr( "Est. pedestal=%.0lf#pm%.1f"
                              , _gain_fit->GetParameter( 0 )
                              , _gain_fit->GetParError( 0 ) ) );

  c.SaveAsPDF( output );
}

void
SiPMLowLightFit::PlotWidthFit( const std::string& output )
{
  usr::plt::Simple1DCanvas c;

  c.PlotGraph( _width_graph.get(),
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::TrackY( usr::plt::tracky::both ),
    usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
    usr::plt::MarkerSize( 0.5 ) );
  c.PlotFunc( _width_fit.get(),
    usr::plt::PlotType( usr::plt::simplefunc ),
    usr::plt::LineColor( usr::plt::col::red ) );

  c.Xaxis().SetTitle( "Discharge peak" );
  c.Yaxis().SetTitle( "Fitted peak width [ADC #times ns]" );
  c.DrawCMSLabel( "Spectral Fit" );
  c.Pad().WriteLine( usr::fstr( "Est. s_{0}=%.0lf#pm%.1f"
                              , _width_fit->GetParameter( 0 )
                              , _width_fit->GetParError( 0 ) ) );
  c.Pad().WriteLine( usr::fstr( "Est. s_{1}=%.0lf#pm%.1f"
                              , _width_fit->GetParameter( 1 )
                              , _width_fit->GetParError( 1 ) ) );

  c.SaveAsPDF( output );
}

void
SiPMLowLightFit::PlotPoissonFit( const std::string& output )
{
  usr::plt::Simple1DCanvas c;

  TH1D hist( usr::RandomString( 5 ).c_str(), "",
             _height_graph->GetN(), -0.5, _height_graph->GetN()-0.5 );

  for( int i = 0; i < _height_graph->GetN(); ++i ){
    hist.Fill( i, _height_fit->Eval( i ) );
  }

  c.PlotHist( hist,
    usr::plt::PlotType( usr::plt::hist ),
    usr::plt::LineColor( usr::plt::col::red ),
    usr::plt::EntryText( "Poisson fit" ) );

  c.PlotGraph( _height_graph.get(),
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::EntryText( "Local gaussian norm" ),
    usr::plt::TrackY( usr::plt::tracky::both ),
    usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
    usr::plt::MarkerSize( 1 ) );


  c.Xaxis().SetTitle( "Discharge peak" );
  c.Yaxis().SetTitle( "Fitted peak height [Events]" );
  c.DrawCMSLabel( "Spectral Fit" );
  c.Pad().WriteLine( usr::fstr( "Est. mean=%.2lf#pm%.3f"
                              , _height_fit->GetParameter( 1 )
                              , _height_fit->GetParError( 1 ) ) );
  c.Pad().WriteLine( usr::fstr( "Est. #lambda=%.3lf#pm%.4f"
                              , _height_fit->GetParameter( 2 )
                              , _height_fit->GetParError( 2 ) ) );

  c.SaveAsPDF( output );
}

void
SiPMLowLightFit::PlotSpectrumFit( const std::string& output )
{
  usr::plt::Ratio1DCanvas c( x() );

  auto& fitgraph = c.PlotPdf( _pdf.get(),
    RooFit::Normalization( _data->sumEntries(), RooAbsReal::NumEvent ),
    usr::plt::EntryText( "Fit" ),
    usr::plt::PlotType( usr::plt::simplefunc ),
    usr::plt::LineColor( usr::plt::col::blue ) );

  auto& datgraph = c.PlotData( _data.get(),
    usr::plt::EntryText( "Data" ),
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::MarkerSize( 0.2 ) );

  c.PlotScale( fitgraph, fitgraph,
    usr::plt::PlotType( usr::plt::simplefunc ) );

  c.PlotScale( datgraph, fitgraph,
    usr::plt::PlotType( usr::plt::scatter ),
    usr::plt::MarkerSize( 0.2 ) );

  c.TopPad().DrawLuminosity( _lumitype );
  c.TopPad().DrawCMSLabel( "", "Spectral fit" );

  // Writing the typical stuff
  if( _sipmtype != "" ){
    c.TopPad().WriteLine( _sipmtype );
  }

  if( _biasv != ""  ){
    c.TopPad().WriteLine( _biasv );
  }

  c.BottomPad().Yaxis().SetTitle( "Data/Global fit" );
  c.SetLogy( true );
  c.TopPad().SetYaxisMax( c.TopPad().GetYaxisMax() );

  c.SaveAsPDF( output );
}

void
SiPMLowLightFit::PrintRaw( std::ostream& sout ) const
{
  static const std::string rawfmt( "%10s %10.5lf %10.5lf\n" );

  auto make_rooline = [&rawfmt]( const std::string& name, const RooRealVar& x )->std::string {
                        return usr::fstr( rawfmt, name
                                        , x.getVal(), x.getError() );
                      };
  auto make_usrline = [&rawfmt]( const std::string& name
                               , const usr::Measurement& x )->std::string {
                        return usr::fstr( rawfmt, name
                                        , x.CentralValue(), x.AbsAvgError() );
                      };

  sout << "Raw Fit parameter results\n";
  sout << make_rooline( "ped",     ped() );
  sout << make_rooline( "gain",    gain() );
  sout << make_rooline( "s0",      s0() );
  sout << make_rooline( "s1",      s1() );
  sout << make_rooline( "mean",    mean() );
  sout << make_rooline( "lambda",  lambda() );
  sout << make_rooline( "alpha",   alpha() );
  sout << make_rooline( "beta",    beta() );
  sout << make_rooline( "dcfrac",  dcfrac() );
  sout << make_rooline( "epsilon", eps() );

  // "Nominal" method for measuring ENF
  const usr::Measurement enf = ExcessNoiseFactor( MeanPhotons() );
  sout << usr::fstr( rawfmt, "ENF",     enf.CentralValue(), enf.AbsAvgError() );

  // "Conventional" method for measuring ENF
  const usr::Measurement m    = MeanPhotonsFromFrac();
  const usr::Measurement enff = ExcessNoiseFactor( m );
  sout << make_usrline( "mfrac",  m     );
  sout << make_usrline( "ENFfrac", enff );
}

void
SiPMLowLightFit::PrintTable( std::ostream& sout ) const
{
  static const std::string fmt    = "%50s & = & %20s & [\\text{%s}] \\\\\n";
  static const std::string mtitle = "\\left\\langle N_{\\gamma}\\right\\rangle";

  auto mkstr = []( const usr::Measurement& x ){
                 return usr::fstr( "%.2lf \\pm %.3lf"
                                 , x.CentralValue()
                                 , x.AbsAvgError() );
               };

  sout << "Human readable fit results\n";
  sout << usr::fstr( fmt, mtitle,                mkstr( MeanPhotons() ),                     "Photons" );
  sout << usr::fstr( fmt, "\\text{Gain}",        mkstr( Gain() ),                            "mV-ns"   );
  sout << usr::fstr( fmt, "\\sigma_\\text{com}", mkstr( CommonNoise() ),                     "mV-ns"   );
  sout << usr::fstr( fmt, "\\sigma_\\text{pix}", mkstr( PixelNoise() ),                      "mV-ns"   );
  sout << usr::fstr( fmt, "P_\\text{ct}",        mkstr( ProbCrosstalk()*100.0 ),             "\\%"     );
  sout << usr::fstr( fmt, "P_\\text{ap}",        mkstr( ProbAfterpulse()*100.0 ),            "\\%"     );
  sout << usr::fstr( fmt, "\\tau_\\text{ap}",    mkstr( AfterpulseTimeNS() ),                "ns"      );
  sout << usr::fstr( fmt, "\\tau_\\text{dc}",    mkstr( DarkcurrentTimeNS() / 1000.0 ),      "\\mus"   );
  sout << usr::fstr( fmt, "ENF",                 mkstr( ExcessNoiseFactor( MeanPhotons() ) ),""        );

}
