#include "SiPMCalib/Common/interface/WaveFormat.hpp"

#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/Common/interface/Maths.hpp"
#include "UserUtils/PlotUtils/interface/Flat2DCanvas.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

static void MakeRawWaveform( const WaveFormat&,
                             const usr::ArgumentExtender& );

static void MakeWaveform( const WaveFormat&,
                          const usr::ArgumentExtender& );

static void MakeIntegrated( const WaveFormat&,
                            const usr::ArgumentExtender& );

static void MakePedestalPlot( const WaveFormat&,
                              const usr::ArgumentExtender& );

static void MakeOnePlot( const WaveFormat&,
                         const usr::ArgumentExtender& );

int
main( int argc, char* argv[] )
{
  usr::po::options_description desc(
    "Options for plotting a SiPM waveform from our hexadecimal data format" );
  desc.add_options()
    ( "data", usr::po::reqvalue<std::string>(), "Data .txt file" )
    ( "rawout", usr::po::value<std::string>(),
    "Output file for name raw waveform display, leave blank to skip" )
    ( "waveout", usr::po::value<std::string>(),
    "Output file name of the waveform display, leave blank to skip" )
    ( "sumout", usr::po::value<std::string>(),
    "Output file name of the integral display, leave blank to skip" )
    ( "oneout", usr::po::value<std::string>(),
    "Output filename of the single wavefor file, leave blank to skip" )
    ( "pedout", usr::po::value<std::string>(),
    "Output filename for the Pedestal value plot, leave blank to skip" )
    ( "start", usr::po::value<unsigned>()->default_value( 0 ),
    "starting sample of the waveform display" )
    ( "stop", usr::po::value<unsigned>()->default_value( 2147483647 ),
    "ending sample of the waveform display" )
    ( "pedstart", usr::po::value<unsigned>()->default_value( -1 ),
    "starting sample of the pedestal subtraction window" )
    ( "pedstop", usr::po::value<unsigned>()->default_value( -1 ),
    "ending sample of the pedestal subtraction window" )
    ( "pedrms", usr::po::value<double>()->default_value( 1000 ),
    "maximum RMS [mV] allowed in pedestal, discarding events otherwise" )
    ( "intstart", usr::po::value<unsigned>()->default_value( 0 ),
    "starting sample of the integration window" )
    ( "intstop", usr::po::value<unsigned>()->default_value( 2147483647 ),
    "ending sample of the intergration window" )
    ( "sumbinwidth", usr::po::value<double>()->default_value( 1.6 ),
    "Width of the sum histogram bins [mV-ns]" )
    ( "oneidx", usr::po::value<unsigned>()->default_value( 0 ),
    "Index of the single waveform output" );
  ;

  usr::ArgumentExtender args;
  args.AddOptions( desc );
  args.ParseOptions( argc, argv );

  // Making sure that there is some data to generate.
  if( !args.CheckArg( "data" ) ){
    std::cerr << "Please provide input data" << std::endl;
    return 1;
  }

  // Making the raw data format container
  WaveFormat wformat( args.Arg<std::string>( "data" ) );

  // Running the separate sub commands
  if( args.CheckArg( "rawout" ) ){
    MakeRawWaveform( wformat, args );
  }

  if( args.CheckArg( "waveout" ) ){
    MakeWaveform( wformat, args );
  }

  if( args.CheckArg( "sumout" ) ){
    MakeIntegrated( wformat, args );
  }

  if( args.CheckArg( "pedout" ) ){
    MakePedestalPlot( wformat, args );
  }

  if( args.CheckArg( "oneout" ) ){
    MakeOnePlot( wformat, args );
  }
  return 0;

}

void
MakeRawWaveform( const WaveFormat&            wformat,
                 const usr::ArgumentExtender& args )
{
  // Adding the additional parsing arguments
  const unsigned start = std::min( args.Arg<unsigned>( "start" )
                                 , wformat.NSamples() );
  const unsigned stop = std::min( args.Arg<unsigned>( "stop" )
                                , wformat.NSamples() );
  const std::string output = args.Arg<std::string>( "rawout" );

  int16_t ymin = 0;
  int16_t ymax = 0;

  for( unsigned i = 0; i < wformat.NWaveforms(); ++i ){
    const auto waveform = wformat.WaveformRaw( i );

    for( unsigned j = start; j < stop; ++j ){
      ymin = std::min( waveform.at( j ), ymin );
      ymax = std::max( waveform.at( j ), ymax );
    }
  }

  if( ( ymax - ymin )%2 == 0 ){ ymax++;}

  const unsigned timebins = ( stop - start )/2;
  const unsigned ybins    = ( ymax - ymin )/2;

  // Making the histogram
  TH2D hist( "hist", "hist",
             timebins, start, stop,
             ybins, ymin, ymax );

  // Filling in the histogram
  for( unsigned i = 0; i < wformat.NWaveforms(); ++i ){
    const auto waveform = wformat.WaveformRaw( i );

    for( unsigned j = start; j < stop; ++j ){
      hist.Fill( j, waveform[j] );
    }
  }

  // Setting the zero bins for aesthetics purposes
  for( int i = 0; i < hist.GetNcells(); ++i ){
    if( hist.GetBinContent( i ) == 0 ){
      hist.SetBinContent( i, 0.3 );// Half a order of magnitude smaller
    }
  }

  // Plotting this histogram;
  usr::plt::Flat2DCanvas c;
  c.PlotHist( hist,
    usr::plt::Plot2DF( usr::plt::heat ) );

  c.DrawCMSLabel( "Preliminary", "HGCal" );
  c.DrawLuminosity( "Pulser Setup" );

  c.Pad().SetTextCursor( 0.05, 0.9, usr::plt::font::top_left );

  c.Xaxis().SetTitle(
    usr::fstr( "Sample Index [%.2lfns]", wformat.Time() ).c_str() );
  c.Yaxis().SetTitle(
    usr::fstr( "Readout [%.2lfmV]", wformat.ADC() ).c_str() );
  c.Zaxis().SetTitle( "Number of data points" );

  c.SaveAsPDF( output );
}

void
MakeWaveform( const WaveFormat&            wformat,
              const usr::ArgumentExtender& args )
{
  // Adding the additional parsing arguments
  const unsigned start = std::min( args.Arg<unsigned>( "start" )
                                 , wformat.NSamples() );
  const unsigned stop = std::min( args.Arg<unsigned>( "stop" )
                                , wformat.NSamples() );
  const unsigned pedstart  = args.Arg<unsigned>( "pedstart" );
  const unsigned pedstop   = args.Arg<unsigned>( "pedstop" );
  const double pedrms      = args.Arg<double>( "pedrms" );
  const std::string output = args.Arg<std::string>( "waveout" );

  // To avoid binning artifacts, the histogram will be displayed in raw sample
  // index and ADC bins
  double ymin = 0;
  double ymax = 0;

  for( unsigned i = 0; i < wformat.NWaveforms(); ++i ){
    const auto waveform = wformat.WaveformRaw( i );

    for( unsigned j = start; j < stop; ++j ){
      ymin = std::min( waveform[j] * wformat.ADC(), ymin );
      ymax = std::max( waveform[j] * wformat.ADC(), ymax );
    }
  }

  // Making the axis based on the result parsing. Making the binning scheme
  // slightly different to the actual bin width to avoid rounding artifacts.
  const unsigned timebins = ( stop - start ) /2;
  const unsigned ybins    = ( ymax-ymin ) / wformat.ADC() / 2.0;

  // Making the histogram
  TH2D hist( "hist", "hist",
             timebins, start, stop,
             ybins, ymin, ymax );

  // Filling in the histogram
  for( unsigned i = 0; i < wformat.NWaveforms(); ++i ){
    // Skip events with large pedestal RMS
    if( wformat.PedRMS( i, pedstart, pedstop ) > pedrms ){
      continue;
    }

    const auto waveform = wformat.Waveform( i, pedstart, pedstop );
    double sum          = 0;

    for( unsigned j = start; j < stop; ++j ){
      sum += waveform[j];
      hist.Fill( j, waveform[j] );
    }
  }

  // Setting the zero bins for aesthetics purposes
  for( int i = 0; i < hist.GetNcells(); ++i ){
    if( hist.GetBinContent( i ) == 0 ){
      hist.SetBinContent( i, 0.3 );// Half a order of magnitude smaller
    }
  }

  // Plotting this histogram;
  usr::plt::Flat2DCanvas c;
  c.PlotHist( hist,
    usr::plt::Plot2DF( usr::plt::heat ) );

  c.DrawCMSLabel( "Preliminary", "HGCal" );
  c.DrawLuminosity( "Pulser Setup" );

  c.Pad().SetTextCursor( 0.05, 0.9, usr::plt::font::top_left );

  c.Xaxis().SetTitle(
    usr::fstr( "Sample index [%.3lfns]", wformat.Time() ).c_str() );
  c.Yaxis().SetTitle( "Readout [mV]" );
  c.Zaxis().SetTitle( "Number of data points" );

  // c.SetLogz( 1 );

  c.SaveAsPDF( output );
}

static void
MakeOnePlot( const WaveFormat&            wformat,
             const usr::ArgumentExtender& args )
{
  // Adding the additional parsing arguments
  const unsigned start = std::min( args.Arg<unsigned>( "start" )
                                 , wformat.NSamples() );
  const unsigned stop = std::min( args.Arg<unsigned>( "stop" )
                                , wformat.NSamples() );
  const unsigned pedstart = args.Arg<unsigned>( "pedstart" );
  const unsigned pedstop  = args.Arg<unsigned>( "pedstop" );
  const unsigned idx      = args.Arg<unsigned>( "oneidx" );
  // const unsigned intstart  = args.Arg<unsigned>( "intstart" );
  // const unsigned intstop   = args.Arg<unsigned>( "intstop" );
  const std::string output = args.Arg<std::string>( "oneout" );


  // Making the waveform graph
  TGraph graph( stop-start );

  // Filling in the histogram
  const auto waveform = wformat.Waveform( idx, pedstart, pedstop );

  for( unsigned j = start; j < stop; ++j ){
    graph.SetPoint( j,
      j*wformat.Time(),
      waveform[j] );
  }

  // Plotting this histogram;
  usr::plt::Simple1DCanvas c;
  c.PlotGraph( graph,
    usr::plt::PlotType( usr::plt::simplefunc ) );

  c.DrawCMSLabel( "Preliminary", "HGCal" );
  c.DrawLuminosity( "Laser setup" );

  c.Xaxis().SetTitle( "Time since trigger [ns]" );
  c.Yaxis().SetTitle( "Readout [mV]" );

  c.SaveAsPDF( output );
}

static void
MakeIntegrated( const WaveFormat&            wformat,
                const usr::ArgumentExtender& args )
{
  // Adding the additional parsing arguments
  const unsigned pedstart  = args.Arg<unsigned>( "pedstart" );
  const unsigned pedstop   = args.Arg<unsigned>( "pedstop" );
  const double pedrms      = args.Arg<double>( "pedrms" );
  const unsigned intstart  = args.Arg<unsigned>( "intstart" );
  const unsigned intstop   = args.Arg<unsigned>( "intstop" );
  const std::string output = args.Arg<std::string>( "sumout" );
  const double binwidth    = args.Arg<double>( "sumbinwidth" );


  double min = 0;
  double max = 0;
  std::vector<double> vals;
  unsigned discarded = 0;

  for( unsigned i = 0; i < wformat.NWaveforms(); ++i ){
    if( wformat.PedRMS( i, pedstart, pedstop ) > pedrms ){
      discarded++;
      continue;
    }

    vals.push_back( wformat.WaveformSum( i, intstart, intstop
                                       , pedstart, pedstop ) );
    min = std::min( vals.back(), min );
    max = std::max( vals.back(), max );
  }

  const double xmin  = usr::RoundDown( min, binwidth );
  const double xmax  = usr::RoundUp( max, binwidth );
  const double nbins = ( xmax - xmin ) / binwidth;

  TH1D hist( "hist1d", "hist1d", nbins, xmin, xmax );

  for( const auto v : vals ){
    hist.Fill( v );
  }

  // Simple message for discarded values
  if( discarded > 0 ){
    usr::fout( "Discarded %d out of %d events\n"
             , discarded, wformat.NWaveforms() );
  }

  usr::plt::Simple1DCanvas c;
  c.PlotHist( hist, usr::plt::PlotType( usr::plt::hist ) );
  c.Xaxis().SetTitle( "Readout [mV-ns]" );
  c.Yaxis().SetTitle( "Events" );

  c.SaveAsPDF( output );
}

static void
MakePedestalPlot( const WaveFormat&            wformat,
                  const usr::ArgumentExtender& args )
{
  // Adding the additional parsing arguments
  const unsigned pedstart  = args.Arg<unsigned>( "pedstart" );
  const unsigned pedstop   = args.Arg<unsigned>( "pedstop" );
  const double pedrms      = args.Arg<double>( "pedrms" );
  const std::string output = args.Arg<std::string>( "pedout" );

  TH1D hist( "pedhist", "pedhist", 25, -5, 5 );

  for( unsigned i = 0; i < wformat.NWaveforms(); ++i ){
    if( wformat.PedRMS( i, pedstart, pedstop ) > pedrms ){
      continue;
    }
    const auto w = wformat.Waveform( i );

    for( unsigned j = pedstart; j < pedstop; ++j  ){
      hist.Fill( w.at( j ) );
    }
  }

  usr::plt::Simple1DCanvas c;
  c.PlotHist( hist
            , usr::plt::PlotType( usr::plt::hist )
            , usr::plt::EntryText( usr::fstr( "RMS %.2lf", hist.GetRMS() ) ) );
  c.Xaxis().SetTitle( "Readout value in pedestal window [mV]" );
  c.Yaxis().SetTitle( "Number of cells" );

  c.SaveAsPDF( output );

}
