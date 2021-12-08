#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/Common/interface/Maths.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

#include <sstream>
#include <string>
#include <vector>

#include "TGraphErrors.h"
#include "TH2D.h"
#include "TProfile.h"

typedef std::vector<std::pair<double, double> > PairList;


struct ProfileMgr
{
  double _bias_min;
  double _bias_max;
  double _sipm_min;
  double _sipm_max;
  double _pulse_min;
  double _pulse_max;
  std::string _name;
  void Fill( double, double, double, double );

  TProfile*_profile;
};

std::vector<ProfileMgr> MakeProfileList( const usr::ArgumentExtender& args,
                                         const PairList&              bias_pair,
                                         const PairList&              sipm_pair,
                                         const PairList&              pulse_pair );


PairList MakePairs( const usr::ArgumentExtender& args,
                    const std::string&           arg_name );

int
main( int argc, char*argv[] )
{
  usr::po::options_description desc( "Plotting 1D stability profiles\n"
                                     "Leave one of biasrange, temppulser, or tempsipm empty for it to be used"
                                     " as the x axis" );
  desc.add_options()
    ( "data", usr::po::value<std::string>(), "Input data .txt file" )
    ( "outfile",
    usr::po::defvalue<std::string>( "profile1D.pdf" ),
    "Prefix of output file name" )
    ( "biasrange",
    usr::po::multivalue<double>(),
    "Selection range on the bias voltage" )
    ( "tempsipm",
    usr::po::multivalue<double>(),
    "Selection range on the sipm temperature" )
    ( "temppulser",
    usr::po::multivalue<double>(),
    "Section ranges on the pulser temperature" )
    ( "xmin", usr::po::defvalue<double>( 25 ), "Selection range on the x axis" )
    ( "xmax", usr::po::defvalue<double>( 35 ), "Selection range on the x axis" )
    ( "xbins", usr::po::defvalue<int>( 20 ), "Number of bins on the X axis" )
  ;
  usr::ArgumentExtender args;
  args.AddOptions( desc );
  args.ParseOptions( argc, argv );

  const std::string inputfile = args.Arg( "data" );
  const std::string outfile   = args.Arg( "outfile" );

  PairList bias_pairs   = MakePairs( args, "biasrange" );
  PairList sipm_pairs   = MakePairs( args, "tempsipm" );
  PairList pulser_pairs = MakePairs( args, "temppulser" );

  const int size_count = usr::sgn( bias_pairs.size() )
                         +usr::sgn( sipm_pairs.size() )
                         +usr::sgn( pulser_pairs.size() );
  if( size_count != 2 ){
    throw std::runtime_error(
            "Expect exactly on of the biasrange, tempsipm, or temppulsers to be empty "
            "(that variable will be used as the x axis)" );
  }

  std::vector<ProfileMgr> profile_list = MakeProfileList( args,
                                                          bias_pairs,
                                                          sipm_pairs,
                                                          pulser_pairs );

  // Reading in files
  std::string   line;
  std::ifstream fin( inputfile, std::ios::in );

  while( std::getline( fin, line ) ){
    std::istringstream linestream( line );
    double             time, readout, readouterr;
    double             bias, ledtemp, sipmtemp;
    linestream >> time >> readout >> readouterr >> bias >> ledtemp >> sipmtemp;

    for( auto& pm : profile_list ){
      pm.Fill( readout, bias, ledtemp, sipmtemp );
    }
  }

  std::cout << "Done!" << std::endl;

  const std::vector<int> colors = {
    usr::plt::col::blue, usr::plt::col::red, usr::plt::col::green,
    usr::plt::col::purple };

  usr::plt::Simple1DCanvas c;

  for( unsigned i = 0; i < profile_list.size(); ++i ){
    c.PlotHist( profile_list.at( i )._profile,
                usr::plt::TrackY( usr::plt::tracky::both ),
                usr::plt::PlotType( usr::plt::scatter ),
                usr::plt::MarkerColor( colors.at( i ) ),
                usr::plt::MarkerStyle( usr::plt::sty::mkrcircle ),
                usr::plt::MarkerSize(
                  0.5 ),
                usr::plt::LineColor( colors.at( i ) ),
                usr::plt::EntryText( profile_list.at( i )._name ) );
  }

  c.DrawLuminosity( "Stationary Gantry" );
  c.DrawCMSLabel( "Preliminary", "HGCal" );
  if( bias_pairs.size() == 1 ){
    const double mean = ( bias_pairs.front().first+bias_pairs.front().second )
                        / 2000;
    const double err = ( bias_pairs.front().second-bias_pairs.front().first )
                       / 2000;
    c.WriteLine( usr::fstr( "Bias:%.3lf#pm%.3lfV", mean, err ) );
  }
  if( sipm_pairs.size() == 1 ){
    const double mean = ( sipm_pairs.front().first+sipm_pairs.front().second )
                        / 2;
    const double err = ( sipm_pairs.front().second-sipm_pairs.front().first )
                       / 2;
    c.WriteLine( usr::fstr( "SiPM:%.1lf#pm%.2lf^{#circ}C", mean, err ) );
  }
  if( pulser_pairs.size() == 1 ){
    const double mean =
      ( pulser_pairs.front().first+pulser_pairs.front().second ) / 2;
    const double err =
      ( pulser_pairs.front().second-pulser_pairs.front().first ) / 2;
    c.WriteLine( usr::fstr( "LED:%.1lf#pm%.2lfV", mean, err ) );
  }

  c.SaveAsPDF( outfile );
  return 0;
}


PairList
MakePairs( const usr::ArgumentExtender& args,
           const std::string&           arg_name )
{
  std::vector<std::pair<double, double> > ans;
  if( args.CheckArg( arg_name ) ){
    const std::vector<double> list = args.ArgList<double>( arg_name );
    if( list.size() % 2 == 1 ){
      throw std::runtime_error( "Expected in put arguments to be in pairs!" );
    }

    for( unsigned i = 0; i < list.size() / 2; ++i ){
      ans.emplace_back( list[2 * i], list[2 * i+1] );
    }
  }

  return ans;
}


std::vector<ProfileMgr>
MakeProfileList( const usr::ArgumentExtender& args,
                 const PairList&              bias_pair,
                 const PairList&              sipm_pair,
                 const PairList&              pulse_pair )
{
  std::vector<ProfileMgr> ans;

  if( !args.CheckArg( "biasrange" ) ){
    for( const auto sp : sipm_pair ){
      for( const auto pp : pulse_pair ){
        ans.push_back( ProfileMgr() );
        ans.back()._bias_min  = 0;
        ans.back()._bias_max  = 0;
        ans.back()._sipm_min  = std::min( sp.first, sp.second );
        ans.back()._sipm_max  = std::max( sp.first, sp.second );
        ans.back()._pulse_min = std::min( pp.first, pp.second );
        ans.back()._pulse_max = std::max( pp.first, pp.second );
      }
    }
  } else if( !args.CheckArg( "tempsipm" ) ){
    for( const auto bp : bias_pair ){
      for( const auto pp : pulse_pair ){
        ans.push_back( ProfileMgr() );
        ans.back()._bias_min  = std::min( bp.first, bp.second );
        ans.back()._bias_max  = std::max( bp.first, bp.second );
        ans.back()._sipm_min  = 0;
        ans.back()._sipm_max  = 0;
        ans.back()._pulse_min = std::min( pp.first, pp.second );
        ans.back()._pulse_max = std::max( pp.first, pp.second );
      }
    }
  } else {
    for( const auto bp : bias_pair ){
      for( const auto sp : sipm_pair ){
        ans.push_back( ProfileMgr() );
        ans.back()._bias_min = std::min( bp.first, bp.second );
        ans.back()._bias_max = std::max( bp.first, bp.second );
        ans.back()._sipm_min = std::min( sp.first, sp.second );
        ans.back()._sipm_max = std::max( sp.first, sp.second );
      }
    }
  }

  for( auto& pm : ans ){
    pm._profile = new TProfile( usr::RandomString( 6 ).c_str(),
                                "",
                                args.Arg<int>( "xbins" ),
                                args.Arg<double>( "xmin" ),
                                args.Arg<double>( "xmax" ) );
    pm._profile->BuildOptions( 0, 0, "s" );
    pm._profile->GetYaxis()->SetTitle( "Readout [V-ns]" );
    if( pm._bias_min == pm._bias_max && pm._bias_min == 0 ){
      pm._profile->GetXaxis()->SetTitle( "Bias Voltage [mV]" );
    } else if( pm._sipm_min == pm._sipm_max && pm._sipm_min == 0 ){
      pm._profile->GetXaxis()->SetTitle( "SiPM Temperature [^{#circ}C]" );
    } else {
      pm._profile->GetXaxis()->SetTitle( "LED Temperature [^{#circ}C]" );
    }

    if( bias_pair.size() > 1 ){
      const double mean =  ( pm._bias_min+pm._bias_max ) / 2000;
      const double err  =  ( pm._bias_max-pm._bias_min ) / 2000;
      pm._name += usr::fstr( "  Bias:%.3lf#pm%.3lfV", mean, err );
    }
    if( sipm_pair.size() > 1 ){
      const double mean =  ( pm._sipm_min+pm._sipm_max ) / 2;
      const double err  =  ( pm._sipm_max-pm._sipm_min ) / 2;
      pm._name += usr::fstr( "  SiPM:%.1lf#pm%.2lf^{#circ}C", mean, err );
    }
    if( pulse_pair.size() > 1 ){
      const double mean =  ( pm._pulse_min+pm._pulse_max ) / 2;
      const double err  =  ( pm._pulse_max-pm._pulse_min ) / 2;
      pm._name += usr::fstr( "  LED:%.1lf#pm%.2lf^{#circ}C", mean, err );
    }
  }

  return ans;
}


void
ProfileMgr::Fill( const double readout,
                  const double bias,
                  const double ledtemp,
                  const double sipmtemp )
{
  if( _bias_min == _bias_max && _bias_min == 0 ){
    if( _pulse_min <= ledtemp  && ledtemp <= _pulse_max &&
        _sipm_min <= sipmtemp && sipmtemp <= _sipm_max ){
      _profile->Fill( bias, -readout );
    }
  } else if( _sipm_min == _sipm_max && _sipm_min == 0 ){
    if( _bias_min <= bias  && bias <= _bias_max && _pulse_min <= ledtemp  &&
        ledtemp <= _pulse_max ){
      _profile->Fill( sipmtemp, -readout );
    }
  } else {
    if( _bias_min <= bias  && bias <= _bias_max && _sipm_min <= sipmtemp &&
        sipmtemp <= _sipm_max ){
      _profile->Fill( ledtemp, -readout );
    }
  }
}
