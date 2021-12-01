#include "SiPMCalib/Common/interface/WaveFormat.hpp"

#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/Maths.hpp"
#include "UserUtils/Common/interface/STLUtils/OStreamUtils.hpp"
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"
#else
#include "UserUtils/Common/Maths.hpp"
#include "UserUtils/Common/STLUtils/OStreamUtils.hpp"
#include "UserUtils/Common/STLUtils/StringUtils.hpp"
#endif

#include <fstream>
#include <sstream>

static inline int16_t
hex_to_int( const char x )
{
  return x > 'a' ? 10 + x - 'a' :
         x > 'A' ? 10 + x - 'A' :
         x - '0';
}

static inline int8_t
bit4_to_bit2( const int16_t x )
{
  return x;
}

/**
 * @brief Construction of the waveform array form a given input file.
 *
 * The user can also specify where the waveform should be inverted or not (change
 * downward firing pulse to be upwards firing pulses).
 *
 * For the data collection of the DRS4, there are occasionally bad samples due to
 * bit flips in the ADC chip. Here we filter out these single bit flips by
 * checking the neighboring cells up to 2 samples away in either direction, and
 * check if the maximum variation in the sample value is more than 70 bits
 * (nothing in the system is expected to be this fast). If a sample readout is
 * determined to contain a bit flip, then the readout is replace with the average
 * of the readout before and after the bad sample.
 */
WaveFormat::WaveFormat( const std::string& file, const bool invert )
{
  const unsigned factor = invert ? -1 : 1;
  std::string line;
  std::ifstream fin( file, std::ios::in );

  if( !fin.is_open() ){
    usr::log::PrintLog( usr::log::FATAL,// Throws exception
      usr::fstr( "Input file [%s] cannot be opened!", file ) );
  }

  // Getting first line
  std::getline( fin, line );
  std::istringstream linestream( line );
  linestream >> time >> nbits >> adc;

  // Getting all other lines
  while( std::getline( fin, line ) ){
    _waveforms.push_back( std::vector<int16_t>( line.length() / nbits ) );

    for( unsigned index = 0; index < _waveforms.back().size(); ++index ){
      int16_t value = 0;

      for( unsigned bit = 0; bit < nbits; ++bit ){
        const int16_t bit_value = hex_to_int( line[nbits * index + bit] );
        value = value << 4 | bit_value;
      }

      _waveforms.back()[index] = nbits == 4 ? factor * value :
                                 bit4_to_bit2( value * factor );
    }

    auto is_peak_cell
      = [this]( const unsigned index )->bool {
          const auto& w         = this->_waveforms.back();
          const unsigned diffp1 = index > w.size() - 2 ? 0 :
                                  abs( w.at( index ) - w.at( index+1 ) );
          const unsigned diffp2 = index > w.size() -3 ? 0 :
                                  abs( w.at( index ) - w.at( index+2 ) );
          const unsigned diffm1 = index < 1 ? 0 :
                                  abs( w.at( index ) - w.at( index-1 ) );
          const unsigned diffm2 = index < 2 ? 0 :
                                  abs( w.at( index ) - w.at( index-2 ) );

          const unsigned diffp = std::max( diffp1, diffp2 );
          const unsigned diffm = std::max( diffm1, diffm2 );

          return ( diffp > 70 ) && ( diffm > 70 );
        };


    for( unsigned index = 0; index < _waveforms.back().size(); ++index ){
      if( is_peak_cell( index ) ){
        _waveforms.back()[index]
          = ( _waveforms.back()[index+1] + _waveforms.back()[index-1] ) / 2;
      }
    }
  }
}


WaveFormat::~WaveFormat(){}

/**
 * @brief Getting the original (integer) waveform at some certain index.
 *
 * The user can specify a certain offset to be applied for the entire waveform.
 */
std::vector<int16_t>
WaveFormat::WaveformRaw( const unsigned index, const int16_t offset ) const
{
  std::vector<int16_t> ans;

  for( const auto val : _waveforms.at( index ) ){
    ans.push_back( val - offset );
  }

  return ans;
}

/**
 * @brief Returning the average voltage value of samples within the given
 * pedestal window. If you set this start and stop index values to be -1
 * (default), then the pedestal value would be strictly 0. The value return will
 * be in units of mV.
 */
double
WaveFormat::PedValue( const unsigned index,
                      const unsigned pedstart,
                      const unsigned pedstop ) const
{
  if( pedstart == unsigned(-1) || pedstop == unsigned(-1) ){
    return 0;
  }

  double ped_value = 0;

  for( unsigned i = pedstart; i < pedstop; ++i ){
    ped_value += _waveforms.at( index ).at( i ) * ADC();
  }

  return ped_value / (double)( pedstop - pedstart );
}


/**
 * @brief Returning the standard deviation of the samples at a given event index
 * within the given pedestal window.
 */
double
WaveFormat::PedRMS( const unsigned index,
                    const unsigned pedstart,
                    const unsigned pedstop ) const
{
  // Returns zero is the default empty range is given.
  if( pedstart == unsigned(-1) || pedstop == unsigned(-1) ){
    return 0;
  }

  std::vector<double> list;

  for( unsigned i = pedstart; i < pedstop; ++i ){
    list.push_back( _waveforms.at( index ).at( i ) * ADC() );
  }

  return usr::StdDev( list );
}

/**
 * @brief Getting a stored waveform after converting to readout values.
 *
 * If the pedestal start and stop index is specified, then the pedestal value is
 * subtracted the entire waveform values. See WaveFormat::PedValue to see how the
 * pedestal value is calculated.
 */
std::vector<double>
WaveFormat::Waveform( const unsigned index,
                      const unsigned pedstart,
                      const unsigned pedstop ) const
{
  std::vector<double> ans;
  const double ped_value = PedValue( index, pedstart, pedstop );

  for( unsigned i = 0; i < NSamples(); ++i ){
    ans.push_back( _waveforms.at( index ).at( i ) * ADC() - ped_value );
  }

  return ans;
}

/**
 * @brief Calculating the area of a waveform given some intergration window.
 *
 * The return value would be in units of mV-ns. The temporal time would be
 * multiplied to the summation results. Here the user can also specify the
 * pedestal window for pedestal subtraction if needed.
 */
double
WaveFormat::WaveformSum( const unsigned index,
                         const unsigned intstart,
                         const unsigned intstop,
                         const unsigned pedstart,
                         const unsigned pedstop ) const
{
  double ans             = 0;
  const double ped_value = PedValue( index, pedstart, pedstop );

  const unsigned start = std::max( intstart, (unsigned)0 );
  const unsigned stop  = std::min( intstop, NSamples() );

  for( unsigned i = start; i < stop; ++i ){
    ans += _waveforms.at( index ).at( i )*ADC()  - ped_value;
  }

  ans *= Time();

  return ans;
}

/**
 * @brief Getting all waveform areas given the integration window and optional.
 * pedestal subtraction window.
 *
 * The output will be sorted according to area.
 */
std::vector<double>
WaveFormat::SumList( const unsigned intstart,
                     const unsigned intstop,
                     const unsigned pedstart,
                     const unsigned pedstop ) const
{
  std::vector<double> ans;
  ans.reserve( NWaveforms() );

  for( unsigned i = 0; i < NWaveforms(); ++i ){
    ans.push_back( WaveformSum( intstart, intstop, pedstart, pedstop ) );
  }

  std::sort( ans.begin(), ans.end() );

  return ans;
}
