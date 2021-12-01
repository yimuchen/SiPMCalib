#ifndef SIPMCALIB_SIPMCALC_SIPMWAVEFORMAT_HPP
#define SIPMCALIB_SIPMCALC_SIPMWAVEFORMAT_HPP

#include <string>
#include <vector>

/**
 * @ingroup Common
 * @brief A class for parsing standard waveform formats files.
 *
 * @details
 * The first line of the file contain exactly 3 numbers:
 *
 *  - the time interval per samples [ns]
 *  - The number of bits in the samples a single sample [N]
 *  - The conversion factor of a 1 bit to mV.
 *
 * All other lines are expected to strings of MxN hexidecimal digits representing
 * the waveformat. This calls will effectively store the format as a matrix of
 * LxM integers, representing the L waveforms extracted for the M samples in a
 * single waveform. Additional functions are provided to the common data
 * processing routines.
 */
class WaveFormat
{
public:
  WaveFormat( const std::string& file, const bool invert = true );
  ~WaveFormat();

  /**
   * @brief Getting the time interval of a single waveform sample.
   */
  inline double
  Time() const { return time; }

  /**
   * @brief Getting the number of bits for a single sample.
   */
  inline unsigned
  NBits() const { return nbits; }

  /**
   * @brief Getting the conversion factor of 1 bit to mV
   */
  inline double
  ADC() const { return adc; }

  /**
   * @brief Getting the number of waveforms collected in the file.
   */
  inline unsigned
  NWaveforms() const { return _waveforms.size();  }

  /**
   * @brief Getting the number of samples for a single sample.
   */
  inline unsigned
  NSamples() const { return _waveforms.front().size(); }

  std::vector<int16_t> WaveformRaw( const unsigned index,
                                    const int16_t  offset = 0 ) const;

  std::vector<double> Waveform( const unsigned index,
                                const unsigned pedstart = -1,
                                const unsigned pedstop  = -1 ) const;


  double WaveformSum( const unsigned index,
                      const unsigned intstart = 0,
                      const unsigned intstop  = -1,
                      const unsigned pedstart = -1,
                      const unsigned pedstop  = -1 ) const;

  std::vector<double> SumList( const unsigned intstart = 0,
                               const unsigned intstop  = -1,
                               const unsigned pedstart = -1,
                               const unsigned pedstrop = -1 ) const;

  double PedValue( const unsigned index,
                   const unsigned pedstart,
                   const unsigned pedstop ) const;
  double PedRMS( const unsigned index,
                 const unsigned pedstart,
                 const unsigned pedstop ) const;

private:
  double time;
  unsigned nbits;
  double adc;
  std::vector<std::vector<int16_t> > _waveforms;
};

#endif
