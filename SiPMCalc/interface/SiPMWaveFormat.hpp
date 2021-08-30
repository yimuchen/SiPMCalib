#ifndef SIPMCALIB_SIPMCALC_SIPMWAVEFORMAT_HPP
#define SIPMCALIB_SIPMCALC_SIPMWAVEFORMAT_HPP

#include <string>
#include <vector>

class SiPMWaveFormat
{
public:
  SiPMWaveFormat( const std::string& file );
  ~SiPMWaveFormat();

  inline double
  Time() const { return time; }

  inline unsigned
  NBits() const { return nbits; }

  inline double
  ADC() const { return adc; }

  inline unsigned
  NWaveforms() const { return _waveforms.size();  }

  inline unsigned
  NSamples() const { return _waveforms.front().size(); }

  std::vector<int16_t> WaveformRaw( const unsigned index,
                                    const int16_t  offset = 0 ) const ;

  std::vector<double> Waveform( const unsigned index,
                                const unsigned pedstart = -1,
                                const unsigned pedstop  = -1 ) const ;


  double WaveformSum( const unsigned index,
                      const unsigned intstart = 0 ,
                      const unsigned intstop = -1 ,
                      const unsigned pedstart = -1,
                      const unsigned pedstop = -1 ) const ;

  double PedValue( const unsigned index,
                   const unsigned pedstart,
                   const unsigned pedstop ) const ;
  double PedRMS( const unsigned index,
                 const unsigned pedstart,
                 const unsigned pedstop ) const ;

private:
  double time;
  unsigned nbits;
  double adc;
  std::vector<std::vector<int16_t> > _waveforms;
};




#endif
