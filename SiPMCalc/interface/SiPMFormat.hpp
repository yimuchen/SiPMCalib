#ifndef SIPMCALIB_SIPMCALC_SIPMFORMAT_HPP
#define SIPMCALIB_SIPMCALC_SIPMFORMAT_HPP

#include <memory>
#include <string>
#include <vector>

#include "RooDataHist.h"
#include "RooRealVar.h"

class SiPMFormat
{
public:
  SiPMFormat( const std::string& file,
              const double       binwdith,
              const unsigned     start = 0,// starting samples
              const unsigned     end   = -1,// ending samples,
              const std::string& baseline =""
              );
  ~SiPMFormat();

  inline RooRealVar&
  x(){ return *_x; }
  inline RooAbsData&
  data(){ return *_data; }

  // Parameter Return
  inline unsigned
  TimeInterval() const { return _timeint; }
  inline unsigned
  NCaptures() const { return _ncapture; }
  inline unsigned
  PreSamples() const { return _presample; }
  inline unsigned
  PostSamples() const { return _postsample; }
  inline double
  ADCConversion() const { return _convfactor * 256; }

  // Making the data set
  void TruncateDataSet( const double maxarea = 0 );

  inline unsigned
  NArea() const { return _arealist.size(); }
  inline double
  Area( const unsigned i ){ return _arealist.at( i ) ;  }
  inline std::vector<double>&
  AreaList() { return _arealist; }

private:
  std::vector<double> _arealist;
  std::vector<double> _baseline;
  std::vector<double> _baseline_err;
  std::unique_ptr<RooAbsData> _data;
  std::unique_ptr<RooRealVar> _x;

  double _binwidth;
  unsigned _timeint;
  unsigned _ncapture;
  unsigned _presample;
  unsigned _postsample;
  double _convfactor;

  void MakeBaseLine( const std::string& );
};




#endif
