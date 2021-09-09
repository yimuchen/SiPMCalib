#ifndef SIPMCALIB_COMMON_MAKEROODATA_HPP
#define STPMCALIB_COMMON_MAKEROODATA_HPP

#include <vector>

#include <RooDataHist.h>
#include <RooRealVar.h>

extern void SetRange( RooRealVar&                var,
                      const double               binwidth,
                      const double               maxarea,
                      const std::vector<double>& datapoints );

extern RooDataHist* MakeData( RooRealVar&                var,
                              const std::vector<double>& data,
                              const double               maxarea );

#endif
