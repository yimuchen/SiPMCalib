#ifndef SIPMCALIB_SIPMCALC_SIPMDARKFUNC_HPP
#define SIPMCALIB_SIPMCALC_SIPMDARKFUNC_HPP

#include "Math/Interpolator.h"
#include <cstdint>
#include <memory>
#include <vector>

class MDistro
{
public:
  MDistro();
  MDistro( const double loEdge,
           const double hiEdge,
           const double epsilon,
           const double width );
  ~MDistro();

  double EvaluateAccum( const double x ) const;
  double Evaluate( const double x ) const;
  void   SetParam( const double, const double, const double, const double );

  // Getting fitting parameters
  double xMin() const;
  double xMax() const;
  double EdgeDist() const;

  double MFuncEval( const double x ) const;

  double loEdge;
  double hiEdge;
  double epsilon;
  double width;

private:
  unsigned                 nbins;
  ROOT::Math::Interpolator spline;
  ROOT::Math::Interpolator spline_acc;

  std::vector<double> xArray;
  std::vector<double> convArray;
  std::vector<double> accArray;

  // Delaring here, as we are going to reserve for faster allocation
  std::vector<double> convTempArray;
  std::vector<double> mfuncArray;
  std::vector<double> gaussArray;
  uint64_t            paramHash;

  void ParamHash();
  void MakeFFTArray();
};

#endif
