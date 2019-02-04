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
  MDistro(
    const double loEdge,
    const double hiEdge,
    const double epsilon,
    const double width
    );
  ~MDistro();

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
  unsigned nbins;
  std::unique_ptr<ROOT::Math::Interpolator> spline;

  std::vector<double> xArray;
  std::vector<double> convArray;
  uint64_t paramHash;

  void ParamHash();
  void MakeFFTArray();
};

class SiPMDarkFunc
{
public:
  SiPMDarkFunc(
    const double ped,
    const double gain,
    const double s0,
    const double s1
    const double dcfraction
    );
  ~SiPMDarkFunc();

  double Evaluate( const double x ) const;
  double EvalM( const double x ) const;

  void SetParam(
    const double ped,
    const double gain,
    const double s0,
    const double s1,
    const double dcfraction
    );
  double ped;
  double dcfraction;

private:
  MDistro _mdistro;
};

#endif
