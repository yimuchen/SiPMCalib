#ifndef SIPMCALIB_SIPMCALC_SIPMDARKFUNC_HPP
#define SIPMCALIB_SIPMCALC_SIPMDARKFUNC_HPP

#include <cstdint>
#include <vector>
#include <memory>
#include "Math/Interpolator.h"

class MDistro {
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
  double xMin() const ;
  double xMax() const ;
  double EdgeDist() const;

  double MFuncEval( const double x ) const ;

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

  void   ParamHash() ;
  void   MakeFFTArray() ;
};

class SiPMDarkFunc
{
public:
  SiPMDarkFunc(
    const double ped,
    const double gain,
    const double s0,
    const double s1,
    const double acfrac1,
    const double acfrac2,
    const double acshift
  );
  ~SiPMDarkFunc();

  double Evaluate( const double x ) const ;
  double EvalM1( const double x ) const ;
  double EvalM2( const double x ) const ;

  void SetParam(
    const double ped,
    const double gain,
    const double s0,
    const double s1,
    const double acfrac1,
    const double acfrac2,
    const double acshift
  );
  double ped;
  double gain;
  double s0;
  double s1;
  double acfrac1;
  double acfrac2;
  double acshift;

private:
  MDistro _m_primary;
  MDistro _m_secondary;

  double low_edge()  const;
  double high_edge() const;
  double epsilon()   const;
  double w1()        const;
};

#endif