#ifndef SIPMCALIB_SIPMCALC_SIPMPDF_FUNCTOR_HPP
#define SIPMCALIB_SIPMCALC_SIPMPDF_FUNCTOR_HPP

#include <Math/IFunction.h>

class SiPMPdf_Functor : public ROOT::Math::IBaseFunctionMultiDim
{
public:
  SiPMPdf_Functor();
  ~SiPMPdf_Functor();
  virtual ROOT::Math::IBaseFunctionMultiDim* Clone() const;
  inline unsigned
  NDim() const { return 9; }

  double g_poisson( const int k, const double* _x ) const;

  double after_pulse_eff(
    const int     k,
    const int     i,
    const double* _x
    ) const;

  double gauss_k( const int k, const double* _x ) const;

private:
  double DoEval( const double* _x ) const;

  static inline double
  x     ( const double* _x ){ return _x[0];}
  static inline double
  ped   ( const double* _x ){ return _x[1];}
  static inline double
  gain  ( const double* _x ){ return _x[2];}
  static inline double
  s0    ( const double* _x ){ return _x[3];}
  static inline double
  s1    ( const double* _x ){ return _x[4];}
  static inline double
  mean  ( const double* _x ){ return _x[5];}
  static inline double
  lambda( const double* _x ){ return _x[6];}
  static inline double
  alpha ( const double* _x ){ return _x[7];}
  static inline double
  beta  ( const double* _x ){ return _x[8];}

};

#endif
