#ifndef SIPMCALIB_SIPMCALC_SIPMSUMGAUS_FUNCTOR_HPP
#define SIPMCALIB_SIPMCALC_SIPMSUMGAUS_FUNCTOR_HPP

#include <Math/IFunction.h>

class SiPMSumGaus_Functor : public ROOT::Math::IBaseFunctionMultiDim
{
public:
  SiPMSumGaus_Functor( const unsigned x );
  ~SiPMSumGaus_Functor();

  virtual ROOT::Math::IBaseFunctionMultiDim* Clone() const;

  inline unsigned
  NDim() const { return _ndim; }
  inline unsigned
  NPeak() const { return _ndim - 5 ; }

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
  amp  ( const double* _x, const int i ){ return _x[5+i];}

  const unsigned _ndim;
};

#endif
