#ifndef SIPMCALIB_SIPMCALC_SIPMPDF_HPP
#define SIPMCALIB_SIPMCALC_SIPMPDF_HPP

#include "SiPMCalib/SiPMCalc/interface/SiPMDarkFunc.hpp"

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooTrace.h"

#include <vector>

class SiPMPdf : public RooAbsPdf
{
public:
  // Full model constructor
  SiPMPdf( const char*, const char*,
           RooAbsReal& x,
           RooAbsReal& ped,
           RooAbsReal& gain,
           RooAbsReal& s0,
           RooAbsReal& s1,
           RooAbsReal& mean,
           RooAbsReal& lambda,
           RooAbsReal& alpha,
           RooAbsReal& beta,
           RooAbsReal& dcfraction,
           RooAbsReal& epsilon
           );

  // No dark current
  SiPMPdf( const char*, const char*,
           RooAbsReal& x,
           RooAbsReal& ped,
           RooAbsReal& gain,
           RooAbsReal& s0,
           RooAbsReal& s1,
           RooAbsReal& mean,
           RooAbsReal& lambda,
           RooAbsReal& alpha,
           RooAbsReal& beta
           );

  // No after pulsing
  SiPMPdf( const char*, const char*,
           RooAbsReal& x,
           RooAbsReal& ped,
           RooAbsReal& gain,
           RooAbsReal& s0,
           RooAbsReal& s1,
           RooAbsReal& mean,
           RooAbsReal& lambda
           );

  SiPMPdf( const SiPMPdf&, const char* name = 0 );

  virtual ~SiPMPdf();

  virtual TObject* clone( const char* name ) const;

  inline double
         Eval() const { return evaluate(); }
  double gen_poisson( const int k ) const;
  double ap_eff( const int k, const int i ) const;
  double gauss_k( const int k  ) const;
  double binomial_prob( const int k, const int i ) const;

  int getAnalyticalIntegral( RooArgSet&  allVars,
                             RooArgSet&  analVars,
                             const char* rangeName = 0 ) const override;
  double analyticalIntegral( int         code,
                             const char* rangeName = 0 ) const override;

  double analyticalIntegral( double x ) const;
  double erf_ap_eff( const double x, const int k, const int i ) const;
  double erf_k( const double x, const int k ) const ;

  inline MDistro& darkdistro() { return mdistro; };


  void RunEstimate( const RooAbsData&, const std::string& plot = "" );

protected:
  RooRealProxy x;
  RooRealProxy ped;
  RooRealProxy gain;
  RooRealProxy s0;
  RooRealProxy s1;
  RooRealProxy mean;
  RooRealProxy lambda;
  RooRealProxy alpha;
  RooRealProxy beta;
  RooRealProxy dcfraction;
  RooRealProxy epsilon;

  mutable MDistro mdistro;

  double evaluate() const;

private:
//  ClassDef(SiPMPdf,1);
};

#endif
