#ifndef SIPMCALIB_INVSQCALC_LEDFORMAT_HPP
#define SIPMCALIB_INVSQCALC_LEDFORMAT_HPP

#include "UserUtils/Common/interface/STLUtils/Filesystem.hpp"

#include "TGraph.h"
#include "TH2D.h"
#include <vector>

class LEDPoint
{
public:
  double x;
  double y;
  double z;
  double lumi;
  double lumierr;
};

class LEDManager
{
public:
  std::vector<LEDPoint> _pointlist;

  LEDManager( const usr::fs::path file );

  double Xmin() const ;
  double Xmax() const ;
  double Ymin() const ;
  double Ymax() const ;
  double Zmin() const ;
  double Zmax() const ;
  double LumiMax() const ;
  double LumiMin() const ;
  void      FindCenter( const double z, double& x, double& y ) const;
  TH2D*     MakeHScanGraph( const double z )  const;
  TGraph*   MakeZScanGraph( const double x, const double y ) const;
  TGraph*   MakeXScanGraph( const double y, const double z ) const;
};

#endif
