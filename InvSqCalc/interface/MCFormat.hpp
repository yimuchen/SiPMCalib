#ifndef SIPMCALIB_INVSQCALC_MCFORMAT_HPP
#define SIPMCALIB_INVSQCALC_MCFORMAT_HPP

#include "UserUtils/MathUtils/interface/Measurement/Measurement.hpp"
#include "UserUtils/Common/interface/STLUtils/Filesystem.hpp"

#include <vector>
#include "TGraphAsymmErrors.h"

class Model {
public:
  const double r ; // Aperture radius
  const double l ; // SiPM edge length
  const double o ; // offset
  const double z ; // z distance
  unsigned total ; // Total number of toys
  unsigned passed ;  // Total number of

  Model( const double, const double, const double , const double );
  ~Model();

  double MaxTheta() const ;
  double Sangle() const ;
  void Run( const unsigned );

  usr::Measurement Lumi() const ;
};


class MCManager {
public:
  std::vector<Model> _modellist;

  MCManager();
  MCManager( const usr::fs::path& ) ; // Initializing from file.
  ~MCManager();

  Model& GetModel(
    const double r, const double l, const double o, const double z );

  TGraphAsymmErrors MakeZScanGraph( const double r, const double l, const double o ) const ;
  TGraphAsymmErrors MakeHScanGraph( const double r, const double l, const double z ) const ;

  void SaveToTXT( const usr::fs::path& );
};

#endif