#ifndef NON_LINEAR_MODEL_HPP
#define NON_LINEAR_MODEL_HPP
#include "UserUtils/MathUtils/interface/Measurement/Measurement.hpp"

extern double LinearModel( const double*, const double* );
extern double LOModel( const double*, const double* );
extern double NLOModel( const double*, const double* );
extern double BiasModel( const double*, const double * );


// Filter wheel filtering values.
extern usr::Measurement A( const unsigned idx );
extern usr::Measurement B( const unsigned idx );


#endif
