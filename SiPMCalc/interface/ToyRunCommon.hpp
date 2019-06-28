#ifndef SIPMCALIB_SIPMCALC_TOYRUNCOMMON
#define SIPMCALIB_SIPMCALC_TOYRUNCOMMON

#include <string>
#include <UserUtils/Common/interface/ArgumentExtender.hpp>
#include <UserUtils/Common/interface/STLUtils/Filesystem.hpp>

// Defining common initial values:
usr::po::options_description ToyOptions();

// Initial values as determined by argument parsing
double Pedestal( const usr::ArgumentExtender& );
double Gain( const usr::ArgumentExtender& );
double ComNoise( const usr::ArgumentExtender& );
double PixNoise( const usr::ArgumentExtender& );
double Mean( const usr::ArgumentExtender& );
double Lambda( const usr::ArgumentExtender& );
double DCFrac( const usr::ArgumentExtender& );
double Alpha( const usr::ArgumentExtender& );
double Beta( const usr::ArgumentExtender& );
double Xmin( const usr::ArgumentExtender& );
double Xmax( const usr::ArgumentExtender& );

usr::fs::path filename( const usr::ArgumentExtender& );
usr::fs::path pullfilename( const usr::ArgumentExtender& );
usr::fs::path pullplotfilename( const usr::ArgumentExtender&,
                                const std::string& );


#endif
