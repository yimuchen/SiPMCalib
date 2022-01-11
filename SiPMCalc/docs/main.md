@defgroup SiPMCalc SiPMCalc

@brief SiPM analysis related computation routines

This is a collection of functions and classes related to the mathematic
representation of functions used for SiPM property analysis. For classes, this
includes the implementation of the related mathematical functions as RooFit
compatible PDF classes, as well as the implementation of "fitting routines as a
class" so that configurations for fitting can be tested and re-used for bulk
analysis where we are expecting the SiPM to have similar properties. For files,
this would mainly be for example fitting routines for the purpose of testing the
"fitting classes", and presenting the results of SiPM analysis while the bulk
calibration routine is still in development.

As a rule of thumb for the fitter routine classes, the helper string within the
code should be short and concise (spanning at most 2 lines on the command line).
The detailed documentation should be provided either in the detailed
documentation of the function used to generate the program options, or in the
detailed documentation of the script file.
