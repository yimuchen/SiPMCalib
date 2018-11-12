#include "SiPMCalib/SiPMCalc/interface/SiPMPdf_Functor.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

#include "RooFunctorBinding.h"
#include <iostream>

int main(int argc, char const *argv[])
{
  SiPMPdf_Functor func;
  double param[9] = {300, 350, 200 , 50 , 10, 1.5,0.01,0.001, 50 };
  double x[1000], y[1000];

  for( int i = 0 ; i < 1000 ; ++i ){
    x[i] = 300 + i * (1500. - 300.) / 1000.;
    param[0] = x[i];
    y[i] = 1e30*func(param);
    std::cout << x[i] << " " << y[i] << std::endl;
  }

  TGraph graph(1000,x,y);

  usr::plt::Simple1DCanvas c(500,500);
  c.PlotGraph( graph, usr::plt::TrackY(usr::plt::TrackY::both) );

  c.SaveAsPNG("test.png");
  c.SetLogy(kTRUE);

  return 0;
}
