#include "SiPMCalib/InvSqCalc/interface/MCFormat.hpp"
#include "UserUtils/PlotUtils/interface/Ratio1DCanvas.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"
#include "UserUtils/Common/interface/ArgumentExtender.hpp"

double
ZProf( const double* vz, const double* param )
{
  const double z = vz[0];
  const double N = param[0];
  return N / ( z*z );
}


int
main( int argc, char** argv )
{

  { // Aperture comparison graphs
    MCManager mc( "data/MC_RComp_Format.txt" );
    TGraphAsymmErrors mc1 = mc.MakeZScanGraph( 1, 1.5, 0 );
    TGraphAsymmErrors mc2 = mc.MakeZScanGraph( 2, 1.5, 0 );
    TGraphAsymmErrors mc3 = mc.MakeZScanGraph( 3, 1.5, 0 );
    TGraphAsymmErrors mc4 = mc.MakeZScanGraph( 4, 1.5, 0 );

    usr::plt::Simple1DCanvas c;
    c.PlotGraph( mc1,
      usr::plt::EntryText("Aperture r = 1mm"),
      usr::plt::TrackY(usr::plt::tracky::both)
      );
    c.PlotGraph( mc2,
      usr::plt::EntryText("Aperture r = 2mm"),
      usr::plt::TrackY(usr::plt::tracky::both)
      );
    c.PlotGraph( mc3,
    usr::plt::EntryText("Aperture r = 3mm"),
      usr::plt::TrackY(usr::plt::tracky::both)
    );
    c.PlotGraph( mc4,
    usr::plt::EntryText("Aperture r = 4mm"),
      usr::plt::TrackY(usr::plt::tracky::both)
    );

    mc1.SetLineColor( kRed );
    mc2.SetLineColor( kBlue );
    mc3.SetLineColor( kGreen );
    mc4.SetLineColor( kViolet );

    c.Pad().SetLogx(1);
    c.Pad().SetLogy(1);
    c.DrawCMSLabel("Simulation", "HCGal" );
    c.DrawLuminosity( "x-offset = 0 [mm]");

    c.Pad().Xaxis().SetTitle("Z distance [mm]");
    c.Pad().Yaxis().SetTitle("Luminosity [A.U.]");

    // c.SaveAsPNG("MC_ApertureCompare.png");
    c.SaveAsPDF("MC_ApertureCompare.pdf");
  }

  {
    MCManager mc("data/MC_R1_Format.txt");
    TGraphAsymmErrors mc0 = mc.MakeZScanGraph( 1, 1.5, 0 );
    TGraphAsymmErrors mc2 = mc.MakeZScanGraph( 1, 1.5, 2 );
    TGraphAsymmErrors mc4 = mc.MakeZScanGraph( 1, 1.5, 4 );
    TGraphAsymmErrors mc8 = mc.MakeZScanGraph( 1, 1.5, 8 );

    usr::plt::Ratio1DCanvas c;
    c.PlotGraph( mc0,
      usr::plt::EntryText("x-offset = 0mm"),
      usr::plt::TrackY(usr::plt::tracky::both)
      );
    c.PlotGraph( mc2,
    usr::plt::EntryText("x-offset = 2 mm"),
      usr::plt::TrackY(usr::plt::tracky::both)
    );
    c.PlotGraph( mc4,
      usr::plt::EntryText("x-offset = 4mm"),
      usr::plt::TrackY(usr::plt::tracky::both)
      );
    c.PlotGraph( mc8,
    usr::plt::EntryText("x-offset = 8mm"),
      usr::plt::TrackY(usr::plt::tracky::both)
    );

    mc0.SetLineColor( kRed );
    mc2.SetLineColor( kViolet );
    mc4.SetLineColor( kBlue );
    mc8.SetLineColor( kGreen );

    c.PlotScale( mc0, mc0 );
    c.PlotScale( mc4, mc0 );
    c.PlotScale( mc8, mc0 );
    c.PlotScale( mc2, mc0 );

    c.TopPad().SetLogx(1);
    c.TopPad().SetLogy(1);
    c.BottomPad().SetLogx(1);
    c.DrawCMSLabel("Simulation", "HCGal" );
    c.DrawLuminosity( "Aperture r = 1 [mm]");

    c.TopPad().Yaxis().SetTitle("Luminosity [A.U.]");
    c.BottomPad().Xaxis().SetTitle("Z distance [mm]");
    c.BottomPad().Yaxis().SetTitle("Lumi / Lumi [0mm]");

    c.SaveAsPNG("MC_OffsetCompare.png");
    // c.SaveAsPDF("MC_OffsetCompare.pdf");

  }

  return 0;
}
