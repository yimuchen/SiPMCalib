#include "TMath.h"

static const unsigned N = 1584;

double
LinearModel( const double* xx, const double* par )
{
  const double x      = xx[0];
  const double gain   = par[0];
  const double offset = par[1];
  const double N      = par[2];
  return gain*x*N + offset;
}

double
LOModel( const double* xx, const double* par )
{
  const double x      = xx[0];
  const double gain   = par[0];
  const double pde    = par[1];
  const double offset = par[2];
  const double N      = par[3];
  return gain * N * ( 1 - TMath::Exp( -pde * x ) ) + offset;
}

double
NLOModel( const double* xx, const double* par )
{
  const double x      = xx[0];
  const double gain   = par[0];
  const double pde    = par[1];
  const double alpha  = par[2];
  const double beta   = par[3];
  const double offset = par[4];
  const double N      = par[5];

  // Calculating number of fired pixel by LO approximation
  const double LOpar[4]  = {1., pde, 0, N};
  const double Nfired_LO = LOModel( xx, LOpar );

  // Simple linear correction
  const double NLOLinear = ( 1-alpha ) * Nfired_LO
                           + alpha * pde * x * N;

  const double nlo_frac   = pde * x * N / Nfired_LO;
  const double NonLinCorr = ( beta + 1 ) / ( beta + nlo_frac );

  return gain * NLOLinear * NonLinCorr + offset;
}
