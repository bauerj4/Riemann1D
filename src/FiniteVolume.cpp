#include "../include/FiniteVolume.h"
#include "../include/Fluxes.h"
#include "../include/Mesh.h"
#include <string>
#include <vector>
#include <math.h>

int FVUpdate(vector<vector<double> > &conserved, vector<vector<double> > &fluxes,
	     mesh_t &mesh, context_t RiemannContext, double &current_time, double smax)
{
  double dx = fabs(RiemannContext.X1 - RiemannContext.X0)/RiemannContext.MESH_RESOLUTION;
  double dt = dx/smax * RiemannContext.CFL_NUMBER;

  current_time += dt;
  return 0;
}
