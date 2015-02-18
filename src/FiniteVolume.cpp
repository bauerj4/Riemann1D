#include "../include/FiniteVolume.h"
#include "../include/Fluxes.h"
#include "../include/Mesh.h"
#include <string>
#include <vector>
#include <math.h>

using namespace std;

int FVUpdate(vector<vector<double> > &conserved, vector<vector<double> > &fluxes,
	     mesh_t &mesh, context_t &RiemannContext, double &current_time, double smax)
{
  printf("smax = %10.10f\n",smax);
  double dx = fabs(RiemannContext.X1 - RiemannContext.X0)/RiemannContext.MESH_RESOLUTION;
  printf("dx = %10.10f\n", dx);
  double dt = dx/smax * RiemannContext.CFL_NUMBER;
  printf("dt = %10.10f\n", dt);

  if(RiemannContext.SOLUTION_METHOD == "HLLC_FLUX_SUPERBEE" ||
     RiemannContext.SOLUTION_METHOD == "HLLC_FLUX_MINMOD")
    {
      printf("Updating conserved quantities \n");
      for(int i = 2; i < mesh.NCells - 2; i++)
        {
	  printf("Updating conserved %d.\n",i+1);

	  if ((fluxes[i-1][0] - fluxes[i-2][0]) != (fluxes[i-1][0] - fluxes[i-2][0]))
	    {
	      printf("ERROR: FLUXES ARE NAN.\n");
	      return 1;
	    }
          conserved[i][0] -= dt/dx * (fluxes[i-1][0] - fluxes[i-2][0]);
          conserved[i][1] -= dt/dx * (fluxes[i-1][1] - fluxes[i-2][1]);
          conserved[i][2] -= dt/dx * (fluxes[i-1][2] - fluxes[i-2][2]);
        }
      printf("Conserved quantities updated.\n");
    }

  else
    {
      for(int i = 1; i < mesh.NCells - 1; i++)
	{
	  conserved[i][0] -= dt/dx * (fluxes[i][0] - fluxes[i-1][0]);
	  conserved[i][1] -= dt/dx * (fluxes[i][1] - fluxes[i-1][1]);
	  conserved[i][2] -= dt/dx * (fluxes[i][2] - fluxes[i-1][2]);
	}
    }
  
  //SPECIFY BCS FOR GHOST CELLS
  if(RiemannContext.SOLUTION_METHOD == "HLLC_FLUX_SUPERBEE" || 
     RiemannContext.SOLUTION_METHOD == "HLLC_FLUX_MINMOD")
    {
      conserved[mesh.NCells - 1][0] = conserved[mesh.NCells- 3][0];
      conserved[mesh.NCells - 1][1] = conserved[mesh.NCells- 3][1];
      conserved[mesh.NCells - 1][2] = conserved[mesh.NCells- 3][2];
      conserved[mesh.NCells - 2][0] = conserved[mesh.NCells- 3][0];
      conserved[mesh.NCells - 2][1] = conserved[mesh.NCells- 3][1];
      conserved[mesh.NCells - 2][2] = conserved[mesh.NCells- 3][2];
    }

  else
    {
      conserved[mesh.NCells - 1][0] = conserved[mesh.NCells- 2][0];
      conserved[mesh.NCells - 1][1] = conserved[mesh.NCells- 2][1];
      conserved[mesh.NCells - 1][2] = conserved[mesh.NCells- 2][2];
    }

  current_time += dt;
  printf("t + dt = %10.10f + %10.10f\n", current_time, dt);
  return 0;
}
