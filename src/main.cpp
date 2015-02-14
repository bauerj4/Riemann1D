#include "../include/Context.h"
#include "../include/Mesh.h"
#include "../include/Fluxes.h"
#include "../include/FiniteVolume.h"
#include <vector>
#include <string>

using namespace std;

int main(int argc, char * argv[])
{

  char * PATH = argv[1];
  
  // Get context.
  
  context_t RiemannContext;
  readContext(PATH, RiemannContext);
  printContext(RiemannContext);

  // Initialize simulation
  
  mesh_t mesh;
  printf("Constructing mesh...\n");
  ConstructMesh(mesh, RiemannContext);
  printf("Mesh constructed.\n");
  
  int snapshot_number = 0;
  PrintDataToFile(mesh, RiemannContext, snapshot_number);

  // Declare flux vector
  vector<vector<double> > fluxes(mesh.NCells - 1, vector<double>(3,0.));
  vector<vector<double> > primatives(mesh.NCells,vector<double> (3,0.));
  vector<vector<double> > conserved(mesh.NCells,vector<double> (3,0.));

  for (int i = 0; i < mesh.NCells; i++)
    {
      primatives[i][0] = mesh.FirstDensityElement[i];
      primatives[i][1] = mesh.FirstVelocityElement[i];
      primatives[i][2] = mesh.FirstPressureElement[i];
      printf("%10.10f, %10.10f, %10.10f", primatives[i][0], primatives[i][1], primatives[i][2]);
    }

  for (int i = 0; i < mesh.NCells; i++)
    {
      conserved[i] = primativeTo1DConservative(primatives[i]);
    }

  double currentTime = 0;
  double SMAX;
  // For time, do:
  while (currentTime < RiemannContext.EVOLVE_TIME)
    {
      // Compute fluxes with TRANSMISSIVE BCs
      HLLC_FLUX(mesh, fluxes, SMAX);
      // Update mesh
      //vector<vector<double> > test(mesh.NCells - 1, vector<double>(3,0.));
      FVUpdate(conserved, fluxes, mesh, RiemannContext, currentTime, SMAX);

      for (int i = 0; i < mesh.NCells; i++)
	{
	  primatives[i] = conservativeTo1DPrimative(conserved[i]);
	  //printf("%10.10f, %10.10f, %10.10f", primatives[i][0], primatives[i][1], primatives[i][2]);
	  mesh.FirstDensityElement[i] = primatives[i][0];
	  mesh.FirstVelocityElement[i] = primatives[i][1];
	  mesh.FirstPressureElement[i] = primatives[i][2];
	}
      SMAX=0;

    }

  // Convert back to primatives
  /*
  for (int i = 0; i < mesh.NCells; i++)
    {
      primatives[i] = conservativeTo1DPrimative(conserved[i]);
      mesh.FirstDensityElement[i] = primatives[i][0];
      mesh.FirstVelocityElement[i] = primatives[i][1];
      mesh.FirstPressureElement[i] = primatives[i][2];
    }
  */


  // Write output.

  
  PrintDataToFile(mesh, RiemannContext, snapshot_number);

  return 0;
  
}
