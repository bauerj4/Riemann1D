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
  vector<double> testVector(3,2.0);
  // Primative should invert conservative
  testVector = primativeTo1DConservative(testVector);
  printf("The conserved test components are: [%10.10f, %10.10f, %10.10f]\n", testVector[0],testVector[1], testVector[2]);
  //
  testVector = conservativeTo1DPrimative(testVector);
  printf("The primative components are: [%10.10f, %10.10f, %10.10f]\n", testVector[0],testVector[1], testVector[2]);
  
  
  // Get context.
  
  context_t RiemannContext;
  readContext(PATH, RiemannContext);
  printContext(RiemannContext);
  
  // Initialize simulation
  mesh_t mesh;
  printf("Constructing mesh...\n");
  ConstructMesh(mesh, RiemannContext);
  printf("Mesh constructed.\n");
  
  // Add more ghost cells for higher order
  /*
    if (RiemannContext.SOLUTION_METHOD == "HLLC_FLUX_SUPERBEE" || RiemannContext.SOLUTION_METHOD == "HLLC_FLUX_SUPERBEE")
    {
    double positions[mesh.NCells + 2];
    double pressures[mesh.NCells + 2];
    double velocities[mesh.NCells + 2];
    double denisities[mesh.NCells + 2];
    mesh.NCells = mesh.NCells + 2;
    for (int i = 0; i < mesh.NCells; i++)/ 	
    {
    if i == mesh.NCells - 1)
    {
    positions[i] = positions[i - 1];
    pressures[i] = pressures[i - 1];
    velocities[i] = velocities[i - 1];
    }
    }
    
    }
   */
  
  int snapshot_number = 0;
  PrintDataToFile(mesh, RiemannContext, snapshot_number);
  //return 1;
  // Declare flux vector
  vector<vector<double> > fluxes(mesh.NCells - 3, vector<double>(3,0.));
  
  if(RiemannContext.SOLUTION_METHOD == "HLLC")
    {
      fluxes.push_back(vector<double>(3,0.));
      fluxes.push_back(vector<double>(3,0.));
    }
  
  vector<vector<double> > primatives(mesh.NCells,vector<double> (3,0.));
  vector<vector<double> > conserved(mesh.NCells,vector<double> (3,0.));
  
  for (int i = 0; i < mesh.NCells; i++)
    {
      primatives[i][0] = mesh.FirstDensityElement[i];
      primatives[i][1] = mesh.FirstVelocityElement[i];
      primatives[i][2] = mesh.FirstPressureElement[i];
      //printf("%10.10f, %10.10f, %10.10f", primatives[i][0], primatives[i][1], primatives[i][2]);
    }
  
  for (int i = 0; i < mesh.NCells; i++)
    {
      conserved[i] = primativeTo1DConservative(primatives[i]);
    }
  
  double currentTime = 0;
  double SMAX;
  int iteration;
  // For time, do:
  while (currentTime < RiemannContext.EVOLVE_TIME)
    {
      // Compute fluxes with TRANSMISSIVE BCs
      if (RiemannContext.SOLUTION_METHOD == "HLLC")
 	{
 	  HLLC_FLUX(mesh, fluxes, SMAX);
 	}
      
      if (RiemannContext.SOLUTION_METHOD=="HLLC_FLUX_SUPERBEE")
 	{
 	  int success = 1;
 	  success = HLLC_FLUX_SUPERBEE(mesh, fluxes, RiemannContext, SMAX);
 	  printf("Flux calculation finished.\n");
 	  if (success != 0)
 	    {
 	      printf("RIEMANN SOLVER RETURNED WITH ERROR.\n");	      
PrintDataToFile(mesh, RiemannContext, snapshot_number);
 	      break;
 	    }
 	  printf("FLUX RETURN CODE %d\n",success);
 	}
       // Update mesh
      //vector<vector<double> > test(mesh.NCells - 1, vector<double>(3,0.));
       printf("Updating Mesh.\n");
      
      int success = FVUpdate(conserved, fluxes, mesh, RiemannContext, currentTime, SMAX, iteration);
if (success != 0)
	{
	  printf("RIEMANN SOLVER RETURNED WITH ERROR.\n");
	  PrintDataToFile(mesh, RiemannContext, snapshot_number);
	  break;
	}


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
