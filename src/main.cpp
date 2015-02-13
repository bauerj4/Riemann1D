#include "../include/Context.h"
#include "../include/Mesh.h"
#include "../include/Fluxes.h"
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
  vector<vector<double> > fluxes(mesh.NCells, vector<double>(3,0));

  // For time, do:
  // Compute fluxes
  HLLC_FLUX(mesh, fluxes);
  // Update mesh
  // Calculate timestep
  // Update timestep

  // Write output.

  
  PrintDataToFile(mesh, RiemannContext, snapshot_number);

  return 0;
  
}
