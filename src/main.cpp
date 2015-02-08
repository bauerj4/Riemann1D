#include "../include/Context.h"
#include "../include/Mesh.h"
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
  ConstructMesh(mesh, RiemannContext);
  int snapshot_number;
  

  // For time, do:
  // Compute fluxes
  // Update mesh
  // Calculate timestep
  // Update timestep

  // Write output.

  
  PrintDataToFile(mesh, RiemannContext, snapshot_number);

  return 0;
  
}
