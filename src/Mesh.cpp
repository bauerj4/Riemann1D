#include "../include/Context.h"
#include <vector>
#include <string>


// The mesh is constructed from context.
int ConstructMesh(mesh_t &mesh, context_t &RiemannContext)
{
  double leftBound = RiemannContext.X0;
  double rightBound = RiemannContext.X1;
  double discontinuityPosition = RiemannContext.INITIAL_DISCONTINUITY;
  double leftDensity = RiemannContext.RHO_L;
  double rightDensity = RiemannContext.RHO_R;
  double leftPressure = RiemannContext.P_L;
  double rightPressure = RiemannContext.P_R;
  double leftVelocity = RiemannContext.U_L;
  double rightVelocity = RiemannContext.U_R;

  // construct bin values
  mesh.NCells = RiemannContext.MESH_RESOLUTION;
  double cellWidth = (rightBound - leftBound) / ((double) (RiemannContext.MESH_RESOLUTION));

  // Construct mesh and find the index best place on the Mesh to split 
  double * density[mesh.NCells];
  double * velocity[mesh.NCells];
  double * pressure[mesh.NCells];
  double * positions[mesh.NCells];
  double difference = rightBound - leftBound + 1.; // Initial difference can't be physical
  int indexOfDiscontinuity = mesh.NCells;

  for (int i = 0; i < mesh.NCells; i++)
    {
      positions[i] = (double) i * cellWidth;
      if ((positions[i] - RiemannContext.INITIAL_DISCONTINUITY) < difference)
	{
	  indexOfDiscontinuity = i;
	}
      difference = positions[i] - RiemannContext.INITIAL_DISCONTINUITY;
    }
  
  for (int i = 0; i < mesh.NCells)
    {
      if (i < indexOfDiscontinuity)
	{
	  density[i] = leftDensity;
	  pressure[i] = leftPressure;
	  velocity[i] = leftVelocity;
	}

      else
	{
	  density[i] = rightDensity;
	  pressure[i] = rightPressure;
	  velocity[i] = rightPressure;
	}
    }
  
  mesh.FirstMeshElement = &positions;
  mesh.FirstDensityElement = &density;
  mesh.FirstPressureElement = &pressure;
  mesh.FirstVelocityElement = &velocity;

  return 0;
}
