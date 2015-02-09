#include "../include/Context.h"
#include "../include/Mesh.h"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
#include <cstdlib>

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
  mesh.NCells = RiemannContext.MESH_RESOLUTION + 2; // For boundaries

  double cellWidth = (rightBound - leftBound) / ((double) (RiemannContext.MESH_RESOLUTION));

  // Construct mesh and find the index best place on the Mesh to split 
  double  density[mesh.NCells];
  double  velocity[mesh.NCells];
  double  pressure[mesh.NCells];
  double  positions[mesh.NCells];
  double difference = rightBound - leftBound + 1.; // Initial difference can't be physical
  int indexOfDiscontinuity = mesh.NCells;

  for (int i = 1; i < mesh.NCells - 1; i++)
    {
      positions[i] = ((double) i - 1) * cellWidth;
      if (fabs(positions[i] - RiemannContext.INITIAL_DISCONTINUITY) < difference)
	{
	  indexOfDiscontinuity = i;
	}
      //printf("The difference is %10.10f\n", difference);

      difference = fabs(positions[i] - RiemannContext.INITIAL_DISCONTINUITY);
    }
  positions[0] = positions[1] - cellWidth;
  positions[mesh.NCells - 1] = positions[mesh.NCells - 2] + cellWidth;
  
  //printf("The discontinuity position is %d\n",indexOfDiscontinuity);
  for (int i = 0; i < mesh.NCells; i++)
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
  
  mesh.FirstMeshElement = (double*)calloc(mesh.NCells, sizeof(double));//positions;
  mesh.FirstDensityElement = (double*)calloc(mesh.NCells, sizeof(double));//density;
  mesh.FirstPressureElement = (double*)calloc(mesh.NCells, sizeof(double));//pressure;
  mesh.FirstVelocityElement = (double*)calloc(mesh.NCells, sizeof(double));//velocity;

  
  for (int k = 0; k < mesh.NCells; k++)
    {
      /* printf("%10.10f, %10.10f, %10.10f, %10.10f\n", mesh.FirstMeshElement[k], mesh.FirstDensityElement[k], 
	     mesh.FirstPressureElement[k], mesh.FirstVelocityElement[k]);
      */
      mesh.FirstMeshElement[k] = positions[k];
      mesh.FirstDensityElement[k] = density[k];
      mesh.FirstPressureElement[k] = pressure[k];
      mesh.FirstVelocityElement[k] = velocity[k];
    }
  
  return 0;
}


int PrintMeshData(mesh_t &mesh)
{
  for (int k = 0; k < mesh.NCells; k++)
    {
      printf("%10.10f, %10.10f, %10.10f, %10.10f\n", mesh.FirstMeshElement[k], mesh.FirstDensityElement[k],
             mesh.FirstPressureElement[k], mesh.FirstVelocityElement[k]);
    }
  return 0;
}

int PrintDataToFile(mesh_t &mesh, context_t &RiemannContext, int snapno)
{
  stringstream cppPath;
  cppPath <<  RiemannContext.SNAPSHOT_PATH << snapno;
  const char * PATH = (cppPath.str()).c_str();
  FILE * f = fopen(PATH, "w");
  
  for (int i = 0; i < mesh.NCells; i++)
    {
      fprintf(f,"%10.5f %10.5f %10.5f %10.5f\n", mesh.FirstMeshElement[i], 
	      mesh.FirstDensityElement[i], mesh.FirstPressureElement[i], mesh.FirstVelocityElement[i]);
    }
  fclose(f);
  
  return 0;
}
