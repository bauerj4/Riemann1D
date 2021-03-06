#include "../include/Context.h"
#include <vector>
#include <string>


/*
  Since we do not require this program to work in parallel,
  a trivial solution to storing mesh information would be to
  store pointers to the heads of the lists which the serial 
  code is perfectly capable of accessing.
*/


// pragma once allows us to define mesh_t in the header
#pragma once
struct mesh_t{
  int NCells;
  double * FirstMeshElement;
  double * FirstDensityElement;
  double * FirstPressureElement;
  double * FirstVelocityElement;
};

int ConstructMesh(mesh_t &mesh, context_t &RiemannContext);
int PrintDataToFile(mesh_t &mesh, context_t &RiemannContext, int snapno);
int PrintMeshData(mesh_t &mesh);
