#include "../include/Fluxes.h"
#include "../include/Context.h"
#include "../include/Mesh.h"
#include <string>
#include <vector>

/*
  Here we implement the finite volume method for updating the computed fluxes.
  This may, in principle, be done with any method.  This will not effect the form
  of finite volume with the methods we've chosen. 
*/

using namespace std;

int FVUpdate(vector<vector<double> > &conserved, vector<vector<double> > &fluxes,
	     mesh_t &mesh, context_t &RiemannContext, double &current_time, double smax, int &i);
