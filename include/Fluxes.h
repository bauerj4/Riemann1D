#include <vector>
#include <string>
#include "Mesh.h"
#include "Context.h"

/*
  This will be the complicated part.  Here, we implement the flux evolution components.
  Having constructed the mesh, these functions will take and evolve the primitive variables
  on the mesh and compute an array of N+1 fluxes which can be used by any approximate
  Riemann solver scheme in conjunction with a finite volume method to update the conserved
  quantities.  These can be recast to primitive quantities fairly simply.  Formulas will 
  be taken from Toro 2009 in the case of HLLC flux.
*/

using namespace std;

// The fluxes need not be printed or converted, so 
// there is no harm in using a vector.

vector<double> primativeTo1DConservative(vector<double> primative);
vector<double> conservativeTo1DPrimative(vector<double> conservative);
int HLLC_FLUX(mesh_t &mesh, vector<vector<double> > &fluxes, double &smax); 
int HLLC_FLUX_SUPERBEE(mesh_t &mesh, vector<vector<double> > &fluxes, context_t &RiemannContext,double &smax);
double SUPERBEE(double r, double c);
double MINBEE(double r, double c);
