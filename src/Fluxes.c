#include <vector>
#include <string>
#include "../include/Mesh.h"
#include "../include/Context.h"

using namespace std;

/*
  For the vector conversion routines,
  use the equation of state:
  
  E = P / (rho (gamma - 1))
*/

int primitiveTo1DConservative(vector<double> &primitive)
{
  if (primitive.size()!=3)
    {
      printf("INCORRECT PRIMITIVE VARIABLE SIZE!\n");
      return 1;
    }
  double gamma = 1.4;
  primitive[1] *= primitive[0];
  // Compute energy from pressure
  double energy = primitive[2] / (primitive[0] * (gamma - 1.));
  primitive[2] = energy * primitive[0];
  return 0;
}

int conservativeTo1DPrimative(vector<double> conservative)
{
  double gamma = 1.4;

  if (conservative.size()!=3)
    {
      printf("INCORRECT CONSERVATIVE VARIABLE SIZE!\n");
      return 1;
    }


  conservative[1] = conservative[1] / conservative[0]; 
  conservative[2] = conservative[2] / conservative[0];
  conservative[2] *= conservative[0] * (gamma - 1.); //from eq of state

  return 0;
}

int HLLC_FLUX(mesh_t &mesh, vector<double> &fluxes)
{
  // From Toro 2009
  double gamma = 1.4; // Given.
  
  
  // Compute quantities in starred region
  
  double p_star, P_L, P_R, rho_bar, rho_L, rho_R, a_bar, a_L, a_R;
  for (int i = 1; i < (mesh.NCells - 1); i++)
    {
      
    }


  return 0;
}





