#include <vector>
#include <string>
#include "../include/Mesh.h"
#include "../include/Context.h"

using namespace std;


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





