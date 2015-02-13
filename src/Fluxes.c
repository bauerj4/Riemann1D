#include <vector>
#include <string>
#include <math.h>
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

int HLLC_FLUX(mesh_t &mesh, vector<vector<double> > &fluxes)
{
  // From Toro 2009
  double gamma = 1.4; // Given.
  
  if (fluxes.size() != mesh.NCells)
    {
      printf("INVALID FLUX VECTOR SIZE!\n");
      return 1;
    }
  // Compute quantities in starred region
  
  double u_L, u_R, p_star, P_bar, P_L, P_R, rho_bar, rho_L, rho_R, 
    a_bar, a_L, a_R, q_L, q_R, S_L, S_R, S_star;

  for (int i = 1; i < (mesh.NCells - 1); i++)
    {
      P_L = mesh.FirstPressureElement[i-1];
      P_R = mesh.FirstPressureElement[i];
      P_bar = 0.5 * (P_L + P_R);

      rho_L = mesh.FirstDensityElement[i-1];
      rho_R = mesh.FirstDensityElement[i];
      rho_bar = 0.5 * (rho_L + rho_R);

      // The a's represent the sound speeds.
      // Compute as consequence of EoS
      a_L = pow((gamma * P_L / rho_L),0.5);
      a_R = pow((gamma * P_R / rho_R),0.5);
      a_bar = 0.5 * (a_L + a_R);

      u_L = mesh.FirstVelocityElement[i-1];
      u_R = mesh.FirstVelocityElement[i];

      p_star = P_bar - 0.5 * (u_R - u_L) * a_bar * rho_bar;

      if (p_star < 0)
	{
	  p_star = 0;
	}

      if (p_star <= P_L)
	{
	  q_L = 1;
	}
      else
	{
	  double num = 1 + ((gamma + 1)/(2 * gamma)) * (p_star/P_L - 1);
	  q_L = pow(num,0.5);
	}

      if (p_star <= P_R)
	{
          q_R =1;
	}
      else
	{
          double num = 1 + ((gamma + 1)/(2 * gamma)) * (p_star/P_R - 1);
          q_R =pow(num,0.5);
	}

      S_L = u_L - a_L * q_L;
      S_R = u_R + a_R * q_R;

      double numerator = P_R - P_L + rho_L * u_L * (S_L - u_L) - rho_R * u_R * (S_R - u_R);
      double denominator = rho_L * (S_L - u_L) - rho_R * (S_R - u_R);
      S_star = numerator / denominator;

      printf("The shock speed in the star region is %10.10f\n", S_star);
      
    }


  return 0;
}





