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

vector<double> primativeTo1DConservative(vector<double> &primative)
{
  if (primative.size()!=3)
    {
      printf("INCORRECT PRIMATIVE VARIABLE SIZE!\n");
      return primative;
    }
  double gamma = 1.4;
  double density = primative[0];
  double velocity = primative[1];
  double pressure = primative[2];
  double conservative_2 = density * velocity;
  //primative[1] *= primative[0];
  // Compute energy from pressure
  double int_energy = pressure / (density * (gamma - 1.));
  double kinetic_energy = 0.5 * velocity * velocity;

  double conservative_3 = density * (int_energy + kinetic_energy);
  vector<double> conservative(3,0.0);
  conservative[0] = density;
  conservative[1] = conservative_2;
  conservative[2] = conservative_3;
  //primative[2] = energy * primative[0] + pow(primative[1]/primative[0],2.0);
  return conservative;
}

vector<double> conservativeTo1DPrimative(vector<double> &conservative)
{
  double gamma = 1.4;

  if (conservative.size()!=3)
    {
      printf("INCORRECT CONSERVATIVE VARIABLE SIZE!\n");
      return conservative;
    }


  double density = conservative[0]; 
  double velocity = conservative[1] / density;
  double spec_energy = conservative[2] / density;
  //from eq of state
  double int_energy = spec_energy - 0.5 * velocity * velocity;
  double pressure = int_energy * density*(gamma -1.);

  vector<double> primative(3,0.0);
  primative[0] = density;
  primative[1] = velocity;
  primative[2] = pressure;

  return primative;
}

int HLLC_FLUX(mesh_t &mesh, vector<vector<double> > &fluxes, double &smax)
{
  // From Toro 2009
  double gamma = 1.4; // Given.
  
  if (fluxes.size() != (mesh.NCells - 1)) // Should confirm that this is the correct number
    {
      printf("INVALID FLUX VECTOR SIZE!\n");
      return 1;
    }
  // Compute quantities in starred region
  
  double u_L, u_R, p_star, P_bar, P_L, P_R, rho_bar, rho_L, rho_R, tempSmax,
    a_bar, a_L, a_R, q_L, q_R, S_L, S_R, S_star, E_L, E_R, prefactor_L, prefactor_R;

  vector<double> Conserved_R(3,0);
  vector<double> Conserved_L(3,0);
  vector<double> Conserved_R_Star(3,0);
  vector<double> Conserved_L_Star(3,0);
  vector<double> Flux_R(3,0);
  vector<double> Flux_L(3,0);
  vector<double> Flux_R_Star(3,0);
  vector<double> Flux_L_Star(3,0);
  
  tempSmax = 0.;
  for (int i = 1; i < (mesh.NCells ); i++)
    {
      //printf("ITERATION %d\n",i);
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
      //printf("SOUND SPEEDS: %10.10f, %10.10f, %10.10f\n", a_L,a_R,a_bar);

      u_L = mesh.FirstVelocityElement[i-1];
      u_R = mesh.FirstVelocityElement[i];

      p_star = P_bar - 0.5 * (u_R - u_L) * a_bar * rho_bar;
      //printf("P STAR IS %10.10f\n", p_star);
      if (p_star < 0)
	{
	  //printf("P star is 0.\n");
	  p_star = 0.;
	}

      if (p_star <= P_L)
	{
	  q_L = 1.;
	}
      else
	{
	  double num = 1 + ((gamma + 1)/(2 * gamma)) * ((p_star/P_L) - 1);
	  q_L = pow(num,0.5);
	}

      if (p_star <= P_R)
	{
          q_R =1.;
	}
      else
	{
          double num = 1. + ((gamma + 1.)/(2. * gamma)) * ((p_star/P_R) - 1.);
          q_R =pow(num,0.5);
	}

      S_L = u_L - a_L * q_L;
      S_R = u_R + a_R * q_R;

      // MAKE SURE THIS IS RIGHT
      if (fabs(u_L) + a_L > tempSmax)
	{
	  tempSmax = fabs(u_L) + a_L;
	}
      if (fabs(u_R) + a_R > tempSmax)
	{
          tempSmax = fabs(u_R) + a_R;
        }



      double numerator = P_R - P_L + rho_L * u_L * (S_L - u_L) - rho_R * u_R * (S_R - u_R);
      double denominator = rho_L * (S_L - u_L) - rho_R * (S_R - u_R);
      S_star = numerator / denominator;
      //printf("THE MAX SIGNAL SPEEDS ARE: %10.10f, %10.10f, %10.10f\n", S_L,S_R, S_star);

      //printf("S_star is %10.10f\n",S_star);

      //printf("The shock speed in the star region is %10.10f\n", S_star);

      // Construct conserved L/R and L/R star vectors from previously defined functions

      vector<double> primative_L(3);
      double primative_L_arr[3] = {rho_L, u_L, P_L};
      primative_L.assign(&primative_L_arr[0], &primative_L_arr[0] + 3);
      //printf("%10.10f, %10.10f,%10.10f\n", primative_L[0], primative_L[1], primative_L[2]);
      //printf("The size of the primative variables is %d\n",(int)primative_L.size());
      vector<double> primative_R(3);
      double primative_R_arr[3] = {rho_R, u_R, P_R};
      primative_R.assign(&primative_R_arr[0], &primative_R_arr[0] + 3);


      Conserved_R = primativeTo1DConservative(primative_R);
      //printf("%10.10f, %10.10f,%10.10f\n", Conserved_R[0], Conserved_R[1], Conserved_R[2]); 

      Conserved_L = primativeTo1DConservative(primative_L);

      prefactor_L = rho_L * (S_L - u_L) / (S_L - S_star);
      prefactor_R = rho_R * (S_R - u_R) / (S_R - S_star);

      E_L = rho_L * (primative_L[2] / (primative_L[0]*(gamma - 1)) + 0.5 * pow(u_L,2.0));
      E_R = rho_R *(primative_R[2] / (primative_R[0]*(gamma - 1)) + 0.5 * pow(u_R,2.0));
      

      Conserved_L_Star[0] = prefactor_L;
      Conserved_L_Star[1] = prefactor_L * S_star;
      Conserved_L_Star[2] = prefactor_L * ((E_L / rho_L) + (S_star - u_L) * (S_star + (P_L)/(rho_L*(S_L - u_L))));

      Conserved_R_Star[0] = prefactor_R;
      Conserved_R_Star[1] = prefactor_R * S_star;
      Conserved_R_Star[2] = prefactor_R * ((E_R / rho_R) + (S_star - u_R) * (S_star + (P_R)/(rho_R*(S_R - u_R))));

      Flux_R[0] = primative_R[0] * primative_R[1];
      Flux_R[1] = primative_R[0] * primative_R[1] * primative_R[1] + P_R;
      Flux_R[2] = primative_R[1] * (E_R + P_R);

      Flux_L[0] = primative_L[0] * primative_L[1];
      Flux_L[1] = primative_L[0] * primative_L[1] * primative_L[1] + P_L;
      Flux_L[2] = primative_L[1] * (E_L + P_L);

      Flux_R_Star[0] = Flux_R[0] + S_R * (Conserved_R_Star[0] - Conserved_R[0]);
      Flux_R_Star[1] = Flux_R[1] + S_R * (Conserved_R_Star[1] - Conserved_R[1]);
      Flux_R_Star[2] = Flux_R[2] + S_R * (Conserved_R_Star[2] - Conserved_R[2]);

      Flux_L_Star[0] = Flux_L[0] + S_L * (Conserved_L_Star[0] - Conserved_L[0]);
      Flux_L_Star[1] = Flux_L[1] + S_L * (Conserved_L_Star[1] - Conserved_L[1]);
      Flux_L_Star[2] = Flux_L[2] + S_L * (Conserved_L_Star[2] - Conserved_L[2]);

      // Construct the HLLC flux
      //printf("Constructing fluxes...\n");
      int debugCount = 0;
      if (S_L >= 0)
	{
	  fluxes[i-1] = Flux_L;
	  //printf("Condition 1 tripped.\n");
	  debugCount++;
	}

      if (S_star >= 0 && S_L <= 0)
	{
	  fluxes[i-1] = Flux_L_Star;
          //printf("Condition 2 tripped.\n");
	  debugCount++;
	}

      if (S_star <= 0 && S_R >= 0)
	{
	  fluxes[i-1] = Flux_R_Star;
          //printf("Condition 3 tripped.\n");
	  debugCount++;
	}

      if (S_R <= 0)
	{
	  fluxes[i-1] = Flux_R;
          //printf("Condition 4 tripped.\n");
	  debugCount++;
	}
      //printf("The debug count is %d\n", debugCount);
      
    }
  smax = tempSmax;
  return 0;
}





