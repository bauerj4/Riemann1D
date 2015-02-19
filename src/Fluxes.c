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
double SUPERBEE(double r, double c)
{
  if (r <= 0)
    {
      return 1.0;
    }
  else if(r > 0 && r <= 0.5)
    {
      return (1. - 2.*(1. - c) * r);
    }
  else if (r > 0.5 && r <= 1.0)
    {
      return c;
    }
  else if (r > 1.0 && r <= 2.0)
    {
      return (1. - (1. - c)*r);
    }
  else if (r > 2.0)
    {
      return ((2. * c) - 1.);
    }
}

double MINBEE(double r, double c)
{
  if (r <= 0.)
    {
      return 1.0;
    }
  else if(r > 0. && r <= 1.) 
    {
      return (1 - (1 - c)*r);
    }
  else if(r > 1)
    {
      return c;
    }
}


vector<double> primativeTo1DConservative(vector<double> primative)
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

vector<double> conservativeTo1DPrimative(vector<double> conservative)
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
  vector<double> Flux_R(3, 0.0);
  vector<double> Flux_L(3, 0.0);
  vector<double> Flux_R_Star(3,0.0);
  vector<double> Flux_L_Star(3,0.0);
  
  tempSmax = 0.;
  //print("COMPUTING FUNDAMENTALS...\n")
  for (int i = 1; i < (mesh.NCells ); i++)
    {
      //printf("ITERATION %d\n",i);
      P_L = mesh.FirstPressureElement[i-1];
      P_R = mesh.FirstPressureElement[i];

      if (P_L < 0)
	{
	  P_L = 0;
	  printf("WARNING: P_L < 0\n");
	}

      if (P_R < 0)
	{
          P_R =0;
          printf("WARNING: P_R < 0\n");
	}

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
      printf("P STAR IS %10.10f\n", p_star);
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
	  double num = 1 + (((gamma + 1)/(2 * gamma)) * ((p_star/P_L) - 1));
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



      double numerator = P_R - P_L + (rho_L * u_L * (S_L - u_L)) - (rho_R * u_R * (S_R - u_R));
      double denominator = (rho_L * (S_L - u_L)) - (rho_R * (S_R - u_R));
      S_star = numerator / denominator;
      //printf("THE MAX SIGNAL SPEEDS ARE: %10.10f, %10.10f, %10.10f\n", S_L,S_R, S_star);

      printf("S_star is %10.10f\n",S_star);

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

      E_L = rho_L *(P_L / (rho_L*(gamma - 1)) + 0.5 * pow(u_L,2.0));
      E_R = rho_R *(P_R / (rho_R*(gamma - 1)) + 0.5 * pow(u_R,2.0));
      

      Conserved_L_Star[0] = prefactor_L;
      Conserved_L_Star[1] = prefactor_L * S_star;
      Conserved_L_Star[2] = prefactor_L * ((E_L / rho_L) + (S_star - u_L) * (S_star + (P_L)/(rho_L*(S_L - u_L))));

      Conserved_R_Star[0] = prefactor_R;
      Conserved_R_Star[1] = prefactor_R * S_star;
      Conserved_R_Star[2] = prefactor_R * ((E_R / rho_R) + (S_star - u_R) * (S_star + (P_R)/(rho_R*(S_R - u_R))));

      Flux_R[0] = rho_R * u_R;
      Flux_R[1] = rho_R * u_R * u_R + P_R;
      Flux_R[2] = u_R * (E_R + P_R);

      Flux_L[0] = rho_L * u_L;
      Flux_L[1] = rho_L * u_L * u_L + P_L;
      Flux_L[2] = u_L * (E_L + P_L);

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
	  printf("Condition 1 tripped.\n");
	  debugCount++;
	}

      else if (S_star >= 0 && S_L <= 0)
	{
	  fluxes[i-1] = Flux_L_Star;
          printf("Condition 2 tripped.\n");
	  debugCount++;
	}

      else if (S_star <= 0 && S_R >= 0)
	{
	  fluxes[i-1] = Flux_R_Star;
          printf("Condition 3 tripped.\n");
	  debugCount++;
	}

      else if (S_R <= 0)
	{
	  fluxes[i-1] = Flux_R;
          printf("Condition 4 tripped.\n");
	  debugCount++;
	}
      //printf("The debug count is %d\n", debugCount);
      printf("flux = [%10.10f, %10.10f, %10.10f]\n", fluxes[i-1][0], fluxes[i-1][1], fluxes[i-1][2]);
  
    }

  //  printf("flux = [%10.10f, %10.10f, %10.10f]", fluxes[i][0], fluxes[i][1], fluxes[i][2]);
  smax = tempSmax;
  return 0;
}




int HLLC_FLUX_SUPERBEE(mesh_t &mesh, vector<vector<double> > &fluxes, context_t &RiemannContext, double &smax)
{
  // From Toro 2009
  double gamma = 1.4; // Given.
  
  if (fluxes.size() != (mesh.NCells - 3)) // Should confirm that this is the correct number
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
  vector<double> SignalSpeeds(mesh.NCells - 3, 0.0);
  vector<double> ck(mesh.NCells - 2, 0.0);
  vector<double> rk(mesh.NCells - 2, 0.0);
  vector<vector<double> > FRs(mesh.NCells - 3, vector<double>(3,0.0));
  vector<vector<double> > FLs(mesh.NCells - 3,vector<double>(3,0.0));
  vector<vector<double> > FRStars(mesh.NCells - 3,vector<double>(3,0.0));
  vector<vector<double> > FLStars(mesh.NCells - 3,vector<double>(3,0.0));



  tempSmax = 0.;
  printf("COMPUTING FUNDAMENTALS...\n");

  //SignalSpeeds[mesh.NCells - 2] = 1.;
  //SignalSpeeds[0] = -1.;
  //SignalSpeeds[1] = -1.;

  for (int i = 2; i < (mesh.NCells - 1); i++)
    {
      //printf("ITERATION %d\n",i);
      P_L = mesh.FirstPressureElement[i-1];
      P_R = mesh.FirstPressureElement[i];

      if (P_L < 0)
        {
          P_L = 0;
          printf("WARNING: P_L < 0\n");
        }

      if (P_R < 0)
        {
          P_R =0;
          printf("WARNING: P_R < 0\n");
	}


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

      if (a_L != a_L)
        {
          printf("A_L IS NAN\n");
          return 1;
        }

      if (a_R != a_R)
        {
          printf("A_R IS NAN AND P_R = %10.10f\n", P_R);
          return 1;
        }



      u_L = mesh.FirstVelocityElement[i-1];
      u_R = mesh.FirstVelocityElement[i];

      p_star = P_bar - (0.5 * (u_R - u_L) * a_bar * rho_bar);
      //printf("P STAR IS %10.10f\n", p_star);
      if (p_star < 0)
	{
	  //printf("P star is 0.\n");
	  p_star = 0.;
	}

      if (p_star != p_star)
	{
          printf("P_STAR IS NAN");
          return 1;
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


      if (q_L != q_L)
        {
          printf("Q_L IS NAN\n");
          return 1;
        }

      if (q_R != q_R)
	{
          printf("Q_R IS NAN\n");
          return 1;
        }


      S_L = u_L - a_L * q_L;
      S_R = u_R + a_R * q_R;
      SignalSpeeds[i-2] = S_R;

      // MAKE SURE THIS IS RIGHT
      if (fabs(u_L) + a_L > tempSmax)
	{
	  tempSmax = fabs(u_L) + a_L;
	}
      if (fabs(u_R) + a_R > tempSmax)
	{
          tempSmax = fabs(u_R) + a_R;
        }


      if(P_R != P_R)
	{
	  printf("P_R is NAN\n");
	  return 1;
	}

      if(P_L != P_L)
	{
          printf("P_L is NAN\n");
          return 1;
        }

      if(u_L != u_L)
	{
          printf("U_L is NAN\n");
          return 1;
        }

      if(u_R != u_R)
        {
          printf("U_R is NAN\n");
          return 1;
        }

      if(rho_L != rho_L)
        {
          printf("rho_L is NAN\n");
          return 1;
        }
      
      if(rho_R != rho_R)
        {
          printf("rho_R is NAN\n");
          return 1;
        }

      if(S_L != S_L)
        {
          printf("S_L is NAN\n");
          return 1;
        }

      if(S_R != S_R)
        {
          printf("S_R is NAN\n");
          return 1;
        }


      double numerator = P_R - P_L + (rho_L * u_L * (S_L - u_L)) - (rho_R * u_R * (S_R - u_R));
      double denominator = (rho_L * (S_L - u_L)) - (rho_R * (S_R - u_R));

      if (numerator != numerator)
	{
	  if (denominator != denominator)
	    {
	      printf("NUMERATOR AND DENOMINATOR ARE NAN\n");
	    }
	  else
	    {
	      printf("JUST NUMERATOR IS NAN");
	    }
	  return 1;
	}

      if (denominator != denominator)
        {
          if (numerator != numerator)
            {
	      printf("NUMERATOR AND DENOMINATOR ARE NAN\n");
            }
          else
            {
	      printf("JUST DENOMINATOR IS NAN");
            }
	  return 1;
        }


      if (denominator == 0)
	{
	  printf("denominator is 0.\n");
	  return 1;
	}
		
      else
	{
	  S_star = numerator / denominator;
	}

      if (S_star != S_star)
	{
	  printf("S_STAR ERROR OTHER THAN DIVISION BY 0.\n");
	  return 1;
	}
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

      if (S_L - S_star == 0 || S_R - S_star == 0)
	{
	  printf("SOUND SPEED DIFFERENCE IS 0;\n");
	  return 1;
	}


      prefactor_L = rho_L * (S_L - u_L) / (S_L - S_star);
      prefactor_R = rho_R * (S_R - u_R) / (S_R - S_star);

      E_L = rho_L *(P_L / (rho_L*(gamma - 1)) + 0.5 * pow(u_L,2.0));
      E_R = rho_R *(P_R / (rho_R*(gamma - 1)) + 0.5 * pow(u_R,2.0));

      
      //      printf("RHO_L = %10.10f, RHO_R = %10.10f, U_L = %10.10f, U_R = %10.10f, E_L = %10.10f,E_R = %10.10f\n", rho_L, rho_R,u_L, u_R, E_L, E_R);
      

      Conserved_L_Star[0] = prefactor_L;
      Conserved_L_Star[1] = prefactor_L * S_star;
      Conserved_L_Star[2] = prefactor_L * ((E_L / rho_L) + (S_star - u_L) * (S_star + (P_L)/(rho_L*(S_L - u_L))));


      if (prefactor_L != prefactor_L)
	{
	  printf("PREFACTOR_L IS NAN\n");
	  return 1;
	}

      else
	{
	  //printf("Prefactor_L = %10.10f\n", prefactor_L);
	}
      

      Conserved_R_Star[0] = prefactor_R;
      Conserved_R_Star[1] = prefactor_R * S_star;
      Conserved_R_Star[2] = prefactor_R * ((E_R / rho_R) + (S_star - u_R) * (S_star + (P_R)/(rho_R*(S_R - u_R))));


      if (prefactor_R != prefactor_R)
	{
          //printf("PREFACTOR_R IS NAN\n");
	}
      else
	{
          //printf("Prefactor_R = %10.10f\n", prefactor_R);
        }


      Flux_R[0] = rho_R * u_R;
      Flux_R[1] = rho_R * u_R * u_R + P_R;
      Flux_R[2] = u_R * (E_R + P_R);
      


      FRs[i-2] = (Flux_R);
      

      Flux_L[0] = rho_L * u_L;
      Flux_L[1] = rho_L * u_L * u_L + P_L;
      Flux_L[2] = u_L * (E_L + P_L);



      FLs[i-2] = (Flux_L);

      if(Conserved_R_Star[0] != Conserved_R_Star[0] || Conserved_R_Star[1] != Conserved_R_Star[1] || Conserved_R_Star[2] != Conserved_R_Star[2])
	{
	  printf("UR STAR IS NAN\n");
	}
      
      Flux_R_Star[0] = Flux_R[0] + S_R * (Conserved_R_Star[0] - Conserved_R[0]);
      Flux_R_Star[1] = Flux_R[1] + S_R * (Conserved_R_Star[1] - Conserved_R[1]);
      Flux_R_Star[2] = Flux_R[2] + S_R * (Conserved_R_Star[2] - Conserved_R[2]);


      FRStars[i-2] = (Flux_R_Star);

      Flux_L_Star[0] = Flux_L[0] + S_L * (Conserved_L_Star[0] - Conserved_L[0]);
      Flux_L_Star[1] = Flux_L[1] + S_L * (Conserved_L_Star[1] - Conserved_L[1]);
      Flux_L_Star[2] = Flux_L[2] + S_L * (Conserved_L_Star[2] - Conserved_L[2]);


      FLStars[i-2] =(Flux_L_Star);
      
      int debugCount = 0;

      
    }
  printf("FUNDAMENTALS COMPUTED.\n");

  smax = tempSmax;
  //SignalSpeeds[0] = SignalSpeeds[1];
  //SignalSpeeds[mesh.NCells - 2] = SignalSpeeds[mesh.NCells - 3];
  ck[0]= -1.;
  ck[mesh.NCells-2] = 1.;
  if ((mesh.FirstDensityElement[3] - mesh.FirstDensityElement[2]) != 0)
    {
      rk[0] = (mesh.FirstDensityElement[2] - mesh.FirstDensityElement[1]) /(mesh.FirstDensityElement[3] - mesh.FirstDensityElement[2]);
    }
  else
    {
      rk[0] = 0;
    }

  if ((mesh.FirstDensityElement[mesh.NCells - 5] - mesh.FirstDensityElement[mesh.NCells-6]) != 0)
    {
      rk[mesh.NCells-2] = (mesh.FirstDensityElement[mesh.NCells - 5] - mesh.FirstDensityElement[mesh.NCells - 4]) /(mesh.FirstDensityElement[mesh.NCells - 5] - mesh.FirstDensityElement[mesh.NCells-6]);
    }

  else
    {
      rk[mesh.NCells-2]=0;
    }
  double dx = fabs(RiemannContext.X1 - RiemannContext.X0)/RiemannContext.MESH_RESOLUTION;
  double dt = RiemannContext.CFL_NUMBER * dx/smax;

  printf("EDGE VALUES SET.\n");
  for (int k = 1; k < mesh.NCells - 2; k++)
    {

      ck[k] = SignalSpeeds[k-1] * dt/dx;
      //printf("THE SIZE OF C IS %d AND THE SIZE OF SPEEDS IS %d\n ",(int)ck.size(),(int)SignalSpeeds.size());

      /*
	For some reason, reversing the sign order fixes the behaviour...
	2/19 JSB
      */
      if (ck[k] > 0 && (mesh.FirstDensityElement[k+1] - mesh.FirstDensityElement[k]) != 0)
	{
	  rk[k] = (mesh.FirstDensityElement[k] - mesh.FirstDensityElement[k-1])/(mesh.FirstDensityElement[k+1] - mesh.FirstDensityElement[k]);
	  printf("RK CONDITION 1 yields %10.10f\n",rk[k]);
	}

      else if (ck[k] <= 0 && (mesh.FirstDensityElement[k+1] - mesh.FirstDensityElement[k]) != 0.) 
	{
	  rk[k] = (mesh.FirstDensityElement[k+2] - mesh.FirstDensityElement[k+1])/(mesh.FirstDensityElement[k+1] - mesh.FirstDensityElement[k]);
          printf("RK CONDITION 2 yields %10.10f\n",rk[k]);

	}
      /*
      else
	{
	  rk[k] = 0.;
          printf("RK CONDITION 3 yields 0.\n");

	  }*/
      //ck[mesh.NCells-5] = 1.;                                                                                                               
      //rk[0] = (mesh.FirstDensityElement[2] - mesh.FirstDensityElement[3]) / (mesh.FirstDensityElement[3] - mesh.FirstDensityElement[4]);    
    }
  printf("CK AND RK ASSIGNED.\n");
  for (int i = 0; i < mesh.NCells - 3; i++)
    {
      double sgnc1, sgnc2, sgnc3, phi1,phi2,phi3;
      //printf("ASSIGNING FLUXES, ITERATION %d / %d\n",i,(int)fluxes.size());
      if (ck[i] != 0  && rk[i] != 0)
	{
	  sgnc1 = ck[ i] / fabs(ck[i ]);
	}
      else
	{
	  sgnc1 = 0.;
	}
      if (ck[i+1] != 0 && rk[i+1] != 0)
	{
	  sgnc2 = ck[i + 1]/fabs(ck[i + 1]);
	}
      else
	{
	  sgnc2 = 0.;
	}
      if(ck[i+1] != 0 && rk[i+2] != 0)
	{
	  sgnc3 = ck[i + 2] / fabs(ck[i + 2]);
	}
      else
	{
	  sgnc3 = 0.;
	}
      if (rk[i] == rk[i])
	{
	  phi1 = SUPERBEE(rk[i],fabs(ck[i]));
	}

      else
	{
	  phi1 = 0.0;
	}
      
      if (rk[i +1] == rk[i + 1])
	{
	  phi2 = SUPERBEE(rk[i + 1],fabs(ck[i + 1]));
	}
      else
	{
	  phi2 = 0.;
	}

      if(rk[i+2] == rk[i + 2])
	{
	  phi3 = SUPERBEE(rk[i + 2],fabs(ck[i + 2] ));
	}
      else
	{
	  phi3 = 0.;
	}
      printf("The limiters are [%10.10f, %10.10f, %10.10f]\n", phi1,phi2,phi3);
      printf("The signal speeds are [%10.10f, %10.10f, %10.10f]\n", ck[i],ck[i+1],ck[i+2]);
      //printf("The signs are [%10.10f, %10.10f, %10.10f]\n", sgnc1, sgnc2, sgnc3);
      if (ck[i] != ck[i] || ck[i+1]!=ck[i+1] || ck[i+2] != ck[i+2])
	{
	  printf("ck at iteration %d is NaN\n",i);
	  return 1;
	}

      if (rk[i] != rk[i] || rk[i+1]!=rk[i+1] || rk[i+2] != rk[i+2])
	{
	  printf("rk at iteration %d is NaN\n",i);
	  return 1;
        }


      if (FLs[i][0] != FLs[i][0] || FLs[i][1] !=FLs[i][1] || FLs[i][2] != FLs[i][2])
	{
	  printf("FL is NaN.\n");
	}

      if (FLStars[i][0] != FLStars[i][0] || FLStars[i][1] !=FLStars[i][1] || FLStars[i][2] != FLStars[i][2])
        {
          printf("FR is NaN.\n");
        }


      if (FRs[i][0] != FRs[i][0] || FRs[i][1] !=FRs[i][1] || FRs[i][2] != FRs[i][2])
	{
          printf("FR is NaN.\n");
	}

      if (FRStars[i][0] != FRStars[i][0] || FRStars[i][1] !=FRStars[i][1] || FRStars[i][2] != FRStars[i][2])
        {
          printf("FR is NaN.\n");
        }

      //printf("r = %10.10f at iteration %d\n", rk[i],i);
      //printf("c = %10.10f at iteration %d\n", ck[i],i);

      fluxes[i][0] = 0.5 * (FLs[i][0] + FRs[i][0]);
      fluxes[i][0] -= 0.5 * sgnc1 * phi1 * (FLStars[i][0] - FLs[i][0]);
      fluxes[i][0] -= 0.5 * sgnc2 * phi2 * (FRStars[i][0] - FLStars[i][0]);
      fluxes[i][0] -= 0.5 * sgnc3 * phi3 * (FRs[i][0] - FRStars[i][0]);


      fluxes[i][1] = 0.5 * (FLs[i][1] + FRs[i][1]);
      fluxes[i][1] -= 0.5 * sgnc1 * phi1 * (FLStars[i][1] - FLs[i][1]);
      fluxes[i][1] -= 0.5 * sgnc2 * phi2 * (FRStars[i][1] - FLStars[i][1]);
      fluxes[i][1] -= 0.5 * sgnc3 * phi3 * (FRs[i][1] - FRStars[i][1]);

      fluxes[i][2] = 0.5 * (FLs[i][2] + FRs[i][2]);
      fluxes[i][2] -= 0.5 * sgnc1 * phi1 * (FLStars[i][2] - FLs[i][2]);
      fluxes[i][2] -= 0.5 * sgnc2 * phi2 * (FRStars[i][2] - FLStars[i][2]);
      fluxes[i][2] -= 0.5 * sgnc3 * phi3 * (FRs[i][2] - FRStars[i][2]);

      printf("flux = [%10.10f, %10.10f, %10.10f]\n", fluxes[i][0], fluxes[i][1], fluxes[i][2]);
    }
    printf("Fluxes assigned.\n");
    return 0;
}


