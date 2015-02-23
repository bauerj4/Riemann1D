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

double DUMMYBEE(double r, double c)
{
  return 1.0;
}

double VLBEE(double r, double c)
{
  if (r <= 0)
    {
      return 1.;
    }
  else 
    {
      return (1. - 2. * r *(1. - c * r)/(1. + r));
    }
}
double SUPERBEE(double r, double c)
{
  if (r <= 0 )
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
  /*
  else if (r!=r)
    {
      return 0;
      }*/
  
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
	  //P_L = 0;
	  printf("WARNING: P_L < 0\n");
	}

      if (P_R < 0)
	{
          //P_R =0;
          printf("WARNING: P_R = %10.10f  < 0\n", P_R);
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
	  //tempSmax = fabs(S_L);
	}
      if (fabs(u_R) + a_R  > tempSmax)
	{
          tempSmax = fabs(u_R) + a_R;
	  //tempSmax = fabs(S_R);
        }



      double numerator = P_R - P_L + (rho_L * u_L * (S_L - u_L)) - (rho_R * u_R * (S_R - u_R));
      double denominator = (rho_L * (S_L - u_L)) - (rho_R * (S_R - u_R));

      if (denominator == 0)
	{
	  printf("DENOMINATOR => DIVISION BY 0\n");
	  //return 1;
	}
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
      
      //int debugCount = 0;
      if (S_L >= 0)
	{
	  fluxes[i-1] = Flux_L;
	  printf("Condition 1 tripped.\n");
	  //debugCount++;
	}

      else if (S_star >= 0 && S_L <= 0)
	{
	  fluxes[i-1] = Flux_L_Star;
          printf("Condition 2 tripped.\n");
	  //debugCount++;
	}

      else if (S_star <= 0 && S_R >= 0)
	{
	  fluxes[i-1] = Flux_R_Star;
          printf("Condition 3 tripped.\n");
	  //debugCount++;
	}

      else if (S_R <= 0)
	{
	  fluxes[i-1] = Flux_R;
          printf("Condition 4 tripped.\n");
	  //debugCount++;
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

  /*
    ON A NEW THOUGHT, WHERE DQ IS ZERO WE SHOULD COMPUTE A FIRST ORDER FLUX.  WILL REDO TOMORROW.
  */

  /*
    REDO WITH WEIGHTED AVERAGE STATE
    Most of the numerical issues I'm seeing are because I am computing fluxes directly.
    It would be preferable to do the one flux calculation as Toro suggests. 
  */

  // From Toro 2009
  double gamma = 1.4; // Given.
  
  if (fluxes.size() != (mesh.NCells - 3)) // Should confirm that this is the correct number
    {
      printf("INVALID FLUX VECTOR SIZE!\n");
      return 1;
    }
  // Compute quantities in starred region
  
  double u_L, u_R, p_star, P_bar, P_L, P_R, rho_bar, rho_L, rho_R, tempSmax,
    a_bar, a_L, a_R, q_L, q_R, S_L, S_R, S_star, E_L, E_R, prefactor_L, prefactor_R,
    P_starR, P_starL, rho_starR, rho_starL, u_star;

  vector<double> Conserved_R(3,0);
  vector<double> Conserved_L(3,0);
  vector<double> Conserved_R_Star(3,0);
  vector<double> Conserved_L_Star(3,0);
  vector<double> Flux_R(3,0);
  vector<double> Flux_L(3,0);
  vector<double> Flux_R_Star(3,0);
  vector<double> Flux_L_Star(3,0);
  vector<double> SignalSpeeds(mesh.NCells - 2, 0.0);
  vector<double> ck(mesh.NCells-2, 0.0);
  vector<double> rk(mesh.NCells-3, 0.0);
  vector<double> dq(mesh.NCells-1,0.0);
  vector<double> PrimativeR(3,0.0);
  vector<double> PrimativeL(3,0.0);
  vector<double> PrimativeRStar(3,0.0);
  vector<double> PrimativeLStar(3,0.0);
  vector<vector<double> > PrimativeRs(mesh.NCells - 2, vector<double>(3,0.0));
  vector<vector<double> > PrimativeLs(mesh.NCells - 2,vector<double>(3,0.0));
  vector<vector<double> > PrimativeRStars(mesh.NCells - 2,vector<double>(3,0.0));
  vector<vector<double> > PrimativeLStars(mesh.NCells - 2,vector<double>(3,0.0));
  vector<vector<double> > PrimativeBar(mesh.NCells - 2, vector<double>(3,0.0));

  tempSmax = 0.;
  printf("COMPUTING FUNDAMENTALS...\n");

  //SignalSpeeds[mesh.NCells - 2] = 1.;
  //SignalSpeeds[0] = -1.;
  //SignalSpeeds[1] = -1.;
  //dq[0] = 0;
  //dq[mesh.NCells - 1] = 0;
  dq[0] = mesh.FirstDensityElement[1] - mesh.FirstDensityElement[0];
  dq[mesh.NCells - 2] = mesh.FirstDensityElement[mesh.NCells - 1] - mesh.FirstDensityElement[mesh.NCells - 2];

  for (int i = 1; i < (mesh.NCells - 1); i++)
    {
      //printf("ITERATION %d\n",i);
      P_L = mesh.FirstPressureElement[i-1];
      P_R = mesh.FirstPressureElement[i];


      P_bar = 0.5 * (P_L + P_R);

      rho_L = mesh.FirstDensityElement[i-1];
      rho_R = mesh.FirstDensityElement[i];
      rho_bar = 0.5 * (rho_L + rho_R);
      

      if(rho_L < 0 || rho_L != rho_L)
	{
	  printf("UNPHYSICAL RHO_L %10.10f.\n", rho_L);
	  //rho_L = fabs(rho_L)/2.;
	  return 1;
	  //rho_L = rho_bar;
	}

      if(rho_R < 0 || rho_R != rho_R)
	{
          printf("UNPHYSICAL RHO_R %10.10f.\n", rho_R);
	  //rho_R = fabs(rho_R)/2.;
	  return 1;
          //rho_R = rho_bar;
        }

      if(P_L < 0 || P_L != P_L)
        {
          printf("UNPHYSICAL P_L %10.10f.\n", P_L);
	  return 1;
          //P_L = fabs(P_L) / 2.;
	  //if (P_R > 0 && P_R < P_L)
	  //  {
	  //    P_L = P_R;
	  //  }
        }


      if(P_R < 0 || P_R != P_R)
        {
          printf("UNPHYSICAL P_R %10.10f.\n", P_R);
	  return 1;
          //P_R = fabs(P_R)/2.;

          //if (P_L > 0 && P_L < P_R)
          //  {
	  //    P_R = P_L;
	  // }

        }



      //rho_bar = 0.5 * (rho_L + rho_R);
      u_L = mesh.FirstVelocityElement[i-1];
      u_R = mesh.FirstVelocityElement[i];

      if (rho_L != 0 && P_L > 0)
	{
	  a_L = pow((gamma * P_L / rho_L),0.5);
	}
      else
	{

	  a_L = 0;
	}

      if (rho_R != 0 && P_R > 0)
	{
	  a_R = pow((gamma * P_R / rho_R),0.5);
	}
      else
	{
	  a_R = 0;
	}

      a_bar = 0.5 * (a_L + a_R);


      p_star = P_bar - (0.5 * (u_R - u_L) * a_bar * rho_bar);

      if (p_star < 0 || p_star != p_star)
        {  
	  printf("UNPHYSICAL P_STAR %10.10f\n",p_star);
          p_star = 0.;
        }


      dq[i] = rho_R - rho_L;
      printf("DQ %d = %10.10f\n",i-1,dq[i-1]);

      // The a's represent the sound speeds.
      // Compute as consequence of EoS

      if (p_star <= P_L || P_L == 0)
	{
	  q_L = 1.;
	}

      else if (P_L > p_star)
	{
	  double num = 1 + ((gamma + 1)/(2 * gamma)) * ((p_star/P_L) - 1);
	  q_L = pow(num,0.5);
	}

      if (p_star <= P_R || P_R == 0)
	{
          q_R =1.;
	}
      else if(P_R > p_star)
	{
          double num = 1. + ((gamma + 1.)/(2. * gamma)) * ((p_star/P_R) - 1.);
          q_R =pow(num,0.5);
	}

      // WAVE SPEED ESTIMATES
      // I have chosen to deviate from Toro because of numerical instabilities.
      // 
     
      //double R = pow(rho_R/rho_L, 0.5);
      // double u_tilde = (u_L +  R * u_R)/ (1 + R);

      S_L = u_L - a_L * q_L;
      S_R = u_R + a_R * q_R;
      P_starL = P_L + rho_L *(S_L - u_L) * (S_star - u_L);
      P_starR =P_R + rho_R *(S_L - u_R) * (S_star - u_R);


      printf("(S_L, S_R) = (%10.10f, %10.10f)\n", S_L, S_R );
      SignalSpeeds[i-1] = S_R;
      // if(S_R < S_L)
      //{
      //  SignalSpeeds[i-2] = S_R;
      //}

      // Compute max signal speed
      /*
      if (fabs(u_L) + a_L > tempSmax)
	{
	  tempSmax = fabs(u_L) + a_L;
	}
      if (fabs(u_R) + a_R > tempSmax)
	{
          tempSmax = fabs(u_R) + a_R;
        }
      */
      if (fabs(u_L) + a_L > tempSmax)
        {
          tempSmax = fabs(u_L) + a_L;
          //tempSmax = fabs(S_L);                                                                                                                 
        }
      if (fabs(u_R) + a_R  > tempSmax)
        {
          tempSmax = fabs(u_R) + a_R;
          //tempSmax = fabs(S_R);                                                                                                                 
        }



      double numerator = P_R - P_L + (rho_L * u_L * (S_L - u_L)) - (rho_R * u_R * (S_R - u_R));
      double denominator = (rho_L * (S_L - u_L)) - (rho_R * (S_R - u_R));
      if (denominator == 0)
	{
	  printf("DENOMINATOR = 0\n");
	  S_star = u_R;
	}
      else
	{
	  S_star = numerator / denominator;
	}
      //SignalSpeeds[i-2] = S_star;

      if (S_star != S_star)
	{
	  printf("S_STAR ERROR OTHER THAN DIVISION BY 0.\n");
	}

      u_star = S_star;

      rho_starR = (S_R - u_R)/(S_R - S_star) * rho_R; 
      rho_starL = (S_L - u_L)/(S_L - S_star) * rho_L;


      PrimativeR[0] = rho_R;
      PrimativeR[1] = u_R;
      PrimativeR[2] = P_R;

      PrimativeL[0] = rho_L;
      PrimativeL[1] = u_L;
      PrimativeL[2] = P_L;

      PrimativeRStar[0] = rho_starR;
      PrimativeRStar[1] = u_star;
      PrimativeRStar[2] = P_starR;

      PrimativeLStar[0] = rho_starL;
      PrimativeLStar[1] = u_star;
      PrimativeLStar[2] = P_starL;

      
      PrimativeRs[i-1] = PrimativeR;
      PrimativeLs[i-1] = PrimativeL;
      PrimativeRStars[i-1] = PrimativeRStar;
      PrimativeLStars[i-1] = PrimativeLStar;


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
      
    }
  printf("FUNDAMENTALS COMPUTED.\n");

  smax = tempSmax;
  ck[0]= -1.;
  ck[mesh.NCells-2] = 1.;

  double dx = fabs(RiemannContext.X1 - RiemannContext.X0)/RiemannContext.MESH_RESOLUTION;
  double dt = RiemannContext.CFL_NUMBER * dx/smax;

  printf("EDGE VALUES SET.\n");
  printf("ASSIGNING WEIGHTS...\n");

  for (int k =0; k < mesh.NCells - 2; k++)
    {
      ck[k+1] = SignalSpeeds[k] * dt/dx;
    }



  for (int k = 0; k < mesh.NCells - 3; k++)
    {
	
      printf("DQ != 0\n");
      /*
      if (k!=0)
	{
	  ck[k] = SignalSpeeds[k] * dt/dx;
	}
	
      ck[k+1] = SignalSpeeds[k+1] * dt/dx;
      if (k != mesh.NCells - 1)
	{
	  ck[k+2] = SignalSpeeds[k+2] * dt/dx;
	  }*/


      //ck[k+3] = SignalSpeeds[k+3] * dt/dx;
      
      //printf("THE SIZE OF C IS %d AND THE SIZE OF SPEEDS IS %d\n ",(int)ck.size(),(int)SignalSpeeds.size());
      
      /*
	For some reason, reversing the sign order fixes the behaviour...
      */
      if(ck[k+1] >= 0 && dq[k + 2] != 0)
	{
	  rk[k] = dq[k+1]/dq[k+2];
	  //printf("RK CONDITION 1 yields %10.10f\n",rk[k+1]);
	}
      
      else if (ck[k+1] < 0 && dq[k+2] != 0.) 
	{
	  rk[k] = dq[k+3]/dq[k+2];
	  // printf("RK CONDITION 2 yields %10.10f\n",rk[k+1]);
	  
	}
      
      else
	{
	  rk[k] = 0.;
	  printf("RK CONDITION 3 yields 0.\n");
	  
	}
    }

  // Check in morning.  Issue is most likely here

  for(int k = 0; k < mesh.NCells - 3; k++)
    {
      
      //ck[mesh.NCells-5] = 1.;                                                                                                               
      //rk[0] = (mesh.FirstDensityElement[2] - mesh.FirstDensityElement[3]) / (mesh.FirstDensityElement[3] - mesh.FirstDensityElement[4]);    
      
      printf("CK AND RK ASSIGNED.\n");
      //int i = k + 2;
      
      double sgnc1, sgnc2, sgnc3, phi1,phi2,phi3;
      //printf("ASSIGNING FLUXES, ITERATION %d / %d\n",i,(int)fluxes.size());
      
      
      if (ck[k+1] != 0  /*&& rk[i] != 0*/)
	{
	  sgnc1 = ck[k+1] / fabs(ck[k+1]);
	}
      else if (ck[k+1] == 0)
	{
	  sgnc1 = 1.;
	}
      else
	{
	  sgnc1 = 0;
	}
      if (ck[k+2] != 0 /*&& rk[i+1] != 0*/)
	{
	  sgnc2 = ck[k+2]/fabs(ck[k+2]);
	}
      else if (ck[k+1] == 0)
	{
	  sgnc2 = 1.;
	}
      else
	{
	  sgnc2 = 0;
	}
      if(ck[k+3] != 0 /*&& rk[i+2] != 0*/)
	{
	  sgnc3 = ck[k+3 ] / fabs(ck[k+3]);
	}
      else if(ck[k+2] == 0)
	{ 
	  sgnc3 = 1.;
	}
      else
	{
	  sgnc3 = 0;
	}
      if (rk[k+1] == rk[k+1])
	{
	  phi1 = SUPERBEE(rk[k+1],fabs(ck[k+1]));
	}
	  
      else
	{
	  phi1 = 1.0;
	}
      
      if (rk[k+1] == rk[k+1])
	{
	  phi2 = SUPERBEE(rk[k+1],fabs(ck[k+ 2]));
	}
      else
	{
	  phi2 = 1.;
	}
      
      if(rk[k+1] == rk[k+1])
	{
	  phi3 = SUPERBEE(rk[k+1],fabs(ck[k+3] ));
	}
      else
	{
	  phi3 = 1.;
	}
      /*
      double beta, sum;
      for (int j = 0; j < mesh.NCells-3; j++)
	{
	   beta = ck[j + 1] - ck[j];
	   sum += beta;
	}
      sum = sum/2.;
      if (sum != )
      */
      PrimativeBar[k][0] = 0.5 * (PrimativeRs[k+1][0] + PrimativeLs[k+1][0]);
      PrimativeBar[k][0] -= 0.5 * sgnc1 * phi1 * (PrimativeLStars[k+1][0] - PrimativeLs[k+1][0]);
      PrimativeBar[k][0] -= 0.5 * sgnc2 * phi2 * (PrimativeRStars[k+1][0] - PrimativeLStars[k+1][0]);
      PrimativeBar[k][0] -= 0.5 * sgnc3 * phi3 * (PrimativeRs[k+1][0] - PrimativeRStars[k+1][0]);

      PrimativeBar[k][1] = 0.5 * (PrimativeRs[k+1][1] + PrimativeLs[k+1][1]);
      PrimativeBar[k][1] -= 0.5 * sgnc1 * phi1 * (PrimativeLStars[k+1][1] - PrimativeLs[k+1][1]);
      PrimativeBar[k][1] -= 0.5 * sgnc2 * phi2 * (PrimativeRStars[k+1][1] - PrimativeLStars[k+1][1]);
      PrimativeBar[k][1] -= 0.5 * sgnc3 * phi3 * (PrimativeRs[k+1][1] - PrimativeRStars[k+1][1]);


      PrimativeBar[k][2] = 0.5 * (PrimativeRs[k+1][2] + PrimativeLs[k+1][2]);
      PrimativeBar[k][2] -= 0.5 * sgnc1 * phi1 * (PrimativeLStars[k+1][2] - PrimativeLs[k+1][2]);
      PrimativeBar[k][2] -= 0.5 * sgnc2 * phi2 * (PrimativeRStars[k+1][2] - PrimativeLStars[k+1][2]);
      PrimativeBar[k][2] -= 0.5 * sgnc3 * phi3 * (PrimativeRs[k+1][2] - PrimativeRStars[k+1][2]);


      double int_energy = PrimativeBar[k][2] / (PrimativeBar[k][0] * (gamma - 1.));

      if (int_energy!=int_energy)
	{
	  int_energy = 0;
	}

      double kinetic_energy = 0.5 * PrimativeBar[k][1] * PrimativeBar[k][1];
      double energy = PrimativeBar[k][0] * (int_energy + kinetic_energy);


      fluxes[k][0] = PrimativeBar[k][0] * PrimativeBar[k][1];
      fluxes[k][1] = PrimativeBar[k][0] * PrimativeBar[k][1] * PrimativeBar[k][1] + PrimativeBar[k][2];
      fluxes[k][2] = PrimativeBar[k][1] * (energy + PrimativeBar[k][2]);

      if (fluxes[k][1] < 0)
	{
	  printf("UNPHYSICAL FLUX.\n");
	  //fluxes[k][1] = 0;
	}

      if (fluxes[k][2] < 0)
	{ 
	  printf("UNPHYSICAL FLUX.\n");
	  //fluxes[k][2] = 0;
	}

      
      if (fluxes[k][1] < 0 )
	{
	  if (S_L >= 0)
	    {
	      fluxes[k] = Flux_L;
	      printf("Condition 1 tripped.\n");
	      //debugCount++;                                                                                                                          
	    }

	  else if (S_star >= 0 && S_L <= 0)
	    {
	      fluxes[k] = Flux_L_Star;
	      printf("Condition 2 tripped.\n");
	      //debugCount++;                                                                                                                          
	    }

	  else if (S_star <= 0 && S_R >= 0)
	    {
	      fluxes[k] = Flux_R_Star;
	      printf("Condition 3 tripped.\n");
	      //debugCount++;                                                                                                                          
	    }

	  else if (S_R <= 0)
	    {
	      fluxes[k] = Flux_R;
	      printf("Condition 4 tripped.\n");
	      //debugCount++;                                                                                                                          
	    }

	}
	 
    }
  for (int i = 0; i < (int)fluxes.size(); i++)
    {
      printf("flux = [%10.10f, %10.10f, %10.10f]\n", fluxes[i][0], fluxes[i][1], fluxes[i][2]);
      printf("cks = [%10.10f, %10.10f, %10.10f]\n", ck[i], ck[i+1], ck[i+2]);
      printf("rks = [%10.10f, %10.10f, %10.10f]\n", rk[i], rk[i], rk[i]);
    }
  printf("%d fluxes assigned.\n", (int)fluxes.size());                                                                                                              
  double beta, sum;
  sum = 0;
  beta = 0;

  for (int j = 0; j < mesh.NCells - 1; j++)
    {
      beta = ck[j + 1] - ck[j];
      printf("c_%d = %10.10f\n", j-2, ck[j]);
      sum += beta;
    }
  //sum /=2.0; //sum;// / 2.;

  if (sum != 1)
    {
      printf("IMPROPER WEIGHTING %10.10f.\n", sum);
      //return 1;
    }

  smax = tempSmax;
  return 0;
  
}


