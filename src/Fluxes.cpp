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
  else if (c!=c)
    {
      return 0;
    }
  */
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
      if (/*fabs(u_L) + a_L*/ fabs(S_L) > tempSmax)
	{
	  //tempSmax = fabs(u_L) + a_L;
	  tempSmax = fabs(S_L);
	}
      if (/*fabs(u_R) + a_R*/fabs(S_R)  > tempSmax)
	{
          //tempSmax = fabs(u_R) + a_R;
	  tempSmax = fabs(S_R);
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
  vector<double> SignalSpeeds(mesh.NCells - 2, 0.0);
  vector<double> ck(mesh.NCells - 2, 0.0);
  vector<double> rk(mesh.NCells - 2, 0.0);
  vector<double> dq(mesh.NCells - 2,0.0);
  vector<vector<double> > FRs(mesh.NCells - 3, vector<double>(3,0.0));
  vector<vector<double> > FLs(mesh.NCells - 3,vector<double>(3,0.0));
  vector<vector<double> > FRStars(mesh.NCells - 3,vector<double>(3,0.0));
  vector<vector<double> > FLStars(mesh.NCells - 3,vector<double>(3,0.0));



  tempSmax = 0.;
  printf("COMPUTING FUNDAMENTALS...\n");

  //SignalSpeeds[mesh.NCells - 2] = 1.;
  //SignalSpeeds[0] = -1.;
  //SignalSpeeds[1] = -1.;
  //dq[0] = 0;
  //dq[mesh.NCells - 1] = 0;

  for (int i = 2; i < (mesh.NCells - 1); i++)
    {
      //printf("ITERATION %d\n",i);
      P_L = mesh.FirstPressureElement[i-1];
      P_R = mesh.FirstPressureElement[i];


      P_bar = 0.5 * (P_L + P_R);

      rho_L = mesh.FirstDensityElement[i-1];
      rho_R = mesh.FirstDensityElement[i];
      rho_bar = 0.5 * (rho_L + rho_R);
      u_L = mesh.FirstVelocityElement[i-1];
      u_R = mesh.FirstVelocityElement[i];

      a_L = pow((gamma * P_L / rho_L),0.5);
      a_R = pow((gamma * P_R / rho_R),0.5);
      a_bar = 0.5 * (a_L + a_R);


      p_star = P_bar - (0.5 * (u_R - u_L) * a_bar * rho_bar);


      if (p_star < 0)
        {  
          p_star = 0.;
        }


      dq[i -1] = rho_R - rho_L;
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


      S_L = u_L - a_L * q_L;
      S_R = u_R + a_R * q_R;
      SignalSpeeds[i-2] = S_L;

      // Compute max sound speed

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
	
      if (S_star != S_star)
	{
	  printf("S_STAR ERROR OTHER THAN DIVISION BY 0.\n");
	}

      // Construct conserved L/R and L/R star vectors from previously defined functions

      vector<double> primative_L(3);
      double primative_L_arr[3] = {rho_L, u_L, P_L};
      primative_L.assign(&primative_L_arr[0], &primative_L_arr[0] + 3);
      vector<double> primative_R(3);
      double primative_R_arr[3] = {rho_R, u_R, P_R};
      primative_R.assign(&primative_R_arr[0], &primative_R_arr[0] + 3);


      Conserved_R = primativeTo1DConservative(primative_R);

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
      


      FRs[i-2] = (Flux_R);
      

      Flux_L[0] = rho_L * u_L;
      Flux_L[1] = rho_L * u_L * u_L + P_L;
      Flux_L[2] = u_L * (E_L + P_L);



      FLs[i-2] = (Flux_L);

      
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
  ck[0]= -1.;
  ck[mesh.NCells-2] = 1.;

  double dx = fabs(RiemannContext.X1 - RiemannContext.X0)/RiemannContext.MESH_RESOLUTION;
  double dt = RiemannContext.CFL_NUMBER * dx/smax;

  printf("EDGE VALUES SET.\n");

  for (int k = 0; k < mesh.NCells - 3; k++)
    {
      if (/*dq[k+1] != 0 &&*/ FRs[0][k] == FRs[0][k] &&  FRs[1][k] == FRs[1][k] && FRs[2][k] == FRs[2][k]
	  && FRStars[0][k] == FRStars[0][k] && FRStars[1][k] == FRStars[1][k] && FRStars[2][k] == FRStars[2][k]
	  && FLStars[0][k] == FLStars[0][k] && FLStars[1][k] == FLStars[1][k] && FLStars[2][k] == FLStars[2][k])
	{
	  printf("DQ != 0\n");
	  ck[k+1] = SignalSpeeds[k+1] * dt/dx;
	  //printf("THE SIZE OF C IS %d AND THE SIZE OF SPEEDS IS %d\n ",(int)ck.size(),(int)SignalSpeeds.size());
	  
	  /*
	    For some reason, reversing the sign order fixes the behaviour...
	  */
	  if(ck[k+1] <= 0 && dq[k + 1] != 0)
	    {
	      rk[k+1] = dq[k]/dq[k+1];
	      printf("RK CONDITION 1 yields %10.10f\n",rk[k + 1]);
	    }
	  
	  else if (ck[k+1] > 0 && dq[k+1] != 0.) 
	    {
	      rk[k+1] = dq[k+2]/dq[k+1];
	      printf("RK CONDITION 2 yields %10.10f\n",rk[k+1]);
	      
	    }
      
	  else
	    {
	      rk[k+1] = 0.;
	      printf("RK CONDITION 3 yields 0.\n");
	      
	    }
	
	  //ck[mesh.NCells-5] = 1.;                                                                                                               
	  //rk[0] = (mesh.FirstDensityElement[2] - mesh.FirstDensityElement[3]) / (mesh.FirstDensityElement[3] - mesh.FirstDensityElement[4]);    
	
	  printf("CK AND RK ASSIGNED.\n");
	  //int i = k + 2;
	      
	  double sgnc1, sgnc2, sgnc3, phi1,phi2,phi3;
	  //printf("ASSIGNING FLUXES, ITERATION %d / %d\n",i,(int)fluxes.size());
	  
	  
	  if (ck[k] != 0  /*&& rk[i] != 0*/)
	    {
	      sgnc1 = ck[k] / fabs(ck[k]);
	    }
	  else
	    {
	      sgnc1 = 1.;
	    }
	  if (ck[k+1] != 0 /*&& rk[i+1] != 0*/)
	    {
	      sgnc2 = ck[k+1]/fabs(ck[k+1]);
	    }
	  else
	    {
	      sgnc2 = 1.;
	    }
	  if(ck[k+2] != 0 /*&& rk[i+2] != 0*/)
	    {
	      sgnc3 = ck[k+2 ] / fabs(ck[k+2]);
	    }
	  else
	    {
	      sgnc3 = 1.;
	    }
	  if (rk[k] == rk[k])
	    {
	      phi1 = MINBEE(rk[k],fabs(ck[k]));
	    }
	  
	  else
	    {
	      phi1 = 1.0;
	    }
	  
	  if (rk[k+1] == rk[k+1])
	    {
	      phi2 = MINBEE(rk[k+1],fabs(ck[k+ 1]));
	    }
	  else
	    {
	      phi2 = 1.;
	    }
	  
	  if(rk[k+2] == rk[k+2])
	    {
	      phi3 = MINBEE(rk[k+2],fabs(ck[k+2 ] ));
	    }
	  else
	    {
	      phi3 = 1.;
	    }
	  //printf("The limiters are [%10.10f, %10.10f, %10.10f]\n", phi1,phi2,phi3);
	  //printf("The signal speeds are [%10.10f, %10.10f, %10.10f]\n", ck[i],ck[i+1],ck[i+2]);
	  //printf("The signs are [%10.10f, %10.10f, %10.10f]\n", sgnc1, sgnc2, sgnc3);
	  /*
	  if (ck[i-2] != ck[i-2] || ck[i-1]!=ck[i-1] || ck[i] != ck[i])
	    {
	      printf("ck at iteration %d is NaN\n",i-2);
	      printf("The signal speeds are [%10.10f, %10.10f, %10.10f]\n", ck[i],ck[i+1],ck[i+2]); 
	      return 1;
	    }
	  
	  if (rk[i-2] != rk[i-2] || rk[i-1]!=rk[i-1] || rk[i] != rk[i])
	    {
	      printf("rk at iteration %d is NaN\n",i-2);
	      return 1;
	    }
	  
	      
	  if (FLs[i-2][0] != FLs[i-2][0] || FLs[i-2][1] !=FLs[i-2][1] || FLs[i-2][2] != FLs[i-2][2])
	    {
	      printf("FL is NaN.\n");
	      
	    }
	  
	  if (FLStars[i-2][0] != FLStars[i-2][0] || FLStars[i-2][1] !=FLStars[i-2][1] || FLStars[i-2][2] != FLStars[i-2][2])
	    {
	      printf("FLStars is NaN.\n");
	      //FLStars[i-2][0] = 0;//FLs[i-2][0];
	      //FLStars[i-2][1] = 0;//FLs[i-2][1];
	      //FLStars[i-2][2] = 0;//FLs[i-2][2];
	      //sgnc1 = 1;
	      //sgnc2 = 1;
	      //sgnc3 = 1;
	      //phi1 = 1;
	      //phi2 = 1;
	      //phi3 = 1;
	      
	    }
	  
	  
	  if (FRs[i-2][0] != FRs[i-2][0] || FRs[i-2][1] !=FRs[i-2][1] || FRs[i-2][2] != FRs[i-2][2])
	    {
	      printf("FR is NaN.\n");
	    }
	  
	  if (FRStars[i-2][0] != FRStars[i-2][0] || FRStars[i-2][1] !=FRStars[i-2][1] || FRStars[i-2][2] != FRStars[i-2][2])
	    {
	      printf("FR Stars is NaN.\n");
	      //FRStars[i-2][0] = 0;//FRs[i-2][0];
	      //FRStars[i-2][1] = 0; //FRs[i-2][1];
	      //FRStars[i-2][2] = 0;//FRs[i-2][2];
	      //sgnc1= 1;
	      //sgnc2= 1;
	      //sgnc3= 1;
	      //phi1 = 1;
	      //phi2 = 1;
	      //phi3 = 1;
	      
	      
	    }
	  */
	  printf("r = %10.10f at iteration %d\n", rk[k],k);
	  printf("c = %10.10f at iteration %d\n", ck[k],k);
	  
	  fluxes[k][0] = 0.5 * (FLs[k][0] + FRs[k][0]);
	  fluxes[k][0] -= 0.5 * sgnc1 * phi1 * (FLStars[k][0] - FLs[k][0]);
	  fluxes[k][0] -= 0.5 * sgnc2 * phi2 * (FRStars[k][0] - FLStars[k][0]);
	  fluxes[k][0] -= 0.5 * sgnc3 * phi3 * (FRs[k][0] - FRStars[k][0]);
	  
	  
	  fluxes[k][1] = 0.5 * (FLs[k][1] + FRs[k][1]);
	  fluxes[k][1] -= 0.5 * sgnc1 * phi1 * (FLStars[k][1] - FLs[k][1]);
	  fluxes[k][1] -= 0.5 * sgnc2 * phi2 * (FRStars[k][1] - FLStars[k][1]);
	  fluxes[k][1] -= 0.5 * sgnc3 * phi3 * (FRs[k][1] - FRStars[k][1]);
	  
	  fluxes[k][2] = 0.5 * (FLs[k][2] + FRs[k][2]);
	  fluxes[k][2] -= 0.5 * sgnc1 * phi1 * (FLStars[k][2] - FLs[k][2]);
	  fluxes[k][2] -= 0.5 * sgnc2 * phi2 * (FRStars[k][2] - FLStars[k][2]);
	  fluxes[k][2] -= 0.5 * sgnc3 * phi3 * (FRs[k][2] - FRStars[k][2]);

          printf("FR = [%10.10f, %10.10f, %10.10f]\n", FRs[k][0], FRs[k][1], FRs[k][2]);                                                         
          printf("FRStar = [%10.10f, %10.10f, %10.10f]\n", FRStars[k][0], FRStars[k][1], FRStars[k][2]);                                         
          printf("FL = [%10.10f, %10.10f, %10.10f]\n", FLs[k][0], FLs[k][1], FLs[k][2]);                                                         
          printf("FL = [%10.10f, %10.10f, %10.10f]\n", FLs[k][0], FLs[k][1], FLs[k][2]);                                                         
          printf("LIMITER = [%10.10f. %10.10f, %10.10f]\n", phi1, phi2, phi3);                                                                   
          printf("SGN = [%10.10f. %10.10f, %10.10f]\n", sgnc1, sgnc2, sgnc3);     
	  
	  /*
	    if (fluxes[i][0] < 0 || fluxes[i-2][1] < 0 || fluxes[i][2] < 0)
	    {
	    printf("WARNING: FLUX 2ND ELEMENT LESS THAN ZERO, ELEMENTS ARE:\n");
	    printf("FR = [%10.10f, %10.10f, %10.10f]\n", FRs[i][0], FRs[i][1], FRs[i][2]);
	    printf("FRStar = [%10.10f, %10.10f, %10.10f]\n", FRStars[i][0], FRStars[i][1], FRStars[i][2]);
	    printf("FL = [%10.10f, %10.10f, %10.10f]\n", FLs[i][0], FLs[i][1], FLs[i][2]);
	    printf("FL = [%10.10f, %10.10f, %10.10f]\n", FLs[i][0], FLs[i][1], FLs[i][2]);
	    printf("LIMITER = [%10.10f. %10.10f, %10.10f]\n", phi1, phi2, phi3);
	    printf("SGN = [%10.10f. %10.10f, %10.10f]\n", sgnc1, sgnc2, sgnc3);
	    
	    //fluxes[i][1] = 0;
	    }
	  */
	  //printf("flux = [%10.10f, %10.10f, %10.10f]\n", fluxes[i-2][0], fluxes[i-2][1], fluxes[i-2][2]);
	}
      
      else
	{
	  if (S_L >= 0)
	    {
	      fluxes[k] = Flux_L;
	      printf("Condition 1 tripped.\n");
	    }

	  else if (S_star >= 0 && S_L <= 0)
	    {
	      fluxes[k] = Flux_L_Star;
	      printf("Condition 2 tripped.\n");
	    }

	  else if (S_star <= 0 && S_R >= 0)
	    {
	      fluxes[k] = Flux_R_Star;
	      printf("Condition 3 tripped.\n");
	    }

	  else if (S_R <= 0)
	    {
	      fluxes[k] = Flux_R;
	      printf("Condition 4 tripped.\n");
	    }
	  //printf("flux = [%10.10f, %10.10f, %10.10f]\n", fluxes[k][0], fluxes[k][1], fluxes[k][2]);

	}
    
      //printf("Fluxes assigned.\n");
      //return 0;
    }
  for (int i = 0; i < (int)fluxes.size(); i++)
    {
     printf("flux = [%10.10f, %10.10f, %10.10f]\n", fluxes[i][0], fluxes[i][1], fluxes[i][2]);
   }
  printf("%d fluxes assigned.\n", (int)fluxes.size());                                                                                                              
  smax = tempSmax;
  return 0;
}


