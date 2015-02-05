                                                                                                         
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

#pragma once
struct context_t
{
  // domain information
  int MESH_RESOLUTION;
  double X0;
  double X1;

  // initial conditions
  double RHO_L;
  double RHO_R;
  double U_L;
  double U_R;
  double P_L;
  double P_R;
  double INITIAL_DISCONTINUITY;

  // Simulation information

  double CFL_NUMBER;
  double EVOLVE_TIME;
  string SOLUTION_METHOD;
  string SNAPSHOT_PATH;
};

int readContext(char * PATH, context_t &context);
int printContext(context_t &context);
