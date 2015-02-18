#include "../include/Context.h"
#include <vector>
#include <string>
#include <fstream>
#include <stdlib.h>

/*
  We should define methods for reading context from a file that should
  be given as an argument to the simulator.  Define similar methods as 
  the ones defined in Body.cpp, but we need not use a binary format since
  the number of entries SHOULD be small.
*/

using namespace std;

int delimitString(vector<string> &tokens, string line, string delimiter)
{
  int pos = 0;
  int newpos = 0;
  string token;
  while(1)
    {
      newpos = line.find(' ');
      token = line.substr(0,newpos);
      line.erase(0,newpos+1);

      tokens.push_back(token);
      if (newpos == string::npos)
        {
          break;
        }
    }
}



int readContext(char * PATH, context_t &context)
{
  ifstream file; 
  file.open(PATH,ifstream::in);
  vector<string> tokens;

  if (file.is_open())
    {
      int i=0;
      string line;
      while(getline(file,line))
        {
	  delimitString(tokens,line, " = ");
	  i++;
	}

      for (int j = 0; j < tokens.size(); j++)
	{
	  if(tokens[j] == "MESH_RESOLUTION")
            {
              context.MESH_RESOLUTION =  atoi(tokens[j + 2].c_str());
            }

	  if(tokens[j] == "RHO_L")
	    {
	      context.RHO_L = atof(tokens[j + 2].c_str());
	    }
         
	  if(tokens[j] == "RHO_R")
            {
              context.RHO_R = atof(tokens[j + 2].c_str());
            }

	  if (tokens[j] == "U_L")
	    {
	      context.U_L = atof(tokens[j + 2].c_str());
	    }

          if (tokens[j] == "U_R")
            {
              context.U_R = atof(tokens[j + 2].c_str());
            }


          if (tokens[j] == "P_L")
            {
              context.P_L = atof(tokens[j + 2].c_str());
            }

          if (tokens[j] == "P_R")
            {
              context.P_R = atof(tokens[j + 2].c_str());
            }


          if (tokens[j] == "CFL_NUMBER")
            {
              context.CFL_NUMBER = atof(tokens[j + 2].c_str());
            }

          if (tokens[j] == "EVOLVE_TIME")
            {
              context.EVOLVE_TIME = atof(tokens[j + 2].c_str());
            }

          if (tokens[j] == "X0")
            {
              context.X0 = atof(tokens[j + 2].c_str());
            }


          if (tokens[j] == "X1")
            {
              context.X1 = atof(tokens[j + 2].c_str());
            }

          if (tokens[j] == "INITIAL_DISCONTINUITY")
            {
              context.INITIAL_DISCONTINUITY = atof(tokens[j + 2].c_str());
            }



          if(tokens[j] == "SOLUTION_METHOD")
            {
              context.SOLUTION_METHOD = tokens[j + 2];
            }


          if(tokens[j] == "SNAPSHOT_PATH")
            {
              context.SNAPSHOT_PATH = tokens[j + 2];
            }

	  //cout << "Read token is " + tokens[j] << endl;
	}
    }
  file.close();

  return 0;
}


int printContext(context_t &context)
{
  printf("THE CONTEXT IS: \n");
  printf("MESH_RESOLUTION = %d\n",context.MESH_RESOLUTION);
  printf("X0 = %10.10f\n", context.X0);
  printf("X1 = %10.10f\n", context.X1);
  printf("INITIAL_DISCONTINUITY = %10.10f\n", context.INITIAL_DISCONTINUITY);
  printf("RHO_L = %10.10f\n", context.RHO_L);
  printf("P_L = %10.10f\n", context.P_L);
  printf("U_L = %10.10f\n", context.U_L);
  printf("RHO_R = %10.10f\n", context.RHO_R);
  printf("P_R = %10.10f\n", context.P_R);
  printf("U_R = %10.10f\n", context.U_R);

  printf("EVOLVE_TIME = %10.10f\n", context.EVOLVE_TIME);
  printf("CFL_NUMBER = %10.10f\n", context.CFL_NUMBER);

  cout << "SNAPSHOT_PATH = " << context.SNAPSHOT_PATH << '\n';
  cout << "SOLUTION_METHOD = " << context.SOLUTION_METHOD << '\n';

}
//delimitString(vector<string> &tokens, string line, string delimiter);
