#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<math.h>
#include<assert.h>

#include "Omega_h.hpp"
#include "Omega_h_math.hpp"

using namespace std;
int main(int argc, char* argv[])
{
  int i,j,k, np;
  double dummy;
  string line;
  fstream inmesh, outmesh;
  int nv;

  //opening netgen mesh file
  OMEGA_H_CHECK(argc == 2);
  inmesh.open(argv[1], ios::in);
  while(true)
  {
    getline(inmesh,line);
    if (line == "surfaceelements") break;
  }
  int mtris, nvmax;
  inmesh >> mtris;

  int tripoints[mtris][3];
  nvmax=1;
  i = 0;
  while (i < mtris)
  {
    inmesh >> dummy >> dummy >> dummy >> dummy >> dummy >> tripoints[i][0] >> tripoints[i][1] >> tripoints[i][2];
    nvmax=max(max(tripoints[i][0], tripoints[i][1]), max(tripoints[i][2],nvmax));

    i++;
  }

  while(true)
  {
    getline(inmesh,line);
    if (line == "edgesegmentsgi2") break;
  }
  int medges;
  inmesh >> medges;
  int epoints[medges][2], ednr1[medges], ednr2[medges], surfid[medges];
  double dist1[medges], dist2[medges];
  int v_surfid[nvmax];
  for (i = 0; i < nvmax; ++i) v_surfid[i] = 0;
  {
    i=0;
  }
  while (i < medges)
  {
    inmesh >> surfid[i] >> dummy >> epoints[i][0] >> epoints[i][1] >> dummy >> dummy >> dummy >> dummy >> ednr1[i] >> dist1[i] >> ednr2[i] >> dist2[i];
    if (dist1[i] == 0 )
    {
       v_surfid[epoints[i][0] - 1] = epoints[i][0];
    }
    else
    {
       v_surfid[epoints[i][0] - 1] = surfid[i];
    }
    i++;
  }

  while(true)
  {
    getline(inmesh,line);
    if (line == "points") break;
  }
  inmesh >> nv;
  double mpoints[nv][2];
  i=0;
  while(i < nv)
  {
    inmesh >> mpoints[i][0] >> mpoints[i][1] >> dummy;
    i++;
  }
  inmesh.close();
}
