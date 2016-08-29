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

  auto lib = Omega_h::Library(&argc, &argv);

  //opening netgen mesh file
  OMEGA_H_CHECK(argc == 4);
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

  /* copy triangle connectivity */
  Omega_h::Write<Omega_h::LO> tv2v_w(mtris * 3);
  for (i = 0; i < mtris; ++i)
    for (j = 0; j < 3; ++j)
      tv2v_w[i * 3 + j] = tripoints[i][j] - 1;
  auto tv2v = Omega_h::LOs(tv2v_w);
  /* copy vertex coordinates */
  Omega_h::Write<Omega_h::Real> coords_w(nv * 2);
  for (i = 0; i < nv; ++i)
    for (j = 0; j < 2; ++j)
      coords_w[i * 2 + j] = mpoints[i][j];
  auto coords = Omega_h::Reals(coords_w);
  /* build the basic mesh from triangle connectivity and coordinates */
  Omega_h::Mesh mesh;
  Omega_h::build_from_elems_and_coords(&mesh, lib, 2, tv2v, coords);
/* begin classification work */
  /* do a simple classification of triangles: all to interior #1 */
  mesh.add_tag(Omega_h::TRI, "class_dim", 1, OMEGA_H_INHERIT,
      OMEGA_H_DO_OUTPUT, Omega_h::Read<Omega_h::I8>(mesh.ntris(), 2));
  mesh.add_tag(Omega_h::TRI, "class_id", 1, OMEGA_H_INHERIT,
      OMEGA_H_DO_OUTPUT, Omega_h::Read<Omega_h::I32>(mesh.ntris(), 1));
  /* Omega_h constructs all edges internally. here we find matches
   * between those internally constructed edges and those specified
   * in "edgesegmentsgi2".
   * we assume ednr1=ednr2 and that this is the main indicator
   * of what boundary something is on.
   */
  /* copy "edgesegmentsgi2" connectivity */
  Omega_h::Write<Omega_h::LO> ev2v_w(medges * 2);
  for (i = 0; i < medges; ++i)
    for (j = 0; j < 2; ++j)
      ev2v_w[i * 2 + j] = epoints[i][j] - 1;
  auto ev2v = Omega_h::LOs(ev2v_w);
  /* use the Omega_h match finding function */
  auto v2e = mesh.ask_up(Omega_h::VERT, Omega_h::EDGE);
  Omega_h::LOs e2e;
  Omega_h::Read<Omega_h::I8> e2e_codes;
  Omega_h::find_matches(1, ev2v, mesh.ask_verts_of(Omega_h::EDGE), v2e,
      &e2e, &e2e_codes);
  /* initialize edge dimensions to 2D interior #1 */
  Omega_h::Write<Omega_h::I8> edge_class_dim_w(mesh.nedges(), 2);
  Omega_h::Write<Omega_h::I32> edge_class_id_w(mesh.nedges(), 1);
  /* edges that matched "edgesegmentsgi2" are set to 1D, i.e. boundary,
   * and their class_id is their ednr */
  for (i = 0; i < medges; ++i) {
    auto edge = e2e[i];
    edge_class_dim_w[edge] = 1;
    edge_class_id_w[edge] = ednr1[i];
  }
  Omega_h::Read<Omega_h::I8> edge_class_dim(edge_class_dim_w);
  Omega_h::Read<Omega_h::I32> edge_class_id(edge_class_id_w);
  /* set edge classification */
  mesh.add_tag(Omega_h::EDGE, "class_dim", 1, OMEGA_H_INHERIT,
      OMEGA_H_DO_OUTPUT, edge_class_dim);
  mesh.add_tag(Omega_h::EDGE, "class_id", 1, OMEGA_H_INHERIT,
      OMEGA_H_DO_OUTPUT, edge_class_id);
  Omega_h::Write<Omega_h::I8> vert_class_dim_w(mesh.nverts());
  Omega_h::Write<Omega_h::I32> vert_class_id_w(mesh.nverts());
  /* classify vertices: if it touches boundary edges with all the
   * same ednr, then it is 1D and has that ednr.
   * if it touches two ednrs, then it is a corner, class_dim is 0D,
   * and class_id will be the vertex id, because the initial mesh generator
   * creates corner points first.
   */
  auto v2ve = v2e.a2ab;
  auto ve2e = v2e.ab2b;
  /* for each vertex (i) */
  for (i = 0; i < mesh.nverts(); ++i) {
    Omega_h::I8 vcd = 2; /* initialize dimension to 2D */
    Omega_h::I32 vci = -1;
    /* for each edge touch (ve) */
    for (auto ve = v2ve[i]; ve < v2ve[i + 1]; ++ve) {
      auto e = ve2e[ve]; /* get the touching edge */
      auto ecd = edge_class_dim[e]; /* its classification */
      if (ecd == 2) continue; /* interior edge */
      auto eci = edge_class_id[e]; /* class_id of edge */
      if (vcd == 2) { /* though vertex was interior, but it touches boundary */
        vcd = 1;
        vci = eci;
        continue;
      }
      /* touches two geometric edges: geometric vertex */
      if (vcd == 1 && eci != vci) {
        vcd = 0;
        vci = i + 1;
      }
    }
    vert_class_dim_w[i] = vcd;
    vert_class_id_w[i] = vci;
  }
  Omega_h::Read<Omega_h::I8> vert_class_dim(vert_class_dim_w);
  Omega_h::Read<Omega_h::I32> vert_class_id(vert_class_id_w);
  /* set vertex classification */
  mesh.add_tag(Omega_h::VERT, "class_dim", 1, OMEGA_H_INHERIT,
      OMEGA_H_DO_OUTPUT, vert_class_dim);
  mesh.add_tag(Omega_h::VERT, "class_id", 1, OMEGA_H_INHERIT,
      OMEGA_H_DO_OUTPUT, vert_class_id);
/* end classification work */

  /* reading in the anisotropy data*/
  fstream inmetric;
  inmetric.open(argv[3], ios::in);
  int nv_metric, dim;
  inmetric >> nv_metric >> dim;
  /* Check to make sure the number of data poitns in the anisotropy informaiton is the same as the number of nodes*/
  if(nv_metric!=nv || dim!=3)
  {
    cout<<"Metric data does not correspond to the given mesh!"<<endl;
    exit(1);
  }
  i = 0;
  /* The anisotropy variables are as follows:
  * aa corresponds to XX
  * bb corresponds to XY
  * cc corresponds to YY
  */
  double aa[nv_metric], bb[nv_metric], cc[nv_metric];
  while(i<nv_metric)
  {
    inmetric >> aa[i] >> bb[i] >> cc[i];
    i++;
  }
  /* Anisotropy information read from file */

  /* attach the metric to the Omega_h mesh, as "target_metric" */
  Omega_h::Write<Omega_h::Real> metric_w(nv * 3);
  for (i = 0; i < nv; ++i) {
    Omega_h::Matrix<2, 2> m;
    m[0][0] = aa[i];
    m[0][1] = m[1][0] = bb[i];
    m[1][1] = cc[i];
    Omega_h::set_symm(metric_w, i, m);
  }
  Omega_h::Reals metric(metric_w);
  mesh.add_tag(Omega_h::VERT, "target_metric", 3, OMEGA_H_METRIC,
      OMEGA_H_DO_OUTPUT, metric);

  /* write the converted mesh to VTK file */
  Omega_h::vtk::write_vtu(argv[2], &mesh, 2);
}
