#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <cassert>
#include <cstring>

#include "Omega_h.hpp"
#include "Omega_h_math.hpp"

Omega_h::Mesh read_vol_mesh(Omega_h::Library* lib, const char* vol_filename)
{
  double dummy;
  std::string line;
  std::fstream inmesh;
  int nv;

  inmesh.open(vol_filename, std::ios::in);
  while(true)
  {
    getline(inmesh,line);
    if (line == "surfaceelements") break;
  }
  int mtris;
  inmesh >> mtris;

  int tripoints[mtris][3];
  for (int i = 0; i < mtris; ++i)
  {
    inmesh >> dummy >> dummy >> dummy >> dummy >> dummy >> tripoints[i][0] >> tripoints[i][1] >> tripoints[i][2];
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
  for (int i = 0; i < medges; ++i)
  {
    inmesh >> surfid[i] >> dummy >> epoints[i][0] >> epoints[i][1] >> dummy >> dummy >> dummy >> dummy >> ednr1[i] >> dist1[i] >> ednr2[i] >> dist2[i];
  }

  while(true)
  {
    getline(inmesh,line);
    if (line == "points") break;
  }
  inmesh >> nv;
  double mpoints[nv][2];
  for (int i = 0; i < nv; ++i)
  {
    inmesh >> mpoints[i][0] >> mpoints[i][1] >> dummy;
  }
  inmesh.close();

  /* copy triangle connectivity */
  Omega_h::Write<Omega_h::LO> tv2v_w(mtris * 3);
  for (int i = 0; i < mtris; ++i)
    for (int j = 0; j < 3; ++j)
      tv2v_w[i * 3 + j] = tripoints[i][j] - 1;
  auto tv2v = Omega_h::LOs(tv2v_w);
  /* copy vertex coordinates */
  Omega_h::Write<Omega_h::Real> coords_w(nv * 2);
  for (int i = 0; i < nv; ++i)
    for (int j = 0; j < 2; ++j)
      coords_w[i * 2 + j] = mpoints[i][j];
  auto coords = Omega_h::Reals(coords_w);
  /* build the basic mesh from triangle connectivity and coordinates */
  Omega_h::Mesh mesh(lib);
  Omega_h::build_from_elems_and_coords(&mesh, 2, tv2v, coords);
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
  for (int i = 0; i < medges; ++i)
    for (int j = 0; j < 2; ++j)
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
  for (int i = 0; i < medges; ++i) {
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
  for (int i = 0; i < mesh.nverts(); ++i) {
    Omega_h::I8 vcd = 2; /* initialize dimension to 2D */
    Omega_h::I32 vci = -1;
    /* for each edge touch (ve) on this vertex */
    for (auto ve = v2ve[i]; ve < v2ve[i + 1]; ++ve) {
      auto e = ve2e[ve]; /* get the touching edge */
      auto ecd = edge_class_dim[e]; /* its classification */
      if (ecd == 2) continue; /* interior edge */
      auto eci = edge_class_id[e]; /* class_id of edge */
      if (vcd == 2) { /* touches first geometric edge */
        vcd = 1;
        vci = eci;
        continue;
      }
      /* touches second geometric edge, becomes geometric vertex */
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
  return mesh;
}

void read_and_attach_metric(Omega_h::Mesh* mesh, const char* metric_filename)
{
  /* reading in the anisotropy data*/
  std::fstream inmetric;
  inmetric.open(metric_filename, std::ios::in);
  if (!inmetric.is_open()) {
    std::cout << "could not open " << metric_filename << '\n';
    exit(-1);
  }
  int nv_metric, dim;
  inmetric >> nv_metric >> dim;
  /* Check to make sure the number of data poitns in the anisotropy informaiton is the same as the number of nodes*/
  if(nv_metric!=mesh->nverts() || dim!=3)
  {
    std::cout<<"Metric data does not correspond to the given mesh!"<<std::endl;
    exit(1);
  }
  /* The anisotropy variables are as follows:
  * aa corresponds to XX
  * bb corresponds to XY
  * cc corresponds to YY
  */
  double aa[nv_metric], bb[nv_metric], cc[nv_metric];
  for (int i = 0; i < nv_metric; ++i)
  {
    inmetric >> aa[i] >> bb[i] >> cc[i];
  }
  /* Anisotropy information read from file */

  /* attach the metric to the Omega_h mesh, as "target_metric" */
  Omega_h::Write<Omega_h::Real> metric_w(mesh->nverts() * 3);
  for (int i = 0; i < mesh->nverts(); ++i) {
    Omega_h::Matrix<2, 2> m;
    m[0][0] = aa[i];
    m[0][1] = m[1][0] = bb[i];
    m[1][1] = cc[i];
    Omega_h::set_symm(metric_w, i, m);
  }
  Omega_h::Reals metric(metric_w);
  mesh->add_tag(Omega_h::VERT, "target_metric", 3, OMEGA_H_METRIC,
      OMEGA_H_DO_OUTPUT, metric);
}

/* derive dist based on ednr !
 * This only works for the square model !
 */
double derive_dist(int ednr, Omega_h::Vector<2> x) {
  switch(ednr) {
    case 1: return x[0];
    case 2: return x[1];
    case 3: return 1.0 - x[0];
    case 4: return 1.0 - x[1];
  }
  Omega_h_fail("unknown ednr!");
}

void write_vol_mesh(Omega_h::Mesh* mesh, const char* vol_filename) {
  std::ofstream file(vol_filename);
  file << "mesh3d\ndimension\n2\ngeomtype\n0\n\n";
  file << "# surfnr    bcnr   domin  domout      np      p1      p2      p3\n";
  file << "surfaceelements\n" << mesh->nelems() << '\n';
  auto tv2v = mesh->ask_verts_of(Omega_h::TRI);
  auto w8 = std::setw(8);
  auto w12 = std::setw(12);
  auto w24 = std::setw(24);
  for (int i = 0; i < mesh->ntris(); ++i) {
    file << w8 << 2 << w8 << 1 << w8 << 0 << w8 << 0 << w8 << 3;
    for (int j = 0; j < 3; ++j)
      file << w8 << tv2v[i * 3 + j] + 1;
    file << '\n';
  }
  file << "\n\n#  matnr      np      p1      p2      p3      p4\n";
  file << "volumeelements\n0\n\n\n";
  auto ev2v = mesh->ask_verts_of(Omega_h::EDGE);
  auto e_class_dim = mesh->get_array<Omega_h::I8>(Omega_h::EDGE, "class_dim");
  auto e_class_id = mesh->get_array<Omega_h::I32>(Omega_h::EDGE, "class_id");
  auto coords = mesh->coords();
  int nsurfe = 0;
  for (int i = 0; i < mesh->nedges(); ++i)
    if (e_class_dim[i] == 1)
      nsurfe++;
  file << "# surfid       0      p1      p2";
  file << "    tri1    tri2 surfnr1 surfnr2";
  file << "   ednr1       dist1   ednr2       dist2\n";
  file << "edgesegmentsgi2\n";
  file << nsurfe << '\n';
  for (int i = 0; i < mesh->nedges(); ++i) {
    if (e_class_dim[i] != 1) continue;
    file << w8 << 1 << w8 << 0;
    for (int j = 0; j < 2; ++j)
      file << w8 << ev2v[i * 2 + j] + 1;
    file << w8 << -1 << w8 << -1 << w8 << 1 << w8 << 0;
    auto ednr = e_class_id[i];
    for (int j = 0; j < 2; ++j) {
      auto vert = ev2v[i * 2 + j];
      auto x = Omega_h::get_vector<2>(coords, vert);
      auto dist = derive_dist(ednr, x);
      file << w8 << ednr << w12 << dist;
    }
    file << '\n';
  }
  file << "\n\n#";
  for (int i = 0; i < 22; ++i) file << ' '; file << 'X';
  for (int i = 0; i < 23; ++i) file << ' '; file << 'Y';
  for (int i = 0; i < 23; ++i) file << ' '; file << 'Z';
  file << "\npoints\n" << mesh->nverts() << '\n';
  auto oldprecision = file.precision();
  file.precision(15);
  file << std::fixed;
  for (int i = 0; i < mesh->nverts(); ++i) {
    auto x = Omega_h::get_vector<2>(coords, i);
    double z = 0.0;
    file << w24 << x[0] << w24 << x[1] << w24 << z << '\n';
  }
  file.precision(oldprecision);
  file << std::defaultfloat;
  file << "materials\n1\n1 domain1\n\n\n";
  file << "# Surfnr";
  file << "          Red";
  file << "        Green";
  file << "         Blue";
  file << "\nface_colours\n1\n";
  file << "       1";
  file << "   0.00000000";
  file << "   1.00000000";
  file << "   0.00000000";
  file << "\n\n\nendmesh\n\n";
}

int main(int argc, char* argv[])
{
  auto lib = Omega_h::Library(&argc, &argv);

  const char* input_vol_file = nullptr;
  const char* metric_file = nullptr;
  const char* output_vol_file = nullptr;
  const char* adapt_log_dir = nullptr;
  const char* axes_file = nullptr;
  bool should_limit = false;
  double gradation_rate = 0.0;
  bool should_smooth = false;
  int nsmooth_iters = 0;

  for (int i = 1; i < argc; ++i) {
    if (!strcmp("--adapt_log", argv[i])) {
      if (i == argc - 1) {
        std::cout << "--adapt_log takes an argument\n";
        return -1;
      }
      adapt_log_dir = argv[++i];
    } else if (!strcmp("--axes", argv[i])) {
      if (i == argc - 1) {
        std::cout << "--axes takes an argument\n";
        return -1;
      }
      axes_file = argv[++i];
      std::cout << "axes file " << axes_file << '\n';
    } else if (!strcmp("--limit", argv[i])) {
      if (i == argc - 1) {
        std::cout << "--limit takes an argument\n";
        return -1;
      }
      should_limit = true;
      gradation_rate = atof(argv[++i]);
    } else if (!strcmp("--smooth", argv[i])) {
      if (i == argc - 1) {
        std::cout << "--smooth takes an argument\n";
        return -1;
      }
      should_smooth = true;
      nsmooth_iters = atoi(argv[++i]);
    } else if (!input_vol_file) {
      input_vol_file = argv[i];
      std::cout << "input vol file " << input_vol_file << '\n';
    } else if (!metric_file) {
      metric_file = argv[i];
      std::cout << "metric file " << metric_file << '\n';
    } else if (!output_vol_file) {
      output_vol_file = argv[i];
      std::cout << "output vol file " << output_vol_file << '\n';
    } else {
      std::cout << "unexpected argument " << argv[i] << '\n';
      return -1;
    }
  }
  if (!input_vol_file || !metric_file || !output_vol_file) {
    std::cout << "need to specify input .vol file, .mtr file, and output .vol file\n";
    return -1;
  }

  auto mesh = read_vol_mesh(&lib, input_vol_file);

  read_and_attach_metric(&mesh, metric_file);

  /* Find the "implied" metric: the one that keeps the mesh the same */
  auto implied_metric = find_implied_metric(&mesh);
  mesh.add_tag(Omega_h::VERT, "metric", 3, OMEGA_H_METRIC,
      OMEGA_H_DO_OUTPUT, implied_metric);
  /* Make sure initial qualities and lengths are attached */
  mesh.ask_qualities();
  mesh.ask_lengths();

  auto target = mesh.get_array<double>(Omega_h::VERT, "target_metric");

  if (should_smooth || should_limit) {
    auto elems_per_elem = Omega_h::expected_elems_per_elem_metric(&mesh, target);
    auto target_nelems = Omega_h::repro_sum(elems_per_elem);
    std::cout << "original metric would produce about " << target_nelems << " elements\n";
    mesh.add_tag(Omega_h::VERT, "original_metric", 3,
        OMEGA_H_DONT_TRANSFER, OMEGA_H_DO_OUTPUT, target);
  }

  if (should_smooth) {
    for (int i = 0; i < nsmooth_iters; ++i) {
      target = Omega_h::smooth_metric_once(&mesh, target);
    }
  }

  if (should_smooth && should_limit) {
    mesh.add_tag(Omega_h::VERT, "smoothed_metric", 3,
        OMEGA_H_DONT_TRANSFER, OMEGA_H_DO_OUTPUT, target);
  }

  if (should_limit) {
    target = Omega_h::limit_metric_gradation(&mesh, target, gradation_rate);
  }

  mesh.set_tag(Omega_h::VERT, "target_metric", target);

  auto elems_per_elem = Omega_h::expected_elems_per_elem_metric(&mesh, target);
  auto target_nelems = Omega_h::repro_sum(elems_per_elem);
  std::cerr << "target metric would produce about " << target_nelems << " elements\n";

  if (axes_file) {
    if (should_smooth || should_limit) {
      Omega_h::axes_from_metric_field(&mesh, "original_metric", "original_axis");
    }
    if (should_smooth && should_limit) {
      Omega_h::axes_from_metric_field(&mesh, "smoothed_metric", "smoothed_axis");
    }
    Omega_h::axes_from_metric_field(&mesh, "target_metric", "axis");
    Omega_h::vtk::write_vtu(axes_file, &mesh, 2);
  }

/* Adapt the mesh ! */
  Omega_h::vtk::FullWriter* writer = nullptr;
  if (adapt_log_dir) {
    writer = new Omega_h::vtk::FullWriter(&mesh, adapt_log_dir);
    writer->write();
  }
  Omega_h::AdaptOpts opts(&mesh);
  int n = 0;
  /* move metric closer to target until element quality too low */
  while (Omega_h::approach_size_field(&mesh, opts)) {
    Omega_h::adapt(&mesh, opts);
    if (adapt_log_dir) writer->write(); /* output VTK file */
    ++n;
    if (n == 100) {
      std::cerr << "aborting due to 100 approach iterations!\n";
      break;
    }
  }
  delete writer;

  write_vol_mesh(&mesh, output_vol_file);
}
