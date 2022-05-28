#ifndef _MESHESOPERATIONSHEADER_
#include "MeshesOperations.h"
#endif

// [[Rcpp::export]]
Rcpp::List smoothMeshK(const Rcpp::List rmesh,
                       const double angle,
                       const unsigned niters,
                       const bool triangulate,
                       const bool normals) {
  const Mesh3 mesh = makeSurfMesh<Mesh3, Point3>(rmesh, true);
  Rcpp::IntegerMatrix Edges0;
  Rcpp::NumericMatrix Normals0;
  if(triangulate) {
    Edges0 = getEdges2<K, Mesh3, Point3>(mesh, 0);
    if(normals) {
      Normals0 = getKNormals(mesh);
    }
    const bool success = PMP::triangulate_faces(mesh);
    if(!success) {
      const std::string msg = "Triangulation has failed.";
      Rcpp::stop(msg);
    }
  }
  // remove degenerate faces
  mesh = removeDegenerateFaces<Mesh3>(mesh);
  // Constrain edges with a dihedral angle over the given angle
  typedef boost::property_map<Mesh3, CGAL::edge_is_feature_t>::type EIFMap;
  EIFMap eif = get(CGAL::edge_is_feature, mesh);
  PMP::detect_sharp_edges(mesh, angle, eif);
  int sharp_counter = 0;
  for(Mesh3::Edge_index e : mesh.edges()) {
    if(get(eif, e)) {
      ++sharp_counter;
    }
  }
  std::cout << sharp_counter << " sharp edges" << std::endl;
  std::cout << "Smoothing mesh... (" << niters << " iterations)" << std::endl;
  // Smooth with both angle and area criteria + Delaunay flips
  PMP::smooth_mesh(mesh,
                   PMP::parameters::number_of_iterations(niters)
                       .use_safety_constraints(false)  // authorize all moves
                       .edge_is_constrained_map(eif));
  Rcpp::List routmesh = RSurfKMesh(mesh, normals, 0);
  if(triangulate) {
    routmesh["edges0"] = Edges0;
    if(normals) {
      routmesh["normals0"] = Normals0;
    }
  }
  return routmesh;
}

//{}
