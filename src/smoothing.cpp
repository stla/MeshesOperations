#ifndef _MESHESOPERATIONSHEADER_
#include "MeshesOperations.h"
#endif

// [[Rcpp::export]]
Rcpp::List smoothMeshK(const Rcpp::List rmesh,
                       const double angle,
                       const unsigned niters,
                       const bool triangulate,
                       const bool normals) {
  Mesh3 mesh0 = makeSurfMesh<Mesh3, Point3>(rmesh, true);
  if(triangulate) {
    const bool success = PMP::triangulate_faces(mesh0);
    if(!success) {
      const std::string msg = "Triangulation has failed.";
      Rcpp::stop(msg);
    }
  }
  // remove degenerate faces
  Mesh3 mesh = removeDegenerateFaces<Mesh3>(mesh0);
  // Constrain edges with a dihedral angle over the given angle
  typedef boost::property_map<Mesh3, CGAL::edge_is_feature_t>::type EIFMap;
  EIFMap eif = get(CGAL::edge_is_feature, mesh);
  PMP::detect_sharp_edges(mesh, angle, eif);
  int sharp_counter = 0;
  for(boost::graph_traits<Mesh3>::edge_descriptor e : mesh.edges()) {
    if(get(eif, e)) {
      ++sharp_counter;
    }
  }
  Rcpp::Rcout << sharp_counter << " sharp edges\n";
  Rcpp::Rcout << "Smoothing mesh... (" << niters << " iterations)\n";
  // Smooth with both angle and area criteria + Delaunay flips
  PMP::smooth_mesh(mesh,
                   PMP::parameters::number_of_iterations(niters)
                       .use_area_smoothing(false)
                       .use_Delaunay_flips(false)
                       .use_safety_constraints(false)  // authorize all moves
                       .edge_is_constrained_map(eif));
  Rcpp::List routmesh = RSurfTKMesh(mesh, normals, 0);
  return routmesh;
}

// [[Rcpp::export]]
Rcpp::List smoothShapeK(const Rcpp::List rmesh,
                        const double time,
                        const unsigned niters,
                        const bool triangulate,
                        const bool normals) {
  Message("\u2014 Processing mesh...");
  Mesh3 mesh = makeSurfMesh<Mesh3, Point3>(rmesh, true);
  if(triangulate) {
    Message("Triangulation.");
    const bool success = PMP::triangulate_faces(mesh);
    if(!success) {
      const std::string msg = "Triangulation has failed.";
      Rcpp::stop(msg);
    }
  }
  Message("... done.\n");
  std::set<Mesh3::Vertex_index> constrained_vertices;
  for(Mesh3::Vertex_index v : mesh.vertices()) {
    if(mesh.is_border(v)) {
      constrained_vertices.insert(v);
    }
  }
  const size_t nbv = constrained_vertices.size();
  std::string word;
  if(nbv > 1) {
    word = " border vertices.\n";
  } else {
    word = " border vertex.\n";
  }
  Rcpp::Rcout << "Constraining: " << nbv << word;
  CGAL::Boolean_property_map<std::set<Mesh3::Vertex_index>> vcmap(
      constrained_vertices);
  std::string tail;
  if(niters == 1) {
    tail = " iteration).\n";
  } else {
    tail = " iterations).\n";
  }
  Rcpp::Rcout << "Smoothing shape (" << niters << tail;
  PMP::smooth_shape(
      mesh, time,
      PMP::parameters::number_of_iterations(niters).vertex_is_constrained_map(
          vcmap));
  Rcpp::List routmesh = RSurfTKMesh(mesh, normals, 0);
  return routmesh;
}

//{}
