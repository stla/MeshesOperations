#ifndef _MESHESOPERATIONSHEADER_
#include "MeshesOperations.h"
#endif

// [[Rcpp::export]]
Rcpp::List isotropicRemeshingK(const Rcpp::List rmesh,
                               const double targetEdgeLength,
                               const unsigned niters,
                               const unsigned nrelaxsteps,
                               const bool triangulate,
                               const bool normals) {
  Mesh3 mesh = makeSurfMesh<Mesh3, Point3>(rmesh, true);
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
  PMP::isotropic_remeshing(
      mesh.faces(), targetEdgeLength, mesh,
      PMP::parameters::number_of_iterations(niters).number_of_relaxation_steps(
          nrelaxsteps));
  mesh.collect_garbage();
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
