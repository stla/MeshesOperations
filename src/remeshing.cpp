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
  Message("\u2014 Processing mesh...");
  Mesh3 mesh = makeSurfMesh<Mesh3, Point3>(rmesh, true, triangulate);
  Message("... done.\n");
  PMP::isotropic_remeshing(
    mesh.faces(), targetEdgeLength, mesh,
    PMP::parameters::number_of_iterations(niters)
                    .number_of_relaxation_steps(nrelaxsteps)
  );
  mesh.collect_garbage();
  Rcpp::List routmesh = RSurfTKMesh(mesh, normals);
  return routmesh;
}

//{}
