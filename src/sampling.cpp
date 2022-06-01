#ifndef _MESHESOPERATIONSHEADER_
#include "MeshesOperations.h"
#endif

// [[Rcpp::export]]
Rcpp::NumericMatrix sampleMeshK(const unsigned nsims,
                                const Rcpp::List rmesh,
                                const bool triangulate) {
  Mesh3 mesh = makeSurfMesh<Mesh3, Point3>(rmesh, true);
  if(triangulate) {
    const bool success = PMP::triangulate_faces(mesh);
    if(!success) {
      const std::string msg = "Triangulation has failed.";
      Rcpp::stop(msg);
    }
  }
  std::vector<Point3> sims;
  PMP::sample_triangle_mesh(mesh, std::back_inserter(sims),
                            PMP::parameters::number_of_points_on_faces(nsims));
  Rcpp::NumericMatrix rsims = points3_to_matrix(sims);
  return rsims;
}
