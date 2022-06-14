#ifndef _MESHESOPERATIONSHEADER_
#include "MeshesOperations.h"
#endif

// [[Rcpp::export]]
Rcpp::NumericMatrix sampleMeshK(const unsigned nsims,
                                const Rcpp::List rmesh,
                                const bool triangulate) {
  Message("\u2014 Processing mesh...");
  Mesh3 mesh = makeSurfMesh<Mesh3, Point3>(rmesh, true, triangulate);
  Message("... done.\n");
  std::vector<Point3> sims;
  PMP::sample_triangle_mesh(
    mesh, std::back_inserter(sims),
    PMP::parameters::number_of_points_on_faces(nsims)
                    .do_sample_edges(false)
                    .do_sample_vertices(false)
  );
  Rcpp::NumericMatrix rsims = points3_to_matrix(sims);
  return Rcpp::transpose(rsims);
}
