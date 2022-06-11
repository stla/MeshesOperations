#ifndef _MESHESOPERATIONSHEADER_
#include "MeshesOperations.h"
#endif

// [[Rcpp::export]]
Rcpp::NumericVector distanceK(const Rcpp::List rmesh,
                              const Rcpp::NumericMatrix points,
                              const bool triangulate) {
  Message("\u2014 Processing mesh...");
  Mesh3 mesh = makeSurfMesh<Mesh3, Point3>(rmesh, true);
  if(triangulate) {
    Message("Triangulation.");
    const bool success = PMP::triangulate_faces(mesh);
    if(!success) {
      Rcpp::stop("Triangulation has failed.");
    }
  }
  Message("... done.\n");
  std::vector<Point3> pts = matrix_to_points3<Point3>(points);
  const size_t npoints = points.ncol();
  Rcpp::NumericVector distances(npoints);
  for(size_t i = 0; i < npoints; i++){
    distances(i) = PMP::max_distance_to_triangle_mesh<CGAL::Sequential_tag>(
      pts, mesh
    );
  }
  return distances;
}
