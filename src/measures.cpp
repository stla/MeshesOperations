#ifndef _MESHESOPERATIONSHEADER_
#include "MeshesOperations.h"
#endif

// [[Rcpp::export]]
double meshVolumeK(const Rcpp::List rmesh, const bool triangulate) {
  Message("\u2014 Processing mesh...");
  Mesh3 mesh = makeSurfMesh<Mesh3, Point3>(rmesh, true, triangulate);
  Message("... done.\n");
  if(!CGAL::is_closed(mesh)) {
    Rcpp::stop("The mesh is not closed.");
  }
  const K::FT vol = PMP::volume(mesh);
  return CGAL::to_double<K::FT>(vol);
}

// [[Rcpp::export]]
double meshAreaK(const Rcpp::List rmesh, const bool triangulate) {
  Message("\u2014 Processing mesh...");
  Mesh3 mesh = makeSurfMesh<Mesh3, Point3>(rmesh, true, triangulate);
  Message("... done.\n");
  const K::FT ar = PMP::area(mesh);
  return CGAL::to_double<K::FT>(ar);
}