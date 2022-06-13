#ifndef _MESHESOPERATIONSHEADER_
#include "MeshesOperations.h"
#endif

// [[Rcpp::export]]
double meshVolumeK(const Rcpp::List rmesh, const bool triangulate) {
  Mesh3 mesh = makeSurfMesh<Mesh3, Point3>(rmesh, true, triangulate);
  // if(triangulate) {
  //   const bool success = PMP::triangulate_faces(mesh);
  //   if(!success) {
  //     const std::string msg = "Triangulation has failed.";
  //     Rcpp::stop(msg);
  //   }
  // }
  if(!CGAL::is_closed(mesh)) {
    Rcpp::stop("The mesh is not closed.");
  }
  // std::string msg;
  // const bool bv = PMP::does_bound_a_volume(mesh);
  // if(bv) {
  //   msg = "The mesh bounds a volume.\n";
  // } else {
  //   msg2 = "The mesh does not bound a volume - reorienting.\n";
  //   PMP::orient_to_bound_a_volume(mesh);
  // }
  // Message(msg2);
  const K::FT vol = PMP::volume(mesh);
  return CGAL::to_double<K::FT>(vol);
}

// [[Rcpp::export]]
double meshAreaK(const Rcpp::List rmesh, const bool triangulate) {
  Mesh3 mesh = makeSurfMesh<Mesh3, Point3>(rmesh, true, triangulate);
  // if(triangulate) {
  //   const bool success = PMP::triangulate_faces(mesh);
  //   if(!success) {
  //     const std::string msg = "Triangulation has failed.";
  //     Rcpp::stop(msg);
  //   }
  // }
  const K::FT ar = PMP::area(mesh);
  return CGAL::to_double<K::FT>(ar);
}