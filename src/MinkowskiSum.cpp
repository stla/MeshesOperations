#ifndef _MESHESOPERATIONSHEADER_
#include "MeshesOperations.h"
#endif

// [[Rcpp::export]]
Rcpp::List MinkowskiSumEK(const Rcpp::List rmesh1,
                          const Rcpp::List rmesh2,
                          const bool triangulate,
                          const bool normals) {
  EMesh3 mesh1 = makeSurfMesh<EMesh3, EPoint3>(rmesh1, true);
  EMesh3 mesh2 = makeSurfMesh<EMesh3, EPoint3>(rmesh2, true);
  ENef3 nef1(mesh1);
  ENef3 nef2(mesh2);
  ENef3 nef = CGAL::minkowski_sum_3(nef1, nef2);
  Rcpp::IntegerMatrix Edges0;
  if(triangulate) {
    EMesh3 mesh0;
    CGAL::convert_nef_polyhedron_to_polygon_mesh(nef, mesh0);
    Edges0 = getEdges2<EK, EMesh3, EPoint3>(mesh0, 0);
  }
  EMesh3 mesh;
  CGAL::convert_nef_polyhedron_to_polygon_mesh(nef, mesh, triangulate);
  Rcpp::List rmesh = RSurfEKMesh(mesh, normals, 0);
  if(triangulate) {
    rmesh["edges0"] = Edges0;
  }
  return rmesh;
}