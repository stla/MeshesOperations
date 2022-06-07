#ifndef _MESHESOPERATIONSHEADER_
#include "MeshesOperations.h"
#endif

// [[Rcpp::export]]
Rcpp::List SurfMesh(const Rcpp::List rmesh,
                    const bool isTriangle,
                    const bool triangulate,
                    const bool clean,
                    const bool normals,
                    const double epsilon) {
  Mesh3 mesh = makeSurfMesh<Mesh3, Point3>(rmesh, clean);
  const bool really_triangulate = !isTriangle && triangulate;
  Rcpp::IntegerMatrix Edges0;
  Rcpp::NumericMatrix Normals0;
  if(really_triangulate) {
    Edges0 = getEdges2<K, Mesh3, Point3>(mesh, epsilon);
    if(normals) {
      Normals0 = getKNormals(mesh);
    }
    bool success = CGAL::Polygon_mesh_processing::triangulate_faces(mesh);
    if(!success) {
      Rcpp::stop("Triangulation has failed.");
    }
  }
  Rcpp::List routmesh = RSurfKMesh(mesh, normals, epsilon);
  if(really_triangulate) {
    routmesh["edges0"] = Edges0;
    if(normals) {
      routmesh["normals0"] = Normals0;
    }
  }
  return routmesh;
}

// [[Rcpp::export]]
Rcpp::List SurfEMesh(const Rcpp::List rmesh,
                     const bool isTriangle,
                     const bool triangulate,
                     const bool clean,
                     const bool normals,
                     const double epsilon) {
  EMesh3 mesh = makeSurfMesh<EMesh3, EPoint3>(rmesh, clean);
  const bool really_triangulate = !isTriangle && triangulate;
  Rcpp::IntegerMatrix Edges0;
  Rcpp::NumericMatrix Normals0;
  if(really_triangulate) {
    Edges0 = getEdges2<EK, EMesh3, EPoint3>(mesh, epsilon);
    if(normals) {
      Normals0 = getEKNormals(mesh);
    }
    bool success = CGAL::Polygon_mesh_processing::triangulate_faces(mesh);
    if(!success) {
      Rcpp::stop("Triangulation has failed.");
    }
  }
  Rcpp::List routmesh = RSurfEKMesh(mesh, normals, epsilon);
  if(really_triangulate) {
    routmesh["edges0"] = Edges0;
    if(normals) {
      routmesh["normals0"] = Normals0;
    }
  }
  return routmesh;
}

// [[Rcpp::export]]
Rcpp::List SurfQMesh(const Rcpp::List rmesh,
                     const bool isTriangle,
                     const bool triangulate,
                     const bool clean,
                     const bool normals,
                     const double epsilon) {
  QMesh3 mesh = makeSurfQMesh(rmesh, clean);
  const bool really_triangulate = !isTriangle && triangulate;
  Rcpp::IntegerMatrix Edges0;
  Rcpp::NumericMatrix Normals0;
  if(really_triangulate) {
    Edges0 = getEdges2<QK, QMesh3, QPoint3>(mesh, epsilon);
    if(normals) {
      Normals0 = getQNormals(mesh);
    }
    bool success = CGAL::Polygon_mesh_processing::triangulate_faces(mesh);
    if(!success) {
      Rcpp::stop("Triangulation has failed.");
    }
  }
  Rcpp::List routmesh = RSurfQMesh(mesh, normals, epsilon);
  if(really_triangulate) {
    routmesh["edges0"] = Edges0;
    if(normals) {
      routmesh["normals0"] = Normals0;
    }
  }
  return routmesh;
}
