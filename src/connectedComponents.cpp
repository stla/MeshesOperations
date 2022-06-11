#ifndef _MESHESOPERATIONSHEADER_
#include "MeshesOperations.h"
#endif

// [[Rcpp::export]]
Rcpp::List connectedComponentsK(const Rcpp::List rmesh0,
                                const bool isTriangle,
                                const bool triangulate,
                                const bool clean,
                                const bool normals) {
  Message("\u2014 Processing mesh...");
  const Mesh3 mesh0 = makeSurfMesh<Mesh3, Point3>(rmesh0, clean);
  Message("... done.\n");
  std::vector<Mesh3> cc_meshes;
  PMP::split_connected_components(mesh0, cc_meshes);
  const size_t ncc = cc_meshes.size();
  if(ncc == 1) {
    Message("Only one component found.\n");
  } else {
    const std::string msg = "Found " + std::to_string(ncc) + " components.\n";
    Message(msg);
  }
  const bool really_triangulate = !isTriangle && triangulate;
  Rcpp::List cc_rmeshes(ncc);
  size_t i = 0;
  for(auto mesh = cc_meshes.begin(); mesh != cc_meshes.end(); ++mesh) {
    Rcpp::DataFrame Edges0;
    Rcpp::NumericMatrix Normals0;
    if(really_triangulate) {
      Edges0 = getEdges2<K, Mesh3, Point3>(*mesh);
      if(normals) {
        Normals0 = getKNormals(*mesh);
      }
      const bool success = PMP::triangulate_faces(*mesh);
      if(!success) {
        const std::string msg = "Triangulation has failed (component " +
                                std::to_string(i + 1) + ").";
        Rcpp::stop(msg);
      }
    }
    Rcpp::List cc_rmesh_i = RSurfKMesh(*mesh, normals);
    if(really_triangulate) {
      cc_rmesh_i["edges0"] = Edges0;
      if(normals) {
        cc_rmesh_i["normals0"] = Normals0;
      }
    }
    cc_rmeshes(i) = cc_rmesh_i;
    i++;
  }
  return cc_rmeshes;
}

// [[Rcpp::export]]
Rcpp::List connectedComponentsEK(const Rcpp::List rmesh0,
                                 const bool isTriangle,
                                 const bool triangulate,
                                 const bool clean,
                                 const bool normals) {
  Message("\u2014 Processing mesh...");
  const EMesh3 mesh0 = makeSurfMesh<EMesh3, EPoint3>(rmesh0, clean);
  Message("... done.\n");
  std::vector<EMesh3> cc_meshes;
  PMP::split_connected_components(mesh0, cc_meshes);
  const size_t ncc = cc_meshes.size();
  if(ncc == 1) {
    Message("Only one component found.\n");
  } else {
    const std::string msg = "Found " + std::to_string(ncc) + " components.\n";
    Message(msg);
  }
  const bool really_triangulate = !isTriangle && triangulate;
  Rcpp::List cc_rmeshes(ncc);
  size_t i = 0;
  for(auto mesh = cc_meshes.begin(); mesh != cc_meshes.end(); ++mesh) {
    Rcpp::DataFrame Edges0;
    Rcpp::NumericMatrix Normals0;
    if(really_triangulate) {
      Edges0 = getEdges2<EK, EMesh3, EPoint3>(*mesh);
      if(normals) {
        Normals0 = getEKNormals(*mesh);
      }
      const bool success = PMP::triangulate_faces(*mesh);
      if(!success) {
        const std::string msg = "Triangulation has failed (component " +
                                std::to_string(i + 1) + ").";
        Rcpp::stop(msg);
      }
    }
    Rcpp::List cc_rmesh_i = RSurfEKMesh(*mesh, normals);
    if(really_triangulate) {
      cc_rmesh_i["edges0"] = Edges0;
      if(normals) {
        cc_rmesh_i["normals0"] = Normals0;
      }
    }
    cc_rmeshes(i) = cc_rmesh_i;
    i++;
  }
  return cc_rmeshes;
}

// [[Rcpp::export]]
Rcpp::List connectedComponentsQ(const Rcpp::List rmesh0,
                                const bool isTriangle,
                                const bool triangulate,
                                const bool clean,
                                const bool normals) {
  Message("\u2014 Processing mesh...");
  const QMesh3 mesh0 = makeSurfQMesh(rmesh0, clean);
  Message("... done.\n");
  std::vector<QMesh3> cc_meshes;
  PMP::split_connected_components(mesh0, cc_meshes);
  const size_t ncc = cc_meshes.size();
  if(ncc == 1) {
    Message("Only one component found.\n");
  } else {
    const std::string msg = "Found " + std::to_string(ncc) + " components.\n";
    Message(msg);
  }
  const bool really_triangulate = !isTriangle && triangulate;
  Rcpp::List cc_rmeshes(ncc);
  size_t i = 0;
  for(auto mesh = cc_meshes.begin(); mesh != cc_meshes.end(); ++mesh) {
    Rcpp::DataFrame Edges0;
    Rcpp::NumericMatrix Normals0;
    if(really_triangulate) {
      Edges0 = getEdges2<QK, QMesh3, QPoint3>(*mesh);
      if(normals) {
        Normals0 = getQNormals(*mesh);
      }
      const bool success = PMP::triangulate_faces(*mesh);
      if(!success) {
        const std::string msg = "Triangulation has failed (component " +
                                std::to_string(i + 1) + ").";
        Rcpp::stop(msg);
      }
    }
    Rcpp::List cc_rmesh_i = RSurfQMesh(*mesh, normals);
    if(really_triangulate) {
      cc_rmesh_i["edges0"] = Edges0;
      if(normals) {
        cc_rmesh_i["normals0"] = Normals0;
      }
    }
    cc_rmeshes(i) = cc_rmesh_i;
    i++;
  }
  return cc_rmeshes;
}
