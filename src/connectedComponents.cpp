#ifndef _MESHESOPERATIONSHEADER_
#include "MeshesOperations.h"
#endif

// [[Rcpp::export]]
Rcpp::List connectedComponentsK(const Rcpp::List rmesh0,
                    const bool isTriangle,
                    const bool triangulate,
                    const bool clean,
                    const bool normals,
                    const double epsilon) {
  const Mesh3 mesh0 = makeSurfMesh<Mesh3, Point3>(rmesh0, clean);
  std::vector<Mesh3> cc_meshes;
  PMP::split_connected_components(mesh0, cc_meshes);
  const size_t ncc = cc_meshes.size();
  // message ncc
  const bool really_triangulate = !isTriangle && triangulate;
  Rcpp::List cc_rmeshes(ncc);
  size_t i = 0;
  for(auto mesh = cc_meshes.begin(); mesh != cc_meshes.end(); ++mesh){
    Rcpp::IntegerMatrix Edges0;
    Rcpp::NumericMatrix Normals0;
    if(really_triangulate) {
      Edges0 = getEdges2<K, Mesh3, Point3>(*mesh, epsilon);
      if(normals) {
        Normals0 = getKNormals(*mesh);
      }
      const bool success = PMP::triangulate_faces(*mesh);
      if(!success) {
        Rcpp::stop("Triangulation has failed.");
      }
    }
    Rcpp::List cc_rmesh_i = RSurfKMesh(*mesh, normals, epsilon);
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
