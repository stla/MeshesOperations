#ifndef _MESHESOPERATIONSHEADER_
#include "MeshesOperations.h"
#endif

// [[Rcpp::export]]
Rcpp::List clipMeshEK(const Rcpp::List rmesh,
                      const Rcpp::List rclipper,
                      const bool clipVolume,
                      const bool triangulate1,
                      const bool triangulate2,
                      const bool normals) {
  Message("Processing mesh to be clipped.\n");
  EMesh3 mesh = makeSurfMesh<EMesh3, EPoint3>(rmesh, true);
  if(triangulate1) {
    const bool success = PMP::triangulate_faces(mesh);
    if(!success) {
      const std::string msg = "Triangulation has failed.";
      Rcpp::stop(msg);
    }
  }
  bool doNotModify = false;
  if(!clipVolume) {
    doNotModify = true;
  } else {
    if(PMP::does_self_intersect(mesh)) {
      Rcpp::stop("The mesh self-intersects.");
    }    
  } 
  Message("Processing the clipping mesh.\n");
  EMesh3 clipper = makeSurfMesh<EMesh3, EPoint3>(rclipper, true);
  if(triangulate2) {
    const bool success = PMP::triangulate_faces(clipper);
    if(!success) {
      const std::string msg = "Triangulation has failed.";
      Rcpp::stop(msg);
    }
  }
  if(PMP::does_self_intersect(clipper)) {
    Rcpp::stop("The clipping mesh self-intersects.");
  }    
  if(!PMP::does_bound_a_volume(clipper)) {
    Rcpp::stop("The clipping mesh does not bound a volume.");
  }    
  Message("Performing clipping.\n");
  const bool clipping = CGAL::Polygon_mesh_processing::clip(
    mesh, clipper, PMP::parameters::clip_volume(clipVolume),
    PMP::parameters::clip_volume(clipVolume).do_not_modify(doNotModify)
  );
  if(!clipping) {
    Rcpp::stop("Clipping has failed.");
  }
  mesh.collect_garbage();
  Rcpp::List routmesh = RSurfTEKMesh(mesh, normals, 0);
  Rcpp::List out;
  if(!doNotModify) {
    clipper.collect_garbage();
    Rcpp::List routclipper = RSurfTEKMesh(clipper, normals, 0);
    out["mesh"] = routmesh;
    out["clipper"] = routclipper;
  } else {
    out = routmesh;
  }
  return out;
}
