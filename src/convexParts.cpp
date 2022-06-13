#ifndef _MESHESOPERATIONSHEADER_
#include "MeshesOperations.h"
#endif

// [[Rcpp::export]]
Rcpp::List convexDecomposition(Rcpp::List rmesh, const bool triangulate){
  Message("\u2014 Processing mesh...");
  EMesh3 mesh = makeSurfMesh<EMesh3, EPoint3>(rmesh, true, true);
  Message("... done.\n");
  NefPol nef(mesh);
  CGAL::convex_decomposition_3(nef);
  std::list<EMesh3> convex_parts;
  // the first volume is the outer volume, which is
  // ignored in the decomposition
  NefPol::Volume_const_iterator ci = ++nef.volumes_begin();
  for( ; ci != nef.volumes_end(); ++ci) {
    if(ci->mark()) {
      EPolyhedron pol;
      nef.convert_inner_shell_to_polyhedron(ci->shells_begin(), pol);
      EMesh3 cmesh;
      CGAL::copy_face_graph(pol, cmesh);
      convex_parts.push_back(cmesh);
    }
  }
  const size_t ncp = convex_parts.size();
  std::string msg;
  if(ncp == 1) {
    msg = "Only one convex part.";
  } else {
    msg = "Found " + std::to_string(ncp) + " convex parts.";
  }
  Message(msg);
  Rcpp::List out(ncp);
  size_t i = 0;
  for(EMesh3 cmesh : convex_parts){
    if(triangulate && !CGAL::is_triangle_mesh(cmesh)){
      Message(
        "Triangulation of convex part n\u00b0" + std::to_string(i+1) + "."
      );
      if(!PMP::triangulate_faces(cmesh)) {
        Rcpp::stop("Triangulation has failed.");
      }
    }
    out(i) = RSurfEKMesh(cmesh, false);
    i++;
  }
  return out;
}

//{}