#ifndef _MESHESOPERATIONSHEADER_
#include "MeshesOperations.h"
#endif

Mesh3 Triangulation(const Rcpp::List rmesh) {
  EMesh3 emesh = makeSurfMesh<EMesh3, EPoint3>(rmesh, true);
  Message("Triangulation.");
  const bool success = PMP::triangulate_faces(emesh);
  if(!success) {
    Rcpp::stop("Triangulation has failed.");
  }
  const size_t nvertices = emesh.number_of_vertices();
  const size_t nedges    = emesh.number_of_edges();
  const size_t nfaces    = emesh.number_of_faces();
  Mesh3 mesh;
  mesh.reserve(nvertices, nedges, nfaces);
  for(EMesh3::Vertex_index vd : emesh.vertices()) {
    const EPoint3 vertex = emesh.point(vd);
    const double x = CGAL::to_double<EK::FT>(vertex.x());
    const double y = CGAL::to_double<EK::FT>(vertex.y());
    const double z = CGAL::to_double<EK::FT>(vertex.z());
    mesh.add_vertex(Point3(x, y ,z));
  }
  for(EMesh3::Face_index fd : emesh.faces()) {
    std::vector<int> face;
    for(EMesh3::Vertex_index vd :
        vertices_around_face(emesh.halfedge(fd), emesh)) {
      face.push_back(vd);
    }
    mesh.add_face(
      CGAL::SM_Vertex_index(face[0]), 
      CGAL::SM_Vertex_index(face[1]), 
      CGAL::SM_Vertex_index(face[2])
    );
  }
  return mesh;
}

// [[Rcpp::export]]
Rcpp::NumericVector distanceK(const Rcpp::List rmesh,
                              const Rcpp::NumericMatrix points,
                              const bool triangulate) {
  Message("\u2014 Processing mesh...");
  Mesh3 mesh;
  if(triangulate) {
    mesh = Triangulation(rmesh);
  } else {
    mesh = makeSurfMesh<Mesh3, Point3>(rmesh, true);
  }
  Message("... done.\n");
  const size_t npoints = points.ncol();
  Rcpp::NumericVector distances(npoints);
  for(size_t i = 0; i < npoints; i++){
    Rcpp::NumericVector point_i = points(Rcpp::_, i);
    std::vector<Point3> pt = {Point3(point_i(0), point_i(1), point_i(2))};
    distances(i) = PMP::max_distance_to_triangle_mesh<CGAL::Sequential_tag>(
      pt, mesh
    );
  }
  return distances;
}
