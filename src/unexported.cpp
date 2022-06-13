#ifndef _MESHESOPERATIONSHEADER_
#include "MeshesOperations.h"
#endif

void Message(std::string msg) {
  SEXP rmsg = Rcpp::wrap(msg);
  Rcpp::message(rmsg);
}

std::string q2str(CGAL::Gmpq r) {
  CGAL::Gmpz numer = r.numerator();
  CGAL::Gmpz denom = r.denominator();
  size_t n = mpz_sizeinbase(numer.mpz(), 10) + 2;
  size_t d = mpz_sizeinbase(denom.mpz(), 10) + 2;
  char* cnumer = new char[n];
  char* cdenom = new char[d];
  cnumer = mpz_get_str(cnumer, 10, numer.mpz());
  cdenom = mpz_get_str(cdenom, 10, denom.mpz());
  std::string snumer = cnumer;
  std::string sdenom = cdenom;
  delete[] cnumer;
  delete[] cdenom;
  return snumer + "/" + sdenom;
}

Rcpp::NumericMatrix points3_to_matrix(std::vector<Point3> points) {
  const size_t nvertices = points.size();
  Rcpp::NumericMatrix Vertices(3, nvertices);
  for(size_t i = 0; i != nvertices; i++) {
    Rcpp::NumericVector col_i(3);
    const Point3 vertex = points[i];
    col_i(0) = vertex.x();
    col_i(1) = vertex.y();
    col_i(2) = vertex.z();
    Vertices(Rcpp::_, i) = col_i;
  }
  return Vertices;
}

template <typename PointT>
std::vector<PointT> matrix_to_points3(const Rcpp::NumericMatrix M) {
  const size_t npoints = M.ncol();
  std::vector<PointT> points;
  points.reserve(npoints);
  for(size_t i = 0; i < npoints; i++) {
    const Rcpp::NumericVector pt = M(Rcpp::_, i);
    points.emplace_back(PointT(pt(0), pt(1), pt(2)));
  }
  return points;
}

template std::vector<Point3> matrix_to_points3<Point3>(
    const Rcpp::NumericMatrix);
template std::vector<EPoint3> matrix_to_points3<EPoint3>(
    const Rcpp::NumericMatrix);

std::vector<QPoint3> matrix_to_qpoints3(const Rcpp::CharacterMatrix M) {
  const size_t npoints = M.ncol();
  std::vector<QPoint3> points;
  points.reserve(npoints);
  for(size_t i = 0; i < npoints; i++) {
    const Rcpp::CharacterVector pt = M(Rcpp::_, i);
    CGAL::Gmpq qpt0(CGAL::Gmpq(Rcpp::as<std::string>(pt(0))));
    CGAL::Gmpq qpt1(CGAL::Gmpq(Rcpp::as<std::string>(pt(1))));
    CGAL::Gmpq qpt2(CGAL::Gmpq(Rcpp::as<std::string>(pt(2))));
    points.emplace_back(QPoint3(qpt0, qpt1, qpt2));
  }
  return points;
}

std::vector<std::vector<int>> matrix_to_Tfaces(
    const Rcpp::IntegerMatrix Faces) {
  const size_t nfaces = Faces.ncol();
  std::vector<std::vector<int>> faces;
  faces.reserve(nfaces);
  for(size_t i = 0; i < nfaces; i++) {
    const Rcpp::IntegerVector face_rcpp = Faces(Rcpp::_, i);
    std::vector<int> face = {face_rcpp(0), face_rcpp(1), face_rcpp(2)};
    faces.emplace_back(face);
  }
  return faces;
}

std::vector<std::vector<int>> list_to_faces(const Rcpp::List L) {
  const size_t nfaces = L.size();
  std::vector<std::vector<int>> faces;
  faces.reserve(nfaces);
  for(size_t i = 0; i < nfaces; i++) {
    Rcpp::IntegerVector face_rcpp = Rcpp::as<Rcpp::IntegerVector>(L(i));
    std::vector<int> face(face_rcpp.begin(), face_rcpp.end());
    faces.emplace_back(face);
  }
  return faces;
}

template <typename MeshT, typename PointT>
MeshT soup2mesh(std::vector<PointT> points,
                std::vector<std::vector<int>> faces,
                const bool clean,
                const bool triangulate) {
  bool success = PMP::orient_polygon_soup(points, faces);
  if(success) {
    Message("Successful polygon orientation.");
  } else {
    Message("Polygon orientation failed.");
  }
  if(clean) {
    PMP::repair_polygon_soup(points, faces);
  }
  MeshT mesh;
  PMP::polygon_soup_to_polygon_mesh(points, faces, mesh);
  bool isTriangle;
  if(triangulate) {
    Message("Triangulation.");
    const bool success = PMP::triangulate_faces(mesh);
    if(!success) {
      Rcpp::stop("Triangulation has failed.");
    }
    isTriangle = true;
  } else {
    isTriangle = CGAL::is_triangle_mesh(mesh);
  } 
  if(isTriangle) {
    Message("The mesh is triangle.");
  } else {
    Message(
        "The mesh is not triangle; no way to ensure it bounds a volume "
        "and whether it is outward oriented.");
  }
  if(CGAL::is_closed(mesh)) {
    Message("The mesh is closed.");
    if(isTriangle) {
      if(!PMP::is_outward_oriented(mesh)) {
        PMP::reverse_face_orientations(mesh);
      }
      const bool bv = PMP::does_bound_a_volume(mesh);
      std::string msg2;
      if(bv) {
        msg2 = "The mesh bounds a volume.";
      } else {
        msg2 = "The mesh does not bound a volume - reorienting.";
        PMP::orient_to_bound_a_volume(mesh);
      }
      Message(msg2);
    }
  } else {
    Message("The mesh is not closed.");
  }
  return mesh;
}

template <typename MeshT, typename PointT>
MeshT makeSurfMesh(
  const Rcpp::List rmesh, const bool clean, const bool triangulate
) {
  const Rcpp::NumericMatrix vertices =
      Rcpp::as<Rcpp::NumericMatrix>(rmesh["vertices"]);
  const Rcpp::List rfaces = Rcpp::as<Rcpp::List>(rmesh["faces"]);
  std::vector<PointT> points = matrix_to_points3<PointT>(vertices);
  std::vector<std::vector<int>> faces = list_to_faces(rfaces);
  return soup2mesh<MeshT, PointT>(points, faces, clean, triangulate);
}

template Mesh3 makeSurfMesh<Mesh3, Point3>(
  const Rcpp::List, const bool, const bool
);
template EMesh3 makeSurfMesh<EMesh3, EPoint3>(
  const Rcpp::List, const bool, const bool
);

template <typename MeshT, typename PointT>
MeshT makeSurfTMesh(
  const Rcpp::List rmesh, const bool clean, const bool triangulate
) {
  const Rcpp::NumericMatrix vertices =
      Rcpp::as<Rcpp::NumericMatrix>(rmesh["vertices"]);
  const Rcpp::IntegerMatrix rfaces =
      Rcpp::as<Rcpp::IntegerMatrix>(rmesh["faces"]);
  std::vector<PointT> points = matrix_to_points3<PointT>(vertices);
  std::vector<std::vector<int>> faces = matrix_to_Tfaces(rfaces);
  return soup2mesh<MeshT, PointT>(points, faces, clean, triangulate);
}

template Mesh3 makeSurfTMesh<Mesh3, Point3>(
  const Rcpp::List, const bool, const bool
);
template EMesh3 makeSurfTMesh<EMesh3, EPoint3>(
  const Rcpp::List, const bool, const bool
);

QMesh3 makeSurfQMesh(
  const Rcpp::List rmesh, const bool clean, const bool triangulate
) {
  const Rcpp::CharacterMatrix vertices =
      Rcpp::as<Rcpp::CharacterMatrix>(rmesh["vertices"]);
  const Rcpp::List rfaces = Rcpp::as<Rcpp::List>(rmesh["faces"]);
  std::vector<QPoint3> points = matrix_to_qpoints3(vertices);
  std::vector<std::vector<int>> faces = list_to_faces(rfaces);
  return soup2mesh<QMesh3, QPoint3>(points, faces, clean, triangulate);
}

QMesh3 makeSurfTQMesh(
  const Rcpp::List rmesh, const bool clean, const bool triangulate
) {
  const Rcpp::CharacterMatrix vertices =
      Rcpp::as<Rcpp::CharacterMatrix>(rmesh["vertices"]);
  const Rcpp::IntegerMatrix rfaces =
      Rcpp::as<Rcpp::IntegerMatrix>(rmesh["faces"]);
  std::vector<QPoint3> points = matrix_to_qpoints3(vertices);
  std::vector<std::vector<int>> faces = matrix_to_Tfaces(rfaces);
  return soup2mesh<QMesh3, QPoint3>(points, faces, clean, triangulate);
}

Rcpp::NumericMatrix getVertices_K(Mesh3 mesh) {
  const size_t nvertices = mesh.number_of_vertices();
  Rcpp::NumericMatrix Vertices(3, nvertices);
  {
    size_t i = 0;
    for(Mesh3::Vertex_index vd : mesh.vertices()) {
      Rcpp::NumericVector col_i(3);
      const Point3 vertex = mesh.point(vd);
      col_i(0) = vertex.x();
      col_i(1) = vertex.y();
      col_i(2) = vertex.z();
      Vertices(Rcpp::_, i) = col_i;
      i++;
    }
  }
  return Vertices;
}

Rcpp::NumericMatrix getVertices_EK(EMesh3 mesh) {
  const size_t nvertices = mesh.number_of_vertices();
  Rcpp::NumericMatrix Vertices(3, nvertices);
  {
    size_t i = 0;
    for(EMesh3::Vertex_index vd : mesh.vertices()) {
      Rcpp::NumericVector col_i(3);
      const EPoint3 vertex = mesh.point(vd);
      col_i(0) = CGAL::to_double<EK::FT>(vertex.x());
      col_i(1) = CGAL::to_double<EK::FT>(vertex.y());
      col_i(2) = CGAL::to_double<EK::FT>(vertex.z());
      Vertices(Rcpp::_, i) = col_i;
      i++;
    }
  }
  return Vertices;
}

Rcpp::CharacterMatrix getVertices_QK(QMesh3 mesh) {
  const size_t nvertices = mesh.number_of_vertices();
  Rcpp::CharacterMatrix Vertices(3, nvertices);
  {
    size_t i = 0;
    for(QMesh3::Vertex_index vd : mesh.vertices()) {
      Rcpp::CharacterVector col_i(3);
      const QPoint3 vertex = mesh.point(vd);
      col_i(0) = q2str(vertex.x());
      col_i(1) = q2str(vertex.y());
      col_i(2) = q2str(vertex.z());
      Vertices(Rcpp::_, i) = col_i;
      i++;
    }
  }
  return Vertices;
}

template <typename KernelT, typename MeshT, typename PointT>
Rcpp::DataFrame getEdges2(MeshT mesh) {
  const size_t nedges = mesh.number_of_edges();
  Rcpp::IntegerVector I1(nedges);
  Rcpp::IntegerVector I2(nedges);
  Rcpp::NumericVector Length(nedges);
  Rcpp::NumericVector Angle(nedges);
  Rcpp::LogicalVector Exterior(nedges);
  Rcpp::LogicalVector Coplanar(nedges);
  {
    size_t i = 0;
    for(typename MeshT::Edge_index ed : mesh.edges()) {
      typename MeshT::Vertex_index s = source(ed, mesh);
      typename MeshT::Vertex_index t = target(ed, mesh);
      I1(i) = (int)s + 1;
      I2(i) = (int)t + 1;
      std::vector<PointT> points(4);
      points[0] = mesh.point(s);
      points[1] = mesh.point(t);
      typename MeshT::Halfedge_index h0 = mesh.halfedge(ed, 0);
      points[2] = mesh.point(mesh.target(mesh.next(h0)));
      typename MeshT::Halfedge_index h1 = mesh.halfedge(ed, 1);
      points[3] = mesh.point(mesh.target(mesh.next(h1)));
      typename KernelT::FT angle = CGAL::abs(CGAL::approximate_dihedral_angle(
          points[0], points[1], points[2], points[3]));
      Angle(i) = CGAL::to_double(angle);
      Exterior(i) = angle < 179.0 || angle > 181.0;
      Coplanar(i) = CGAL::coplanar(points[0], points[1], points[2], points[3]);
      typename KernelT::FT el = PMP::edge_length(h0, mesh);
      Length(i) = CGAL::to_double(el);
      i++;
    }
  }
  Rcpp::DataFrame Edges = Rcpp::DataFrame::create(
    Rcpp::Named("i1")       = I1,
    Rcpp::Named("i2")       = I2,
    Rcpp::Named("length")   = Length,
    Rcpp::Named("angle")    = Angle,
    Rcpp::Named("exterior") = Exterior,
    Rcpp::Named("coplanar") = Coplanar
  );
  return Edges;
}

template Rcpp::DataFrame getEdges2<K, Mesh3, Point3>(Mesh3);
template Rcpp::DataFrame getEdges2<EK, EMesh3, EPoint3>(EMesh3);
template Rcpp::DataFrame getEdges2<QK, QMesh3, QPoint3>(QMesh3);

template <typename MeshT>
Rcpp::List getFaces(MeshT mesh) {
  const size_t nfaces = mesh.number_of_faces();
  Rcpp::List Faces(nfaces);
  {
    size_t i = 0;
    for(typename MeshT::Face_index fd : mesh.faces()) {
      Rcpp::IntegerVector col_i;
      for(typename MeshT::Vertex_index vd :
          vertices_around_face(mesh.halfedge(fd), mesh)) {
        col_i.push_back(vd + 1);
      }
      Faces(i) = col_i;
      i++;
    }
  }
  return Faces;
}

template <typename MeshT>
Rcpp::IntegerMatrix getTFaces(MeshT mesh) {
  const size_t nfaces = mesh.number_of_faces();
  Rcpp::IntegerMatrix Faces(3, nfaces);
  {
    size_t i = 0;
    for(typename MeshT::Face_index fd : mesh.faces()) {
      Rcpp::IntegerVector col_i;
      for(typename MeshT::Vertex_index vd :
          vertices_around_face(mesh.halfedge(fd), mesh)) {
        col_i.push_back(vd + 1);
      }
      Faces(Rcpp::_, i) = col_i;
      i++;
    }
  }
  return Faces;
}

Rcpp::NumericMatrix getKNormals(Mesh3 mesh) {
  const size_t nvertices = mesh.number_of_vertices();
  Rcpp::NumericMatrix Normals(3, nvertices);
  auto vnormals = mesh.add_property_map<Mesh3::Vertex_index, Vector3>(
                          "v:normals", CGAL::NULL_VECTOR)
                      .first;
  auto fnormals = mesh.add_property_map<Mesh3::Face_index, Vector3>(
                          "f:normals", CGAL::NULL_VECTOR)
                      .first;
  PMP::compute_normals(mesh, vnormals, fnormals);
  {
    size_t i = 0;
    for(Mesh3::Vertex_index vd : vertices(mesh)) {
      Rcpp::NumericVector col_i(3);
      const Vector3 normal = vnormals[vd];
      col_i(0) = normal.x();
      col_i(1) = normal.y();
      col_i(2) = normal.z();
      Normals(Rcpp::_, i) = col_i;
      i++;
    }
  }
  return Normals;
}

Rcpp::NumericMatrix getEKNormals(EMesh3 mesh) {
  const size_t nvertices = mesh.number_of_vertices();
  Rcpp::NumericMatrix Normals(3, nvertices);
  auto vnormals = mesh.add_property_map<EMesh3::Vertex_index, EVector3>(
                          "v:normals", CGAL::NULL_VECTOR)
                      .first;
  auto fnormals = mesh.add_property_map<EMesh3::Face_index, EVector3>(
                          "f:normals", CGAL::NULL_VECTOR)
                      .first;
  PMP::compute_normals(mesh, vnormals, fnormals);
  {
    size_t i = 0;
    for(EMesh3::Vertex_index vd : vertices(mesh)) {
      Rcpp::NumericVector col_i(3);
      const EVector3 normal = vnormals[vd];
      col_i(0) = CGAL::to_double<EK::FT>(normal.x());
      col_i(1) = CGAL::to_double<EK::FT>(normal.y());
      col_i(2) = CGAL::to_double<EK::FT>(normal.z());
      Normals(Rcpp::_, i) = col_i;
      i++;
    }
  }
  return Normals;
}

Rcpp::NumericMatrix getQNormals(QMesh3 mesh) {
  const size_t nvertices = mesh.number_of_vertices();
  Rcpp::NumericMatrix Normals(3, nvertices);
  auto vnormals = mesh.add_property_map<QMesh3::Vertex_index, QVector3>(
                          "v:normals", CGAL::NULL_VECTOR)
                      .first;
  auto fnormals = mesh.add_property_map<QMesh3::Face_index, QVector3>(
                          "f:normals", CGAL::NULL_VECTOR)
                      .first;
  PMP::compute_normals(mesh, vnormals, fnormals);
  {
    size_t i = 0;
    for(QMesh3::Vertex_index vd : vertices(mesh)) {
      Rcpp::NumericVector col_i(3);
      const QVector3 normal = vnormals[vd];
      col_i(0) = normal.x().to_double();
      col_i(1) = normal.y().to_double();
      col_i(2) = normal.z().to_double();
      Normals(Rcpp::_, i) = col_i;
      i++;
    }
  }
  return Normals;
}

Rcpp::List RSurfKMesh(Mesh3 mesh, const bool normals) {
  Rcpp::DataFrame Edges = getEdges2<K, Mesh3, Point3>(mesh);
  Rcpp::NumericMatrix Vertices = getVertices_K(mesh);
  Rcpp::List Faces = getFaces<Mesh3>(mesh);
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("vertices") = Vertices,
                                      Rcpp::Named("edges") = Edges,
                                      Rcpp::Named("faces") = Faces);
  if(normals) {
    Rcpp::NumericMatrix Normals = getKNormals(mesh);
    out["normals"] = Normals;
  }
  return out;
}

Rcpp::List RSurfEKMesh(EMesh3 mesh, const bool normals) {
  Rcpp::DataFrame Edges = getEdges2<EK, EMesh3, EPoint3>(mesh);
  Rcpp::NumericMatrix Vertices = getVertices_EK(mesh);
  Rcpp::List Faces = getFaces<EMesh3>(mesh);
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("vertices") = Vertices,
                                      Rcpp::Named("edges") = Edges,
                                      Rcpp::Named("faces") = Faces);
  if(normals) {
    Rcpp::NumericMatrix Normals = getEKNormals(mesh);
    out["normals"] = Normals;
  }
  return out;
}

Rcpp::List RSurfQMesh(QMesh3 mesh, const bool normals) {
  Rcpp::DataFrame Edges = getEdges2<QK, QMesh3, QPoint3>(mesh);
  Rcpp::CharacterMatrix Vertices = getVertices_QK(mesh);
  Rcpp::List Faces = getFaces<QMesh3>(mesh);
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("vertices") = Vertices,
                                      Rcpp::Named("edges") = Edges,
                                      Rcpp::Named("faces") = Faces);
  if(normals) {
    Rcpp::NumericMatrix Normals = getQNormals(mesh);
    out["normals"] = Normals;
  }
  return out;
}

Rcpp::List RSurfTKMesh(Mesh3 mesh, const bool normals) {
  Rcpp::DataFrame Edges = getEdges2<K, Mesh3, Point3>(mesh);
  Rcpp::NumericMatrix Vertices = getVertices_K(mesh);
  Rcpp::IntegerMatrix Faces = getTFaces<Mesh3>(mesh);
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("vertices") = Vertices,
                                      Rcpp::Named("edges") = Edges,
                                      Rcpp::Named("faces") = Faces);
  if(normals) {
    Rcpp::NumericMatrix Normals = getKNormals(mesh);
    out["normals"] = Normals;
  }
  return out;
}

Rcpp::List RSurfTEKMesh(EMesh3 mesh, const bool normals) {
  Rcpp::DataFrame Edges = getEdges2<EK, EMesh3, EPoint3>(mesh);
  Rcpp::NumericMatrix Vertices = getVertices_EK(mesh);
  Rcpp::IntegerMatrix Faces = getTFaces<EMesh3>(mesh);
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("vertices") = Vertices,
                                      Rcpp::Named("edges") = Edges,
                                      Rcpp::Named("faces") = Faces);
  if(normals) {
    Rcpp::NumericMatrix Normals = getEKNormals(mesh);
    out["normals"] = Normals;
  }
  return out;
}

Rcpp::List RSurfTQMesh(QMesh3 mesh, const bool normals) {
  Rcpp::DataFrame Edges = getEdges2<QK, QMesh3, QPoint3>(mesh);
  Rcpp::CharacterMatrix Vertices = getVertices_QK(mesh);
  Rcpp::IntegerMatrix Faces = getTFaces<QMesh3>(mesh);
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("vertices") = Vertices,
                                      Rcpp::Named("edges") = Edges,
                                      Rcpp::Named("faces") = Faces);
  if(normals) {
    Rcpp::NumericMatrix Normals = getQNormals(mesh);
    out["normals"] = Normals;
  }
  return out;
}

template <typename MeshT>
MeshT removeDegenerateFaces(MeshT mesh) {  // triangular faces only
  for(typename MeshT::Face_index fd : mesh.faces()) {
    if(PMP::is_degenerate_triangle_face(fd, mesh)) {
      mesh.remove_face(fd);
    }
  }
  const size_t nrmfs = mesh.number_of_removed_faces();
  if(nrmfs > 0) {
    const std::string endmsg =
        nrmfs == 1 ? "one degenerate face"
                   : (std::to_string(nrmfs) + "degenerate faces");
    const std::string msg = "Removed " + endmsg + ".\n";
    Message(msg);
    mesh.collect_garbage();
  }
  return mesh;
}

template Mesh3 removeDegenerateFaces<Mesh3>(Mesh3);

//{}
