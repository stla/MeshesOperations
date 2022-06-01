#ifndef _MESHESOPERATIONSHEADER_
#define _MESHESOPERATIONSHEADER_
#endif

// #define BOOST_THREAD_DONT_USE_MOVE
// #undef BOOST_GCC
// #pragma GCC diagnostic ignored "-Wclass-memaccess"
// [[Rcpp::plugins(cpp14)]]
#include <RcppEigen.h>
#define CGAL_EIGEN3_ENABLED 1

// // Helper macros to disable macros
// #if defined(__clang__) || defined(BOOST_GCC)
// #  define CGAL_PRAGMA_DIAG_PUSH _Pragma("GCC diagnostic push")
// #  define CGAL_PRAGMA_DIAG_POP  _Pragma("GCC diagnostic pop")
// #else
// #  define CGAL_PRAGMA_DIAG_PUSH
// #  define CGAL_PRAGMA_DIAG_POP
// #endif

// #include <CGAL/assertions.h>
// #undef CGAL_error
// #define CGAL_error
// #undef CGAL_error_msg
// #define CGAL_error_msg(msg)

// #include <iostream>

// #include <boost/noncopyable.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Vector_3.h>

#include <CGAL/utility.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>

#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <CGAL/Polygon_mesh_processing/corefinement.h>

#include <CGAL/number_utils.h>

#include <CGAL/Cartesian.h>
//#include <boost/multiprecision/gmp.hpp>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <CGAL/Gmpq.h>
#include <CGAL/Gmpz.h>
#include "gmp.h"

#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <CGAL/Polygon_mesh_processing/smooth_mesh.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>

#include <CGAL/Polygon_mesh_processing/smooth_shape.h>

#include <CGAL/Polygon_mesh_processing/remesh.h>

// to remove degenerate faces:
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>

// for volume:
#include <CGAL/Polygon_mesh_processing/measure.h>

// for sampling:
#include <CGAL/Polygon_mesh_processing/distance.h>

#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>
#include <CGAL/minkowski_sum_3.h>

// #include <algorithm>
// #include <array>
// #include <map>
// #include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_predicates_exact_constructions_kernel EK;

typedef K::Point_3 Point3;
typedef EK::Point_3 EPoint3;
typedef std::vector<Point3> Points3;

typedef CGAL::Surface_mesh<Point3> Mesh3;
typedef CGAL::Surface_mesh<EPoint3> EMesh3;

// typedef boost::graph_traits<Mesh3>::vertex_descriptor
// boost_vertex_descriptor; typedef boost::graph_traits<Mesh3>::face_descriptor
// boost_face_descriptor;

typedef K::Vector_3 Vector3;
typedef EK::Vector_3 EVector3;

typedef CGAL::Cartesian<CGAL::Gmpq> QK;
typedef CGAL::Surface_mesh<QK::Point_3> QMesh3;
typedef QK::Point_3 QPoint3;
typedef QK::Vector_3 QVector3;

typedef CGAL::Nef_polyhedron_3<EK> ENef3;

// -------------------------------------------------------------------------- //
namespace PMP = CGAL::Polygon_mesh_processing;
// namespace mp = boost::multiprecision;

// -------------------------------------------------------------------------- //

Rcpp::NumericMatrix points3_to_matrix(std::vector<Point3>);

template <typename PointT>
std::vector<PointT> matrix_to_points3(const Rcpp::NumericMatrix);

template <typename MeshT, typename PointT>
MeshT makeSurfMesh(const Rcpp::List, const bool);

template <typename MeshT, typename PointT>
MeshT makeSurfTMesh(const Rcpp::List, const bool);

QMesh3 makeSurfQMesh(const Rcpp::List, const bool);
QMesh3 makeSurfTQMesh(const Rcpp::List, const bool);

template <typename KernelT, typename MeshT, typename PointT>
Rcpp::IntegerMatrix getEdges2(MeshT, const double);

Rcpp::NumericMatrix getKNormals(Mesh3);
Rcpp::NumericMatrix getEKNormals(EMesh3);
Rcpp::NumericMatrix getQNormals(QMesh3);

Rcpp::List RSurfKMesh(Mesh3, const bool, const double);
Rcpp::List RSurfEKMesh(EMesh3, const bool, const double);
Rcpp::List RSurfQMesh(QMesh3, const bool, const double);
Rcpp::List RSurfTKMesh(Mesh3, const bool, const double);
Rcpp::List RSurfTEKMesh(EMesh3, const bool, const double);
Rcpp::List RSurfTQMesh(QMesh3, const bool, const double);

std::string q2str(CGAL::Gmpq);

void Message(std::string);

template <typename MeshT>
MeshT removeDegenerateFaces(MeshT);