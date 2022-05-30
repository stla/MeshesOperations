#ifndef _MESHESOPERATIONSHEADER_
#include "MeshesOperations.h"
#endif

#include <CGAL/IO/io.h>
#include <locale>  // tolower

std::string toLower(std::string s) {
  for(char& c : s) {
    c = std::tolower(c);
  }
  return s;
}

// [[Rcpp::export]]
Rcpp::List readFile(const std::string filename) {
  const std::string ext = toLower(filename.substr(filename.length() - 3, 3));
  std::ifstream infile;
  infile.open(filename);
  const bool binary = CGAL::IO::is_binary(infile);
  std::vector<Point3> points;
  std::vector<std::vector<int>> faces;
  bool ok = false;
  if(ext == "ply") {
    ok = CGAL::IO::read_PLY(infile, points, faces,
                            CGAL::parameters::use_binary_mode(binary));
  } else if(ext == "stl") {
    ok = CGAL::IO::read_STL(filename, points, faces,
                            CGAL::parameters::use_binary_mode(binary));
  } else if(ext == "obj") {
    ok = CGAL::IO::read_OBJ(infile, points, faces);
  } else if(ext == "off") {
    ok = CGAL::IO::read_OFF(infile, points, faces);
  } else {
    Rcpp::stop("Unknown file extension.");
  }
  infile.close();
  Rcpp::List out;
  if(ok) {
    const size_t npoints = points.size();
    Rcpp::NumericMatrix Vertices(3, npoints);
    for(size_t i = 0; i < npoints; i++) {
      const Point3 point_i = points[i];
      Rcpp::NumericVector col_i =
          Rcpp::NumericVector::create(point_i.x(), point_i.y(), point_i.z());
      Vertices(Rcpp::_, i) = col_i;
    }
    const size_t nfaces = faces.size();
    Rcpp::List Faces(nfaces);
    for(size_t i = 0; i < nfaces; i++) {
      const std::vector<int> face_i = faces[i];
      Rcpp::IntegerVector col_i(face_i.begin(), face_i.end());
      Faces(i) = col_i + 1;
    }
    out["vertices"] = Rcpp::transpose(Vertices);
    out["faces"] = Faces;
  } else {
    Rcpp::stop("Reading failure.");
  }
  return out;
}

//{}
