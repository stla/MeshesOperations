// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// MinkowskiSumEK
Rcpp::List MinkowskiSumEK(const Rcpp::List rmesh1, const Rcpp::List rmesh2, const bool triangulate, const bool normals, const bool triangulate1, const bool triangulate2);
RcppExport SEXP _MeshesOperations_MinkowskiSumEK(SEXP rmesh1SEXP, SEXP rmesh2SEXP, SEXP triangulateSEXP, SEXP normalsSEXP, SEXP triangulate1SEXP, SEXP triangulate2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh1(rmesh1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh2(rmesh2SEXP);
    Rcpp::traits::input_parameter< const bool >::type triangulate(triangulateSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    Rcpp::traits::input_parameter< const bool >::type triangulate1(triangulate1SEXP);
    Rcpp::traits::input_parameter< const bool >::type triangulate2(triangulate2SEXP);
    rcpp_result_gen = Rcpp::wrap(MinkowskiSumEK(rmesh1, rmesh2, triangulate, normals, triangulate1, triangulate2));
    return rcpp_result_gen;
END_RCPP
}
// readFile
Rcpp::List readFile(const std::string filename);
RcppExport SEXP _MeshesOperations_readFile(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(readFile(filename));
    return rcpp_result_gen;
END_RCPP
}
// writeFile
void writeFile(const std::string filename, const bool binary, const int precision, const Rcpp::NumericMatrix Vertices, const Rcpp::List Faces);
RcppExport SEXP _MeshesOperations_writeFile(SEXP filenameSEXP, SEXP binarySEXP, SEXP precisionSEXP, SEXP VerticesSEXP, SEXP FacesSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< const bool >::type binary(binarySEXP);
    Rcpp::traits::input_parameter< const int >::type precision(precisionSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix >::type Vertices(VerticesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type Faces(FacesSEXP);
    writeFile(filename, binary, precision, Vertices, Faces);
    return R_NilValue;
END_RCPP
}
// clipMeshEK
Rcpp::List clipMeshEK(const Rcpp::List rmesh, const Rcpp::List rclipper, const bool clipVolume, const bool triangulate1, const bool triangulate2, const bool normals);
RcppExport SEXP _MeshesOperations_clipMeshEK(SEXP rmeshSEXP, SEXP rclipperSEXP, SEXP clipVolumeSEXP, SEXP triangulate1SEXP, SEXP triangulate2SEXP, SEXP normalsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh(rmeshSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type rclipper(rclipperSEXP);
    Rcpp::traits::input_parameter< const bool >::type clipVolume(clipVolumeSEXP);
    Rcpp::traits::input_parameter< const bool >::type triangulate1(triangulate1SEXP);
    Rcpp::traits::input_parameter< const bool >::type triangulate2(triangulate2SEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    rcpp_result_gen = Rcpp::wrap(clipMeshEK(rmesh, rclipper, clipVolume, triangulate1, triangulate2, normals));
    return rcpp_result_gen;
END_RCPP
}
// connectedComponentsK
Rcpp::List connectedComponentsK(const Rcpp::List rmesh0, const bool isTriangle, const bool triangulate, const bool clean, const bool normals, const double epsilon);
RcppExport SEXP _MeshesOperations_connectedComponentsK(SEXP rmesh0SEXP, SEXP isTriangleSEXP, SEXP triangulateSEXP, SEXP cleanSEXP, SEXP normalsSEXP, SEXP epsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh0(rmesh0SEXP);
    Rcpp::traits::input_parameter< const bool >::type isTriangle(isTriangleSEXP);
    Rcpp::traits::input_parameter< const bool >::type triangulate(triangulateSEXP);
    Rcpp::traits::input_parameter< const bool >::type clean(cleanSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    Rcpp::traits::input_parameter< const double >::type epsilon(epsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(connectedComponentsK(rmesh0, isTriangle, triangulate, clean, normals, epsilon));
    return rcpp_result_gen;
END_RCPP
}
// connectedComponentsEK
Rcpp::List connectedComponentsEK(const Rcpp::List rmesh0, const bool isTriangle, const bool triangulate, const bool clean, const bool normals, const double epsilon);
RcppExport SEXP _MeshesOperations_connectedComponentsEK(SEXP rmesh0SEXP, SEXP isTriangleSEXP, SEXP triangulateSEXP, SEXP cleanSEXP, SEXP normalsSEXP, SEXP epsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh0(rmesh0SEXP);
    Rcpp::traits::input_parameter< const bool >::type isTriangle(isTriangleSEXP);
    Rcpp::traits::input_parameter< const bool >::type triangulate(triangulateSEXP);
    Rcpp::traits::input_parameter< const bool >::type clean(cleanSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    Rcpp::traits::input_parameter< const double >::type epsilon(epsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(connectedComponentsEK(rmesh0, isTriangle, triangulate, clean, normals, epsilon));
    return rcpp_result_gen;
END_RCPP
}
// connectedComponentsQ
Rcpp::List connectedComponentsQ(const Rcpp::List rmesh0, const bool isTriangle, const bool triangulate, const bool clean, const bool normals, const double epsilon);
RcppExport SEXP _MeshesOperations_connectedComponentsQ(SEXP rmesh0SEXP, SEXP isTriangleSEXP, SEXP triangulateSEXP, SEXP cleanSEXP, SEXP normalsSEXP, SEXP epsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh0(rmesh0SEXP);
    Rcpp::traits::input_parameter< const bool >::type isTriangle(isTriangleSEXP);
    Rcpp::traits::input_parameter< const bool >::type triangulate(triangulateSEXP);
    Rcpp::traits::input_parameter< const bool >::type clean(cleanSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    Rcpp::traits::input_parameter< const double >::type epsilon(epsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(connectedComponentsQ(rmesh0, isTriangle, triangulate, clean, normals, epsilon));
    return rcpp_result_gen;
END_RCPP
}
// meshVolumeK
double meshVolumeK(const Rcpp::List rmesh, const bool triangulate);
RcppExport SEXP _MeshesOperations_meshVolumeK(SEXP rmeshSEXP, SEXP triangulateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh(rmeshSEXP);
    Rcpp::traits::input_parameter< const bool >::type triangulate(triangulateSEXP);
    rcpp_result_gen = Rcpp::wrap(meshVolumeK(rmesh, triangulate));
    return rcpp_result_gen;
END_RCPP
}
// meshAreaK
double meshAreaK(const Rcpp::List rmesh, const bool triangulate);
RcppExport SEXP _MeshesOperations_meshAreaK(SEXP rmeshSEXP, SEXP triangulateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh(rmeshSEXP);
    Rcpp::traits::input_parameter< const bool >::type triangulate(triangulateSEXP);
    rcpp_result_gen = Rcpp::wrap(meshAreaK(rmesh, triangulate));
    return rcpp_result_gen;
END_RCPP
}
// SurfMesh
Rcpp::List SurfMesh(const Rcpp::List rmesh, const bool isTriangle, const bool triangulate, const bool clean, const bool normals, const double epsilon);
RcppExport SEXP _MeshesOperations_SurfMesh(SEXP rmeshSEXP, SEXP isTriangleSEXP, SEXP triangulateSEXP, SEXP cleanSEXP, SEXP normalsSEXP, SEXP epsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh(rmeshSEXP);
    Rcpp::traits::input_parameter< const bool >::type isTriangle(isTriangleSEXP);
    Rcpp::traits::input_parameter< const bool >::type triangulate(triangulateSEXP);
    Rcpp::traits::input_parameter< const bool >::type clean(cleanSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    Rcpp::traits::input_parameter< const double >::type epsilon(epsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(SurfMesh(rmesh, isTriangle, triangulate, clean, normals, epsilon));
    return rcpp_result_gen;
END_RCPP
}
// SurfEMesh
Rcpp::List SurfEMesh(const Rcpp::List rmesh, const bool isTriangle, const bool triangulate, const bool clean, const bool normals, const double epsilon);
RcppExport SEXP _MeshesOperations_SurfEMesh(SEXP rmeshSEXP, SEXP isTriangleSEXP, SEXP triangulateSEXP, SEXP cleanSEXP, SEXP normalsSEXP, SEXP epsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh(rmeshSEXP);
    Rcpp::traits::input_parameter< const bool >::type isTriangle(isTriangleSEXP);
    Rcpp::traits::input_parameter< const bool >::type triangulate(triangulateSEXP);
    Rcpp::traits::input_parameter< const bool >::type clean(cleanSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    Rcpp::traits::input_parameter< const double >::type epsilon(epsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(SurfEMesh(rmesh, isTriangle, triangulate, clean, normals, epsilon));
    return rcpp_result_gen;
END_RCPP
}
// SurfQMesh
Rcpp::List SurfQMesh(const Rcpp::List rmesh, const bool isTriangle, const bool triangulate, const bool clean, const bool normals, const double epsilon);
RcppExport SEXP _MeshesOperations_SurfQMesh(SEXP rmeshSEXP, SEXP isTriangleSEXP, SEXP triangulateSEXP, SEXP cleanSEXP, SEXP normalsSEXP, SEXP epsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh(rmeshSEXP);
    Rcpp::traits::input_parameter< const bool >::type isTriangle(isTriangleSEXP);
    Rcpp::traits::input_parameter< const bool >::type triangulate(triangulateSEXP);
    Rcpp::traits::input_parameter< const bool >::type clean(cleanSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    Rcpp::traits::input_parameter< const double >::type epsilon(epsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(SurfQMesh(rmesh, isTriangle, triangulate, clean, normals, epsilon));
    return rcpp_result_gen;
END_RCPP
}
// Intersection_K
Rcpp::List Intersection_K(const Rcpp::List rmeshes, const bool clean, const bool normals, const Rcpp::LogicalVector triangulate);
RcppExport SEXP _MeshesOperations_Intersection_K(SEXP rmeshesSEXP, SEXP cleanSEXP, SEXP normalsSEXP, SEXP triangulateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmeshes(rmeshesSEXP);
    Rcpp::traits::input_parameter< const bool >::type clean(cleanSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector >::type triangulate(triangulateSEXP);
    rcpp_result_gen = Rcpp::wrap(Intersection_K(rmeshes, clean, normals, triangulate));
    return rcpp_result_gen;
END_RCPP
}
// Intersection_EK
Rcpp::List Intersection_EK(const Rcpp::List rmeshes, const bool clean, const bool normals, const Rcpp::LogicalVector triangulate);
RcppExport SEXP _MeshesOperations_Intersection_EK(SEXP rmeshesSEXP, SEXP cleanSEXP, SEXP normalsSEXP, SEXP triangulateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmeshes(rmeshesSEXP);
    Rcpp::traits::input_parameter< const bool >::type clean(cleanSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector >::type triangulate(triangulateSEXP);
    rcpp_result_gen = Rcpp::wrap(Intersection_EK(rmeshes, clean, normals, triangulate));
    return rcpp_result_gen;
END_RCPP
}
// Intersection_Q
Rcpp::List Intersection_Q(const Rcpp::List rmeshes, const bool clean, const bool normals, const Rcpp::LogicalVector triangulate);
RcppExport SEXP _MeshesOperations_Intersection_Q(SEXP rmeshesSEXP, SEXP cleanSEXP, SEXP normalsSEXP, SEXP triangulateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmeshes(rmeshesSEXP);
    Rcpp::traits::input_parameter< const bool >::type clean(cleanSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector >::type triangulate(triangulateSEXP);
    rcpp_result_gen = Rcpp::wrap(Intersection_Q(rmeshes, clean, normals, triangulate));
    return rcpp_result_gen;
END_RCPP
}
// Difference_K
Rcpp::List Difference_K(const Rcpp::List rmesh1, const Rcpp::List rmesh2, const bool clean, const bool normals, const bool triangulate1, const bool triangulate2);
RcppExport SEXP _MeshesOperations_Difference_K(SEXP rmesh1SEXP, SEXP rmesh2SEXP, SEXP cleanSEXP, SEXP normalsSEXP, SEXP triangulate1SEXP, SEXP triangulate2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh1(rmesh1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh2(rmesh2SEXP);
    Rcpp::traits::input_parameter< const bool >::type clean(cleanSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    Rcpp::traits::input_parameter< const bool >::type triangulate1(triangulate1SEXP);
    Rcpp::traits::input_parameter< const bool >::type triangulate2(triangulate2SEXP);
    rcpp_result_gen = Rcpp::wrap(Difference_K(rmesh1, rmesh2, clean, normals, triangulate1, triangulate2));
    return rcpp_result_gen;
END_RCPP
}
// Difference_EK
Rcpp::List Difference_EK(const Rcpp::List rmesh1, const Rcpp::List rmesh2, const bool clean, const bool normals, const bool triangulate1, const bool triangulate2);
RcppExport SEXP _MeshesOperations_Difference_EK(SEXP rmesh1SEXP, SEXP rmesh2SEXP, SEXP cleanSEXP, SEXP normalsSEXP, SEXP triangulate1SEXP, SEXP triangulate2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh1(rmesh1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh2(rmesh2SEXP);
    Rcpp::traits::input_parameter< const bool >::type clean(cleanSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    Rcpp::traits::input_parameter< const bool >::type triangulate1(triangulate1SEXP);
    Rcpp::traits::input_parameter< const bool >::type triangulate2(triangulate2SEXP);
    rcpp_result_gen = Rcpp::wrap(Difference_EK(rmesh1, rmesh2, clean, normals, triangulate1, triangulate2));
    return rcpp_result_gen;
END_RCPP
}
// Difference_Q
Rcpp::List Difference_Q(const Rcpp::List rmesh1, const Rcpp::List rmesh2, const bool clean, const bool normals, const bool triangulate1, const bool triangulate2);
RcppExport SEXP _MeshesOperations_Difference_Q(SEXP rmesh1SEXP, SEXP rmesh2SEXP, SEXP cleanSEXP, SEXP normalsSEXP, SEXP triangulate1SEXP, SEXP triangulate2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh1(rmesh1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh2(rmesh2SEXP);
    Rcpp::traits::input_parameter< const bool >::type clean(cleanSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    Rcpp::traits::input_parameter< const bool >::type triangulate1(triangulate1SEXP);
    Rcpp::traits::input_parameter< const bool >::type triangulate2(triangulate2SEXP);
    rcpp_result_gen = Rcpp::wrap(Difference_Q(rmesh1, rmesh2, clean, normals, triangulate1, triangulate2));
    return rcpp_result_gen;
END_RCPP
}
// Union_K
Rcpp::List Union_K(const Rcpp::List rmeshes, const bool clean, const bool normals, const Rcpp::LogicalVector triangulate);
RcppExport SEXP _MeshesOperations_Union_K(SEXP rmeshesSEXP, SEXP cleanSEXP, SEXP normalsSEXP, SEXP triangulateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmeshes(rmeshesSEXP);
    Rcpp::traits::input_parameter< const bool >::type clean(cleanSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector >::type triangulate(triangulateSEXP);
    rcpp_result_gen = Rcpp::wrap(Union_K(rmeshes, clean, normals, triangulate));
    return rcpp_result_gen;
END_RCPP
}
// Union_EK
Rcpp::List Union_EK(const Rcpp::List rmeshes, const bool clean, const bool normals, const Rcpp::LogicalVector triangulate);
RcppExport SEXP _MeshesOperations_Union_EK(SEXP rmeshesSEXP, SEXP cleanSEXP, SEXP normalsSEXP, SEXP triangulateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmeshes(rmeshesSEXP);
    Rcpp::traits::input_parameter< const bool >::type clean(cleanSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector >::type triangulate(triangulateSEXP);
    rcpp_result_gen = Rcpp::wrap(Union_EK(rmeshes, clean, normals, triangulate));
    return rcpp_result_gen;
END_RCPP
}
// Union_Q
Rcpp::List Union_Q(const Rcpp::List rmeshes, const bool clean, const bool normals, const Rcpp::LogicalVector triangulate);
RcppExport SEXP _MeshesOperations_Union_Q(SEXP rmeshesSEXP, SEXP cleanSEXP, SEXP normalsSEXP, SEXP triangulateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmeshes(rmeshesSEXP);
    Rcpp::traits::input_parameter< const bool >::type clean(cleanSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector >::type triangulate(triangulateSEXP);
    rcpp_result_gen = Rcpp::wrap(Union_Q(rmeshes, clean, normals, triangulate));
    return rcpp_result_gen;
END_RCPP
}
// isotropicRemeshingK
Rcpp::List isotropicRemeshingK(const Rcpp::List rmesh, const double targetEdgeLength, const unsigned niters, const unsigned nrelaxsteps, const bool triangulate, const bool normals);
RcppExport SEXP _MeshesOperations_isotropicRemeshingK(SEXP rmeshSEXP, SEXP targetEdgeLengthSEXP, SEXP nitersSEXP, SEXP nrelaxstepsSEXP, SEXP triangulateSEXP, SEXP normalsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh(rmeshSEXP);
    Rcpp::traits::input_parameter< const double >::type targetEdgeLength(targetEdgeLengthSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type niters(nitersSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type nrelaxsteps(nrelaxstepsSEXP);
    Rcpp::traits::input_parameter< const bool >::type triangulate(triangulateSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    rcpp_result_gen = Rcpp::wrap(isotropicRemeshingK(rmesh, targetEdgeLength, niters, nrelaxsteps, triangulate, normals));
    return rcpp_result_gen;
END_RCPP
}
// sampleMeshK
Rcpp::NumericMatrix sampleMeshK(const unsigned nsims, const Rcpp::List rmesh, const bool triangulate);
RcppExport SEXP _MeshesOperations_sampleMeshK(SEXP nsimsSEXP, SEXP rmeshSEXP, SEXP triangulateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned >::type nsims(nsimsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh(rmeshSEXP);
    Rcpp::traits::input_parameter< const bool >::type triangulate(triangulateSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleMeshK(nsims, rmesh, triangulate));
    return rcpp_result_gen;
END_RCPP
}
// smoothMeshK
Rcpp::List smoothMeshK(const Rcpp::List rmesh, const double angle, const unsigned niters, const bool triangulate, const bool normals);
RcppExport SEXP _MeshesOperations_smoothMeshK(SEXP rmeshSEXP, SEXP angleSEXP, SEXP nitersSEXP, SEXP triangulateSEXP, SEXP normalsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh(rmeshSEXP);
    Rcpp::traits::input_parameter< const double >::type angle(angleSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type niters(nitersSEXP);
    Rcpp::traits::input_parameter< const bool >::type triangulate(triangulateSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    rcpp_result_gen = Rcpp::wrap(smoothMeshK(rmesh, angle, niters, triangulate, normals));
    return rcpp_result_gen;
END_RCPP
}
// smoothShapeK
Rcpp::List smoothShapeK(const Rcpp::List rmesh, const double time, const unsigned niters, const bool triangulate, const bool normals);
RcppExport SEXP _MeshesOperations_smoothShapeK(SEXP rmeshSEXP, SEXP timeSEXP, SEXP nitersSEXP, SEXP triangulateSEXP, SEXP normalsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh(rmeshSEXP);
    Rcpp::traits::input_parameter< const double >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type niters(nitersSEXP);
    Rcpp::traits::input_parameter< const bool >::type triangulate(triangulateSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    rcpp_result_gen = Rcpp::wrap(smoothShapeK(rmesh, time, niters, triangulate, normals));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MeshesOperations_MinkowskiSumEK", (DL_FUNC) &_MeshesOperations_MinkowskiSumEK, 6},
    {"_MeshesOperations_readFile", (DL_FUNC) &_MeshesOperations_readFile, 1},
    {"_MeshesOperations_writeFile", (DL_FUNC) &_MeshesOperations_writeFile, 5},
    {"_MeshesOperations_clipMeshEK", (DL_FUNC) &_MeshesOperations_clipMeshEK, 6},
    {"_MeshesOperations_connectedComponentsK", (DL_FUNC) &_MeshesOperations_connectedComponentsK, 6},
    {"_MeshesOperations_connectedComponentsEK", (DL_FUNC) &_MeshesOperations_connectedComponentsEK, 6},
    {"_MeshesOperations_connectedComponentsQ", (DL_FUNC) &_MeshesOperations_connectedComponentsQ, 6},
    {"_MeshesOperations_meshVolumeK", (DL_FUNC) &_MeshesOperations_meshVolumeK, 2},
    {"_MeshesOperations_meshAreaK", (DL_FUNC) &_MeshesOperations_meshAreaK, 2},
    {"_MeshesOperations_SurfMesh", (DL_FUNC) &_MeshesOperations_SurfMesh, 6},
    {"_MeshesOperations_SurfEMesh", (DL_FUNC) &_MeshesOperations_SurfEMesh, 6},
    {"_MeshesOperations_SurfQMesh", (DL_FUNC) &_MeshesOperations_SurfQMesh, 6},
    {"_MeshesOperations_Intersection_K", (DL_FUNC) &_MeshesOperations_Intersection_K, 4},
    {"_MeshesOperations_Intersection_EK", (DL_FUNC) &_MeshesOperations_Intersection_EK, 4},
    {"_MeshesOperations_Intersection_Q", (DL_FUNC) &_MeshesOperations_Intersection_Q, 4},
    {"_MeshesOperations_Difference_K", (DL_FUNC) &_MeshesOperations_Difference_K, 6},
    {"_MeshesOperations_Difference_EK", (DL_FUNC) &_MeshesOperations_Difference_EK, 6},
    {"_MeshesOperations_Difference_Q", (DL_FUNC) &_MeshesOperations_Difference_Q, 6},
    {"_MeshesOperations_Union_K", (DL_FUNC) &_MeshesOperations_Union_K, 4},
    {"_MeshesOperations_Union_EK", (DL_FUNC) &_MeshesOperations_Union_EK, 4},
    {"_MeshesOperations_Union_Q", (DL_FUNC) &_MeshesOperations_Union_Q, 4},
    {"_MeshesOperations_isotropicRemeshingK", (DL_FUNC) &_MeshesOperations_isotropicRemeshingK, 6},
    {"_MeshesOperations_sampleMeshK", (DL_FUNC) &_MeshesOperations_sampleMeshK, 3},
    {"_MeshesOperations_smoothMeshK", (DL_FUNC) &_MeshesOperations_smoothMeshK, 5},
    {"_MeshesOperations_smoothShapeK", (DL_FUNC) &_MeshesOperations_smoothShapeK, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_MeshesOperations(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
