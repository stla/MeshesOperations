#ifndef _MESHESOPERATIONSHEADER_
#include "MeshesOperations.h"
#endif

template <typename MeshT>
void checkMesh(MeshT mesh, size_t i) {
  const bool si = PMP::does_self_intersect(mesh);
  if(si) {
    std::string msg = "Mesh n\u00b0" + std::to_string(i) + " self-intersects.";
    Rcpp::stop(msg);
  }
  const bool bv = PMP::does_bound_a_volume(mesh);
  if(!bv) {
    std::string msg =
        "Mesh n\u00b0" + std::to_string(i) + " does not bound a volume.";
    Rcpp::stop(msg);
  }
}

template <typename KernelT, typename MeshT, typename PointT>
MeshT Intersection(const Rcpp::List rmeshes,
                   const bool clean,
                   const bool exact,  // ?
                   const Rcpp::LogicalVector triangulate) {
  const size_t nmeshes = rmeshes.size();
  std::vector<MeshT> meshes(nmeshes);
  Rcpp::List rmesh = Rcpp::as<Rcpp::List>(rmeshes(0));
  Message("Processing mesh n\u00b01...\n");
  MeshT mesh_0 = makeSurfTMesh<MeshT, PointT>(rmesh, clean);
  if(triangulate[0]) {
    const bool success = PMP::triangulate_faces(mesh_0);
    if(!success) {
      Rcpp::stop("Triangulation of mesh n\u00b01 has failed.");
    }
  }
  if(true) {
    checkMesh<MeshT>(mesh_0, 1);
  }
  meshes[0] = mesh_0;
  for(size_t i = 1; i < nmeshes; i++) {
    // if(true) {
    //   checkMesh<MeshT>(meshes[i - 1], i);
    // }
    const std::string meshnum = std::to_string(i + 1);
    Rcpp::List rmesh_i = Rcpp::as<Rcpp::List>(rmeshes(i));
    Message("Processing mesh n\u00b0" + meshnum + "...\n");
    MeshT mesh_i = makeSurfTMesh<MeshT, PointT>(rmesh_i, clean);
    if(triangulate[i]) {
      const bool success = PMP::triangulate_faces(mesh_i);
      if(!success) {
        Rcpp::stop("Triangulation of mesh n\u00b0" + meshnum + " has failed.");
      }
    }
    checkMesh<MeshT>(mesh_i, i + 1);
    bool ok = PMP::corefine_and_compute_intersection(meshes[i - 1], mesh_i,
                                                     meshes[i]);
    if(!ok) {
      Rcpp::stop("Intersection computation has failed.");
    }
  }
  return meshes[nmeshes - 1];
}

// [[Rcpp::export]]
Rcpp::List Intersection_K(const Rcpp::List rmeshes,
                          const bool clean,
                          const bool normals,
                          const Rcpp::LogicalVector triangulate) {
  Mesh3 mesh =
      Intersection<K, Mesh3, Point3>(rmeshes, clean, false, triangulate);
  return RSurfTKMesh(mesh, normals, 0);
}

// [[Rcpp::export]]
Rcpp::List Intersection_EK(const Rcpp::List rmeshes,
                           const bool clean,
                           const bool normals,
                           const Rcpp::LogicalVector triangulate) {
  EMesh3 mesh =
      Intersection<EK, EMesh3, EPoint3>(rmeshes, clean, true, triangulate);
  return RSurfTEKMesh(mesh, normals, 0);
}

// [[Rcpp::export]]
Rcpp::List Intersection_Q(const Rcpp::List rmeshes,  // must be triangles
                          const bool clean,
                          const bool normals,
                          const Rcpp::LogicalVector triangulate) {
  const size_t nmeshes = rmeshes.size();
  std::vector<QMesh3> meshes(nmeshes);
  Rcpp::List rmesh = Rcpp::as<Rcpp::List>(rmeshes(0));
  Message("Processing mesh n\u00b01...\n");
  QMesh3 mesh_0 = makeSurfTQMesh(rmesh, clean);
  if(triangulate[0]) {
    const bool success = PMP::triangulate_faces(mesh_0);
    if(!success) {
      Rcpp::stop("Triangulation of mesh n\u00b01 has failed.");
    }
  }
  checkMesh<QMesh3>(mesh_0, 0);
  meshes[0] = mesh_0;
  for(size_t i = 1; i < nmeshes; i++) {
    const std::string meshnum = std::to_string(i + 1);
    Rcpp::List rmesh_i = Rcpp::as<Rcpp::List>(rmeshes(i));
    Message("Processing mesh n\u00b0" + meshnum + "...\n");
    QMesh3 mesh_i = makeSurfTQMesh(rmesh_i, clean);
    if(triangulate[i]) {
      const bool success = PMP::triangulate_faces(mesh_i);
      if(!success) {
        Rcpp::stop("Triangulation of mesh n\u00b0" + meshnum + " has failed.");
      }
    }
    checkMesh<QMesh3>(mesh_i, i);
    bool ok = PMP::corefine_and_compute_intersection(meshes[i - 1], mesh_i,
                                                     meshes[i]);
    if(!ok) {
      Rcpp::stop("Intersection computation has failed.");
    }
  }
  return RSurfTQMesh(meshes[nmeshes - 1], normals, 0);
}

template <typename KernelT, typename MeshT, typename PointT>
MeshT Difference(const Rcpp::List rmesh1,  // must be triangles
                 const Rcpp::List rmesh2,
                 const bool clean) {
  Message("Processing mesh n\u00b01.\n");
  MeshT smesh1 = makeSurfTMesh<MeshT, PointT>(rmesh1, clean);
  checkMesh<MeshT>(smesh1, 1);
  Message("Processing mesh n\u00b02.\n");
  MeshT smesh2 = makeSurfTMesh<MeshT, PointT>(rmesh2, clean);
  checkMesh<MeshT>(smesh2, 2);
  MeshT outmesh;
  bool ok = PMP::corefine_and_compute_difference(smesh1, smesh2, outmesh);
  if(!ok) {
    Rcpp::stop("Difference computation has failed.");
  }
  return outmesh;
}

// [[Rcpp::export]]
Rcpp::List Difference_K(const Rcpp::List rmesh1,
                        const Rcpp::List rmesh2,
                        const bool clean,
                        const bool normals) {
  Mesh3 mesh = Difference<K, Mesh3, Point3>(rmesh1, rmesh2, clean);
  return RSurfTKMesh(mesh, normals, 0);
}

// [[Rcpp::export]]
Rcpp::List Difference_EK(const Rcpp::List rmesh1,
                         const Rcpp::List rmesh2,
                         const bool clean,
                         const bool normals) {
  EMesh3 mesh = Difference<EK, EMesh3, EPoint3>(rmesh1, rmesh2, clean);
  return RSurfTEKMesh(mesh, normals, 0);
}

// [[Rcpp::export]]
Rcpp::List Difference_Q(const Rcpp::List rmesh1,  // must be triangles
                        const Rcpp::List rmesh2,
                        const bool clean,
                        const bool normals) {
  Message("Processing mesh n\u00b01.\n");
  QMesh3 smesh1 = makeSurfTQMesh(rmesh1, clean);
  checkMesh<QMesh3>(smesh1, 1);
  Message("Processing mesh n\u00b02.\n");
  QMesh3 smesh2 = makeSurfTQMesh(rmesh2, clean);
  checkMesh<QMesh3>(smesh2, 2);
  QMesh3 outmesh;
  bool ok = PMP::corefine_and_compute_difference(smesh1, smesh2, outmesh);
  if(!ok) {
    Rcpp::stop("Difference computation has failed.");
  }
  return RSurfQMesh(outmesh, normals, 0);
}

template <typename KernelT, typename MeshT, typename PointT>
MeshT Union(const Rcpp::List rmeshes,  // must be triangles
            const bool clean,
            const bool exact) {
  const size_t nmeshes = rmeshes.size();
  std::vector<MeshT> meshes(nmeshes);
  Rcpp::List rmesh = Rcpp::as<Rcpp::List>(rmeshes(0));
  Message("Processing mesh n\u00b01.\n");
  meshes[0] = makeSurfTMesh<MeshT, PointT>(rmesh, clean);
  if(exact) {
    checkMesh<MeshT>(meshes[0], 1);
  }
  for(size_t i = 1; i < nmeshes; i++) {
    if(!exact) {
      checkMesh<MeshT>(meshes[i - 1], i);
    }
    Rcpp::List rmesh_i = Rcpp::as<Rcpp::List>(rmeshes(i));
    Message("Processing mesh n\u00b0" + std::to_string(i + 1) + ".\n");
    MeshT mesh_i = makeSurfTMesh<MeshT, PointT>(rmesh_i, clean);
    checkMesh<MeshT>(mesh_i, i + 1);
    bool ok = PMP::corefine_and_compute_union(meshes[i - 1], mesh_i, meshes[i]);
    if(!ok) {
      Rcpp::stop("Union computation has failed.");
    }
  }
  return meshes[nmeshes - 1];
}

// [[Rcpp::export]]
Rcpp::List Union_K(const Rcpp::List rmeshes,
                   const bool clean,
                   const bool normals) {
  Mesh3 mesh = Union<K, Mesh3, Point3>(rmeshes, clean, false);
  return RSurfTKMesh(mesh, normals, 0);
}

// [[Rcpp::export]]
Rcpp::List Union_EK(const Rcpp::List rmeshes,
                    const bool clean,
                    const bool normals) {
  EMesh3 mesh = Union<EK, EMesh3, EPoint3>(rmeshes, clean, true);
  return RSurfTEKMesh(mesh, normals, 0);
}

// [[Rcpp::export]]
Rcpp::List Union_Q(const Rcpp::List rmeshes,  // must be triangles
                   const bool clean,
                   const bool normals) {
  const size_t nmeshes = rmeshes.size();
  std::vector<QMesh3> meshes(nmeshes);
  Rcpp::List rmesh = Rcpp::as<Rcpp::List>(rmeshes(0));
  Message("Processing mesh n\u00b01.\n");
  meshes[0] = makeSurfTQMesh(rmesh, clean);
  checkMesh<QMesh3>(meshes[0], 0);
  for(size_t i = 1; i < nmeshes; i++) {
    Rcpp::List rmesh_i = Rcpp::as<Rcpp::List>(rmeshes(i));
    Message("Processing mesh n\u00b0" + std::to_string(i + 1) + ".\n");
    QMesh3 mesh_i = makeSurfTQMesh(rmesh_i, clean);
    checkMesh<QMesh3>(mesh_i, i);
    bool ok = PMP::corefine_and_compute_union(meshes[i - 1], mesh_i, meshes[i]);
    if(!ok) {
      Rcpp::stop("Union computation has failed.");
    }
  }
  return RSurfTQMesh(meshes[nmeshes - 1], normals, 0);
}
