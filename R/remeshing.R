#' @title Isotropic remeshing
#' @description Isotropically remesh a mesh.
#'
#' @param vertices a numeric matrix with three columns
#' @param faces either an integer matrix (each row provides the vertex indices
#'   of the corresponding face) or a list of integer vectors, each one
#'   providing the vertex indices of the corresponding face
#' @param targetEdgeLength positive number, the target edge length of the 
#'   remeshed mesh
#' @param iterations number of iterations, a positive integer
#' @param relaxSteps number of relaxation steps, a positive integer
#' @param normals Boolean, whether to compute the vertex normals of the 
#'   output mesh
#'
#' @return A triangle mesh represented as the output of the 
#'   \code{\link{Mesh}} function.
#' @export
#'
#' @examples
#' library(MeshesOperations)
#' library(rgl)
isotropicRemesh <- function(
  vertices, faces, targetEdgeLength, 
  iterations = 1, relaxSteps = 1, normals = FALSE
){
  stopifnot(isPositiveNumber(targetEdgeLength))
  stopifnot(isStrictPositiveInteger(iterations))
  stopifnot(isStrictPositiveInteger(relaxSteps))
  stopifnot(isBoolean(normals))
  checkedMesh <- checkMesh(vertices, faces, gmp = FALSE, aslist = TRUE)
  vertices         <- checkedMesh[["vertices"]]
  faces            <- checkedMesh[["faces"]]
  isTriangle       <- checkedMesh[["isTriangle"]]
  rmesh <- list("vertices" = vertices, "faces" = faces)
  triangulate <- !isTriangle
  mesh <- isotropicRemeshingK(
    rmesh, targetEdgeLength, as.integer(iterations), as.integer(relaxSteps), 
    triangulate, normals
  )
  mesh[["vertices"]] <- t(mesh[["vertices"]])
  mesh[["faces"]] <- t(mesh[["faces"]])
  edges <- unname(t(mesh[["edges"]]))
  exteriorEdges <- edges[edges[, 3L] == 1L, c(1L, 2L)]
  mesh[["exteriorEdges"]] <- exteriorEdges
  mesh[["exteriorVertices"]] <- which(table(exteriorEdges) != 2L)
  mesh[["edges"]] <- edges[, c(1L, 2L)]
  if(normals){
    mesh[["normals"]] <- t(mesh[["normals"]])
  }
  attr(mesh, "toRGL") <- 3L
  class(mesh) <- "cgalMesh"
  mesh
}