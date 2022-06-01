#' @title Sampling on a mesh
#' @description Uniformly samples points on a mesh.
#'
#' @param n number of simulations, a positive integer
#' @param mesh xxx
#'
#' @return The simulated points on a matrix with three columns.
#' @export
#'
#' @examples
#' library(MeshesOperations)
#' library(rgl)
#' mesh <- torusMesh(R = 4, r = 2)
#' sims <- sampleOnMesh(200, mesh)
#' open3d(windowRect = c(50, 50, 562, 562))
#' view3d(0, 0, zoom = 0.75)
#' shade3d(mesh, color = "yellow", alpha = 0.3)
#' points3d(sims)
sampleOnMesh <- function(n, mesh){
  stopifnot(isStrictPositiveInteger(n))
  vft <- getVFT(mesh)
  sampleMeshK(as.integer(n), vft[["rmesh"]], !vft[["isTriangle"]])
}
