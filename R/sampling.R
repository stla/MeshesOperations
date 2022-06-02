#' @title Sampling on a mesh
#' @description Uniformly samples points on a mesh.
#'
#' @param n number of simulations, a positive integer
#' @param mesh either a list containing (at least) two fields \code{vertices} 
#'   (numeric matrix with three columns) and \code{faces} (integer matrix or 
#'   list of integer vectors), otherwise a \strong{rgl} mesh (i.e. a 
#'   \code{mesh3d} object)
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
#' shade3d(mesh, color = "yellow")
#' points3d(sims, size = 5)
sampleOnMesh <- function(n, mesh){
  stopifnot(isStrictPositiveInteger(n))
  vft <- getVFT(mesh)
  sampleMeshK(as.integer(n), vft[["rmesh"]], !vft[["isTriangle"]])
}
