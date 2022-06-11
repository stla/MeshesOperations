#' @title Distance to a mesh
#'
#' @description Computes the distances from given points to a mesh.
#'
#' @param mesh a mesh given either as a list containing (at least) the fields
#'   \code{vertices} and \code{faces}, otherwise a \strong{rgl} mesh
#'   (i.e. a \code{mesh3d} object)
#' @param points either one point given as a numeric vector or several points
#'   given as a numeric matrix with three columns
#'
#' @return A numeric vector providing the distances between the given point(s)
#'   to the mesh.
#' @export
#'
#' @examples
#' # cube example ####
#' library(MeshesOperations)
#' mesh <- rgl::cube3d()
#' points <- rbind(
#'   c(0, 0, 0),
#'   c(1, 1, 1)
#' )
#' distancesToMesh(mesh, points) # should be 1 and 0
#'
#' # cyclide example ####
#' library(MeshesOperations)
#' a <- 100; c <- 30; mu <- 80
#' mesh <- cyclideMesh(a, c, mu, nu = 100L, nv = 100L)
#' points <- c(c, 0, 0)
#' # should be a - mu = 20 (see ?cyclideMesh):
#' distancesToMesh(mesh, points)
distancesToMesh <- function(mesh, points){
  if(!is.matrix(points)){
    points <- rbind(points)
  }
  if(ncol(points) != 3L){
    stop(
      "The `points` argument must be a vector of length 3 or ",
      "a matrix with three columns."
    )
  }
  stopifnot(is.numeric(points))
  storage.mode(points) <- "double"
  if(anyNA(points)){
    stop("Found missing values in `points`.")
  }
  tpoints <- t(points)
  if(inherits(mesh, "mesh3d")){
    vft  <- getVFT(mesh, beforeCheck = TRUE)
    mesh <- vft[["rmesh"]]
  }
  vertices <- mesh[["vertices"]]
  faces    <- mesh[["faces"]]
  checkedMesh <- checkMesh(vertices, faces, gmp = FALSE, aslist = TRUE)
  vertices         <- checkedMesh[["vertices"]]
  faces            <- checkedMesh[["faces"]]
  isTriangle       <- checkedMesh[["isTriangle"]]
  triangulate <- !isTriangle
  rmesh <- list("vertices" = vertices, "faces" = faces)
  distanceK(rmesh, tpoints, triangulate)
}
