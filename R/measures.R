#' @title Mesh volume
#' @description Computes the volume bounded by a mesh.
#'
#' @param mesh xxx
#'
#' @return A number, the volume bounded by the mesh.
#' @export
#'
#' @examples
#' library(MeshesOperations)
#' R <- 4; r <- 2
#' mesh <- torusMesh(R, r)
#' meshVolume(mesh)
#' # true volume of the torus: 
#' 2 * pi^2 * R * r^2
meshVolume <- function(mesh){
  if(inherits(mesh, "mesh3d")){
    triangles <- mesh[["it"]]
    if(!is.null(triangles)){
      triangles <- lapply(1L:ncol(triangles), function(i) triangles[, i] - 1L)
    }
    quads <- mesh[["ib"]]
    isTriangle <- is.null(quads)
    if(!isTriangle){
      quads <- lapply(1L:ncol(quads), function(i) quads[, i] - 1L)
    }
    faces <- c(triangles, quads)
    vertices <- mesh[["vb"]][-4L, ]
    rmesh <- list("vertices" = vertices, "faces" = faces)
  }else if(inherits(mesh, "cgalMesh")){
    isTriangle <- attr(mesh, "toRGL") == 3L
    vertices <- mesh[["vertices"]]
    faces <- mesh[["faces"]]
    if(is.matrix(faces)){
      faces <- lapply(1L:nrow(faces), function(i) faces[i, ] - 1L)
    }
    rmesh <- list("vertices" = vertices, "faces" = faces)
  }else if(is.list(mesh)){
    rmesh <- 
      checkMesh(mesh[["vertices"]], mesh[["faces"]], gmp = FALSE, aslist = TRUE)
    isTriangle <- rmesh[["isTriangle"]]
  }else{
    stop("Invalid `mesh` argument.")
  }
  meshVolumeK(rmesh, !isTriangle)
}
