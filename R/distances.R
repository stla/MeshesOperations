#' Title
#'
#' @param mesh x
#' @param points x
#'
#' @return x
#' @export
#'
#' @examples
#'#
distancesToMesh <- function(mesh, points){
  if(!is.matrix(points)){
    points <- rbind(points)
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
