#' @title Minkowski sum of two meshes
#' @description Returns the mesh defined as the Minkowski sum of the 
#'   two input meshes.
#'
#' @param mesh1,mesh2 xxx 
#' @param normals Boolean, whether to compute the vertex normals of the 
#'   output mesh 
#'
#' @return xxx
#' @export
#'
#' @importFrom data.table uniqueN
#'
#' @examples
#' library(MeshesOperations)
#' library(rgl)
#' library(Rvcg)
#' # Stanford bunny (light version)
#' mesh1 <- readMeshFile(
#'   system.file("extdata", "bunny.off", package = "MeshesOperations")
#' )
#' mesh2 <- vcgSphere(subdivision = 2, normals = FALSE)
#' mesh <- MinkowskiSum(mesh1, mesh2)
#' rglmesh <- toRGL(mesh)
#' open3d(windowRect = c(50, 50, 900, 500))
#' view3d(0, 0, zoom = 0.8)
#' shade3d(rglmesh, color = "maroon")
MinkowskiSum <- function(mesh1, mesh2, normals = FALSE){
  stopifnot(isBoolean(normals))
  vft1 <- getVFT(mesh1)
  vft2 <- getVFT(mesh2)
  mesh <- MinkowskiSumEK(vft1[["rmesh"]], vft2[["rmesh"]], normals)
  mesh[["vertices"]] <- t(mesh[["vertices"]])
  toRGL <- FALSE
  faces <- mesh[["faces"]]
  sizes <- lengths(faces)
  usizes <- uniqueN(sizes)
  if(usizes == 1L){
    if(sizes[1L] %in% c(3L, 4L)){
      toRGL <- sizes[1L]
    }
    mesh[["faces"]] <- do.call(rbind, faces)
  }else if(usizes == 2L && all(sizes %in% c(3L, 4L))){
    toRGL <- 34L
  }
  edges <- unname(t(mesh[["edges"]]))
  exteriorEdges <- edges[edges[, 3L] == 1L, c(1L, 2L)]
  mesh[["exteriorEdges"]] <- exteriorEdges
  mesh[["exteriorVertices"]] <- which(table(exteriorEdges) != 2L)
  mesh[["edges"]] <- edges[, c(1L, 2L)]
  if(normals){
    mesh[["normals"]] <- t(mesh[["normals"]])
  }
  attr(mesh, "toRGL") <- toRGL
  class(mesh) <- "cgalMesh"
  mesh
}
