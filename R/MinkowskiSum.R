#' @title Minkowski sum of two meshes
#' @description Returns the mesh defined as the Minkowski sum of the 
#'   two input meshes.
#'
#' @param mesh1,mesh2 xxx 
#' @param triangulate Boolean, whether to triangulate the output mesh (note 
#'   that it is not necessarily triangle when the two input meshes are triangle)
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
#' 
#' # example 1: octahedron + sphere 
#' library(Rvcg)
#' mesh1 <- octahedron3d()
#' mesh2 <- vcgSphere(subdivision = 2, normals = FALSE)
#' mesh <- MinkowskiSum(mesh1, mesh2, normals = TRUE)
#' rglmesh <- toRGL(mesh)
#' open3d(windowRect = c(50, 50, 562, 562))
#' view3d(30, 30, zoom = 0.8)
#' shade3d(rglmesh, color = "maroon")
#' 
#' # example2: truncated icosahedron + tetrahedron
#' library(MeshesOperations)
#' library(rgl)
#' # mesh 1
#' mesh1 <- truncatedIcosahedron
#' # mesh 2: regular tetrahedron
#' a <- 1 / sqrt(3)
#' vertices <- rbind(
#' 		c( a, -a, -a),
#' 		c( a,  a,  a),
#' 		c(-a, -a,  a),
#' 		c(-a,  a, -a)
#' )
#' faces <- rbind(
#' 		c(1L, 2L, 3L),
#' 		c(3L, 2L, 4L),
#' 		c(4L, 2L, 1L),
#' 		c(1L, 3L, 4L)
#' )
#' mesh2 <- list(vertices = vertices, faces = faces)
#' # sum 
#' mesh <- MinkowskiSum(mesh1, mesh2, normals = FALSE)
#' # plot
#' rglmesh <- toRGL(mesh)
#' open3d(windowRect = c(50, 50, 562, 562))
#' view3d(30, 30, zoom = 0.8)
#' shade3d(rglmesh, color = "navy")
#' plotEdges(mesh[["vertices"]], mesh[["edges0"]], color = "yellow")
MinkowskiSum <- function(mesh1, mesh2, triangulate = TRUE, normals = FALSE){
  stopifnot(isBoolean(normals))
  vft1 <- getVFT(mesh1)
  vft2 <- getVFT(mesh2)
  mesh <- 
    MinkowskiSumEK(vft1[["rmesh"]], vft2[["rmesh"]], triangulate, normals)
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
