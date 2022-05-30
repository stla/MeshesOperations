#' @title Read a mesh file
#' @description Read mesh vertices and faces from a file.
#'
#' @param filepath path to the mesh file; supported formats are \code{stl},
#'   \code{ply}, \code{obj} and \code{off}
#'
#' @return A list with two fields: \code{vertices}, a numeric matrix with three 
#'   columns, and \code{faces}, either a list of integer vectors or, in the 
#'   case if all faces have the same number of sides, an integer matrix.
#' @export
#' 
#' @importFrom data.table uniqueN
#'
#' @examples
#' library(MeshesOperations)
#' library(rgl)
#' vf <- readMeshFile(
#'   system.file("extdata", "beethoven.ply", package = "MeshesOperations")
#' )
#' mesh <- Mesh(vf[["vertices"]], vf[["faces"]], normals = TRUE)
#' rglmesh <- toRGL(mesh)
#' open3d(windowRect = c(50, 50, 562, 562))
#' view3d(0, 0, zoom = 0.8)
#' shade3d(rglmesh, color = "palevioletred")
readMeshFile <- function(filepath){
  if(!file.exists(filepath)){
    stop("File not found.")
  }
  mesh <- readFile(filepath)
  faces <- mesh[["faces"]]
  usizes <- uniqueN(lengths(faces))
  if(usizes == 1L){
    mesh[["faces"]] <- do.call(rbind, faces)
  }
  mesh
}