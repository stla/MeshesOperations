#' @title Isotropic remeshing
#' @description Isotropically remesh a mesh.
#'
#' @param vertices a numeric matrix with three columns
#' @param faces either an integer matrix (each row provides the vertex indices
#'   of the corresponding face) or a list of integer vectors, each one
#'   providing the vertex indices of the corresponding face
#' @param mesh if not \code{NULL}, this argument takes precedence over \code{vertices} 
#'   and \code{faces}, and must be either a list containing the fields \code{vertices} 
#'   and \code{faces} (objects as described above), otherwise a \strong{rgl} mesh 
#'   (i.e. a \code{mesh3d} object)
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
#' 
#' theta <- seq(0, 2*pi, length.out = 16)
#' torus <- cylinder3d(
#'   cbind(cos(theta), sin(theta), 0), 
#'   radius = 0.4, closed = TRUE
#' )
#' 
#' mesh <- isotropicRemesh(
#'   mesh = torus,
#'   targetEdgeLength = 0.3,
#'   iterations = 3
#' )
#' rglmesh <- toRGL(mesh)
#' 
#' open3d(windowRect = c(50, 50, 950, 500))
#' mfrow3d(1, 2)
#' view3d(0, 0, zoom = 0.8)
#' wire3d(torus)
#' next3d()
#' view3d(0, 0, zoom = 0.8)
#' wire3d(rglmesh)
isotropicRemesh <- function(
		vertices, faces, mesh = NULL, targetEdgeLength, 
		iterations = 1, relaxSteps = 1, normals = FALSE
){
	stopifnot(isPositiveNumber(targetEdgeLength))
	stopifnot(isStrictPositiveInteger(iterations))
	stopifnot(isStrictPositiveInteger(relaxSteps))
	stopifnot(isBoolean(normals))
	if(!is.null(mesh)){
		if(inherits(mesh, "mesh3d")){
			vft  <- getVFT(mesh, beforeCheck = TRUE)
			mesh <- vft[["rmesh"]]
		}
		vertices <- mesh[["vertices"]]
		faces    <- mesh[["faces"]]
	}
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