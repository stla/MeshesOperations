#' @title Smoothing of the shape of a mesh
#' @description Smooths the overall shape of the mesh by using the mean
#'   curvature flow.
#'
#' @param vertices a numeric matrix with three columns
#' @param faces either an integer matrix (each row provides the vertex indices
#'   of the corresponding face) or a list of integer vectors, each one
#'   providing the vertex indices of the corresponding face
#' @param mesh if not \code{NULL}, this argument takes precedence over \code{vertices} 
#'   and \code{faces}, and must be either a list containing two fields \code{vertices} 
#'   and \code{faces} as described above, otherwise a \strong{rgl} mesh (i.e. a 
#'   \code{mesh3d} object)
#' @param time positive number, a time step that corresponds to the speed by
#'   which the surface is smoothed (the larger the faster); typical values lie
#'   between \code{1e-6} and \code{1}
#' @param iterations number of iterations, a positive integer
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
#' # parabola ####
#' x <- seq(-1, 1, length.out = 30)
#' parabola <- cylinder3d(cbind(x, x^2, 0), radius = 0.2, closed = -2)
#' vertices <- t(parabola$vb[-4L, ])
#' faces <- c(
#'   split(t(parabola$it), 1L:ncol(parabola$it)),
#'   split(t(parabola$ib), 1L:ncol(parabola$ib))
#' )
#' sparabola <- smoothShape(
#'   vertices, faces, time = 0.0005, iterations = 10
#' )
#' sparabola <- toRGL(sparabola)
#' open3d(windowRect = c(50, 50, 950, 500))
#' mfrow3d(1, 2)
#' view3d(0, 0, zoom = 0.9)
#' shade3d(parabola, color = "orange")
#' wire3d(parabola)
#' next3d()
#' view3d(0, 0)
#' shade3d(sparabola, color = "green")
#' wire3d(sparabola)
#'
#' # Stanford bunny (light version)
#' vf <- readMeshFile(
#'   system.file("extdata", "bunny.off", package = "MeshesOperations")
#' )
#' mesh <- Mesh(vf[["vertices"]], vf[["faces"]], normals = TRUE)
#' rglmesh <- toRGL(mesh)
#' smesh <- smoothShape(
#'   mesh = mesh,
#'   time = 0.00001, iterations = 1, normals = TRUE
#' )
#' srglmesh <- toRGL(smesh)
#' open3d(windowRect = c(50, 50, 900, 500))
#' mfrow3d(1, 2)
#' view3d(0, 0, zoom = 0.8)
#' shade3d(rglmesh, color = "purple")
#' next3d()
#' view3d(0, 0, zoom = 0.8)
#' shade3d(srglmesh, color = "violetred")
smoothShape <- function(
		vertices, faces, mesh = NULL, time, iterations = 1, normals = FALSE
){
	stopifnot(isPositiveNumber(time))
	stopifnot(isStrictPositiveInteger(iterations))
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
	mesh <- smoothShapeK(
			rmesh, time, as.integer(iterations), triangulate, normals
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
