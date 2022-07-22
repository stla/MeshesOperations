#' @title Clip a mesh
#' @description Clip a mesh to the volume bounded by another mesh.
#'
#' @param mesh a mesh given either as a list containing (at least) the fields
#'   \code{vertices} and \code{faces}, otherwise a \strong{rgl} mesh
#'   (i.e. a \code{mesh3d} object)
#' @param clipper a mesh given either as a list containing (at least) the fields
#'   \code{vertices} and \code{faces}, otherwise a \strong{rgl} mesh
#'   (i.e. a \code{mesh3d} object)
#' @param clipVolume Boolean, whether the clipping has to be done on the volume
#'   bounded by \code{mesh} rather than on its surface (i.e. \code{mesh} will be
#'   kept closed if it is closed)
#' @param normals Boolean, whether to compute the vertex normals of the
#'   output mesh
#'
#' @return A triangle mesh represented as the output of the
#'   \code{\link{Mesh}} function.
#'
#' @export
#'
#' @note If \code{clipVolume=TRUE}, the mesh to be clipped (\code{mesh})
#'   must be without self-intersection.
#'
#' @examples
#' # cube clipped to sphere
#' library(MeshesOperations)
#' library(rgl)
#' mesh    <- cube3d()
#' clipper <- sphereMesh(r= sqrt(2))
#' clippedMesh <- clipMesh(mesh, clipper)
#' open3d(windowRect = c(50, 50, 562, 562))
#' view3d(zoom = 0.9)
#' shade3d(toRGL(clippedMesh), color = "purple")
#'
#' # Barth sextic ####
#' \donttest{library(MeshesOperations)
#' library(rgl)
#' library(rmarchingcubes)
#' # isosurface function
#' gold <- (1+sqrt(5))/2
#' f <- function(x,y,z){
#' 	x2 <- x*x; y2 <- y*y; z2 <- z*z
#' 	4*(gold^2*x2-y2)*(gold^2*y2-z2)*(gold^2*z2-x2) -
#' 			(1+2*gold)*(x2+y2+z2-1)^2
#' }
#' # grid
#' n <- 200L
#' x <- y <- z <- seq(-sqrt(3), sqrt(3), length.out = n)
#' g <- expand.grid(X = x, Y = y, Z = z)
#' # calculate voxel
#' voxel <- array(with(g, f(X, Y, Z)), dim = c(n, n, n))
#' # calculate isosurface
#' contour_shape <- contour3d(
#' 		griddata = voxel, level = 0, x = x, y = y, z = z
#' )
#' # make rgl mesh (plotted later)
#' mesh <- tmesh3d(
#' 		vertices = t(contour_shape[["vertices"]]),
#' 		indices  = t(contour_shape[["triangles"]]),
#' 		normals  = contour_shape[["normals"]],
#' 		homogeneous = FALSE
#' )
#' # clip to sphere of radius sqrt(3)
#' clipper <- sphereMesh(r = sqrt(3))
#' clippedMesh <- clipMesh(mesh, clipper, clipVolume = FALSE, normals = TRUE)
#' # plot
#' open3d(windowRect = c(50, 50, 950, 500))
#' mfrow3d(1, 2)
#' view3d(zoom = 0.8)
#' shade3d(mesh, color = "darkred")
#' next3d()
#' view3d(zoom = 0.8)
#' shade3d(toRGL(clippedMesh), color = "darkred")}
clipMesh <- function(mesh, clipper, clipVolume = TRUE, normals = FALSE){
	stopifnot(isBoolean(clipVolume))
	stopifnot(isBoolean(normals))
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
	rmesh <- list("vertices" = vertices, "faces" = faces)
	triangulate1 <- !isTriangle
	if(inherits(clipper, "mesh3d")){
		vft  <- getVFT(clipper, beforeCheck = TRUE)
		clipper <- vft[["rmesh"]]
	}
	vertices <- clipper[["vertices"]]
	faces    <- clipper[["faces"]]
	checkedMesh <- checkMesh(vertices, faces, gmp = FALSE, aslist = TRUE)
	vertices         <- checkedMesh[["vertices"]]
	faces            <- checkedMesh[["faces"]]
	isTriangle       <- checkedMesh[["isTriangle"]]
	rclipper <- list("vertices" = vertices, "faces" = faces)
	triangulate2 <- !isTriangle
	mesh <- clipMeshEK(
		rmesh, rclipper, clipVolume, triangulate1, triangulate2, normals
	)
	# if(clipVolume){
	# 	mesh <- clip[["clipper"]]
	# 	mesh[["vertices"]] <- t(mesh[["vertices"]])
	# 	mesh[["faces"]] <- t(mesh[["faces"]])
	# 	edges <- unname(t(mesh[["edges"]]))
	# 	exteriorEdges <- edges[edges[, 3L] == 1L, c(1L, 2L)]
	# 	mesh[["exteriorEdges"]] <- exteriorEdges
	# 	mesh[["exteriorVertices"]] <- which(table(exteriorEdges) != 2L)
	# 	mesh[["edges"]] <- edges[, c(1L, 2L)]
	# 	if(normals){
	# 		mesh[["normals"]] <- t(mesh[["normals"]])
	# 	}
	# 	attr(mesh, "toRGL") <- 3L
	# 	class(mesh) <- "cgalMesh"
	# 	clipper <- mesh
	# 	mesh    <- clip[["mesh"]]
	# }else{
	# 	mesh <- clip
	# }
	mesh[["vertices"]] <- t(mesh[["vertices"]])
	mesh[["faces"]] <- t(mesh[["faces"]])
	edgesDF <- mesh[["edges"]]
	mesh[["edgesDF"]] <- edgesDF
	mesh[["edges"]] <- as.matrix(edgesDF[, c("i1", "i2")])
	exteriorEdges <- as.matrix(subset(edgesDF, exterior)[, c("i1", "i2")])
	mesh[["exteriorEdges"]] <- exteriorEdges
	mesh[["exteriorVertices"]] <- which(table(exteriorEdges) != 2L)
	if(normals){
		mesh[["normals"]] <- t(mesh[["normals"]])
	}
	attr(mesh, "toRGL") <- 3L
	class(mesh) <- "cgalMesh"
	# if(clipVolume){
	# 	list("mesh" = mesh, "clipper" = clipper)
	# }else{
	# 	mesh
	# }
	mesh
}
