#' @title Clip a mesh
#' @description Clip a mesh to the volume bounded by another mesh.
#'
#' @param mesh a mesh given either as a list containing (at least) the fields 
#'   \code{vertices} and \code{faces}, otherwise a \strong{rgl} mesh 
#'   (i.e. a \code{mesh3d} object)
#' @param clipper a mesh given either as a list containing (at least) the fields 
#'   \code{vertices} and \code{faces}, otherwise a \strong{rgl} mesh 
#'   (i.e. a \code{mesh3d} object); if \code{clipVolume=TRUE}, this mesh will be 
#'   modified: it will be refined with the intersection with \code{mesh}
#' @param clipVolume Boolean, whether the clipping has to be done on the volume 
#'   bounded by \code{mesh} rather than on its surface (i.e. \code{mesh} will be 
#'   kept closed if it is closed)
#' @param normals Boolean, whether to compute the vertex normals of the 
#'   output mesh(es)
#'
#' @return If \code{clipVolume=FALSE}, a triangle mesh represented as the output of the 
#'   \code{\link{Mesh}} function, otherwise a list of two such triangle meshes: the 
#'   clipped mesh in the field \code{"mesh"} and the modified clipping mesh in the 
#'   field \code{"clipper"}.
#' @export
#'
#' @note If \code{clipVolume=TRUE}, the mesh to be clipped (\code{"mesh"}) must be 
#'   without self-intersection.
#' 
#' @examples
#' library(MeshesOperations)
#' library(rgl)
#' library(Rvcg)
#' mesh    <- cube3d()
#' clipper <- scale3d(vcgSphere(), sqrt(2), sqrt(2), sqrt(2))
#' clippedMesh <- clipMesh(mesh, clipper, clipVolume = TRUE)
#' open3d(windowRect = c(50, 50, 950, 500))
#' mfrow3d(1, 2)
#' view3d(zoom = 0.9)
#' shade3d(toRGL(clippedMesh[["mesh"]]), color = "purple")
#' next3d()
#' view3d(zoom = 0.9)
#' shade3d(toRGL(clippedMesh[["clipper"]]), color = "khaki")
clipMesh <- function(mesh, clipper, clipVolume, normals = FALSE){
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
	clip <- clipMeshEK(
		rmesh, rclipper, clipVolume, triangulate1, triangulate2, normals
	)
	if(clipVolume){
		mesh <- clip[["clipper"]]
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
		clipper <- mesh
		mesh    <- clip[["mesh"]]
	}else{
		mesh <- clip
	}
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
	if(clipVolume){
		list("mesh" = mesh, "clipper" = clipper)
	}else{
		mesh
	}
}