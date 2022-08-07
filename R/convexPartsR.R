#' @title Decomposition into convex parts
#' @description Decomposition of a mesh into convex parts.
#'
#' @param vertices a numeric matrix with three columns, or a \code{bigq}
#'   matrix with three columns if \code{numbersType="gmp"}
#' @param faces either an integer matrix (each row provides the vertex indices
#'   of the corresponding face) or a list of integer vectors, each one
#'   providing the vertex indices of the corresponding face
#' @param mesh if not \code{NULL}, this argument takes precedence over
#'   \code{vertices} and \code{faces}, and must be either a list containing
#'   the fields \code{vertices} and \code{faces} (objects as described above),
#'   otherwise a \strong{rgl} mesh (i.e. a \code{mesh3d} object)
#' @param triangulate Boolean, whether to triangulate the convex parts
#'
#' @return A list of \code{cgalMesh} lists, each corresponding to a convex part.
#' @export
#'
#' @importFrom data.table uniqueN
#'
#' @examples
#' # a non-convex polyhedron ####
#' \donttest{library(MeshesOperations)
#' library(rgl)
#' library(randomcoloR)
#' meshes <- convexParts(mesh = NonConvexPolyhedron)
#' ncp <- length(meshes)
#' colors <- randomColor(ncp, hue = "random", luminosity = "bright")
#' open3d(windowRect = c(50, 50, 562, 562), zoom = 0.8)
#' for(i in seq_len(ncp)){
#'   shade3d(toRGL(meshes[[i]]), color = colors[i])
#' }
#' plotEdges(
#'   NonConvexPolyhedron[["vertices"]],
#'   NonConvexPolyhedron[["edges"]]
#' )}
#'
#' # pentagrammic prism ####
#' \donttest{library(MeshesOperations)
#' library(rgl)
#' library(randomcoloR)
#' meshes <- convexParts(mesh = pentagrammicPrism)
#' ncp <- length(meshes)
#' colors <- randomColor(ncp, hue = "random", luminosity = "bright")
#' open3d(windowRect = c(50, 50, 562, 562), zoom = 0.8)
#' for(i in seq_len(ncp)){
#'   shade3d(toRGL(meshes[[i]]), color = colors[i])
#' }
#' plotEdges(
#'   pentagrammicPrism[["vertices"]],
#'   pentagrammicPrism[["edges"]],
#'   tubesRadius = 0.01,
#'   spheresRadius = 0.02
#' )}
convexParts <- function(
  vertices, faces, mesh = NULL, triangulate = TRUE
){
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
  rmesh <- list("vertices" = vertices, "faces" = faces)
  cxparts <- convexDecomposition(rmesh, triangulate)
  ncp <- length(cxparts)
  meshes <- vector("list", ncp)
  for(i in seq_len(ncp)){
    mesh <- cxparts[[i]]
    mesh[["vertices"]] <- t(mesh[["vertices"]])
    edgesDF <- mesh[["edges"]]
    mesh[["edgesDF"]] <- edgesDF
    mesh[["edges"]] <- as.matrix(edgesDF[, c("i1", "i2")])
    exteriorEdges <- as.matrix(subset(edgesDF, exterior)[, c("i1", "i2")])
    mesh[["exteriorEdges"]] <- exteriorEdges
    mesh[["exteriorVertices"]] <- which(table(exteriorEdges) != 2L)
    # if(triangulate){
    #   edges0DF <- mesh[["edges0"]]
    #   mesh[["edges0DF"]] <- edges0DF
    #   mesh[["edges0"]] <- as.matrix(edges0DF[, c("i1", "i2")])
    #   if(normals){
    #     mesh[["normals0"]] <- t(mesh[["normals0"]])
    #   }
    # }
    if(triangulate){
      mesh[["faces"]] <- do.call(rbind, mesh[["faces"]])
      toRGL <- 3L
    }else{
      sizes <- lengths(mesh[["faces"]])
      usizes <- uniqueN(sizes)
      if(usizes == 1L){
        if(sizes[1L] %in% c(3L, 4L)){
          toRGL <- sizes[1L]
        }
        mesh[["faces"]] <- do.call(rbind, mesh[["faces"]])
      }else if(usizes == 2L && all(sizes %in% c(3L, 4L))){
        toRGL <- 34L
      }else{
        toRGL <- FALSE
      }
    }
    attr(mesh, "toRGL") <- toRGL
    class(mesh) <- "cgalMesh"
    meshes[[i]] <- mesh
  }
  meshes
}
